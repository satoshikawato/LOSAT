//! EXPERIMENT (LOSAT_X_BXSUBJFLOW[=k], in the strict set: on unless `=0`):
//! the subject loop of a query chunk as a flow of tasks of `k` consecutive subjects
//! (`search::X_BX_CHUNK` = 64, the waves' subjects per worker, when the variable is `1` or not a
//! number from 2 to 65,536), searched ahead by the other threads of the pool and handed to the sink
//! in subject order by the thread that runs the loop. Shadow: `LOSAT_X_BXSUBJFLOWSHADOW`. With
//! `LOSAT_X_BXTASKRES` (in the strict set: on unless `=0`) a task also ends before the subject that
//! would take its residues past `k` times the mean subject length (`XTasks`).
//!
//! The waves of `search::search_core` search `threads x X_BX_CHUNK` (= 64) subjects, wait for
//! the slowest worker of the wave, then hand the wave's subjects to the sink on one thread while
//! the other threads wait for the next wave. Here the helpers of `utils::xahead` search the tasks
//! up to a window ahead of the owner (the thread of the loop) without stopping between waves; the
//! owner takes each task's subjects in order, or searches the task itself when no helper has
//! started it, and makes the same sink calls, in the same order, as the waves.
//!
//! A task's value is `compute_subject` for each of its subjects, which depends on the subject and
//! the query chunk alone, not on the diagonal table it runs on, while no diagonal offset passes
//! 2^30: an entry left by an earlier subject lies at least a window behind the next subject's
//! offsets (`Diagonals::finish_subject`), as in a cleared table. `chosen` keeps the waves when
//! `-window_size` plus twice the longest subject exceeds 2^28, where the bit-field wrap of NCBI's
//! `DiagStruct` could make the table's history matter; otherwise the offset, reset at
//! INT4_MAX / 4, stays below 2^29 + 2^28 < 2^30.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ```
//! NCBI scans the subjects one after another in this loop (with `-num_threads 1`). The flow keeps
//! the loop's order for everything after `compute_subject` (the sink: linking, reaping, the HSP
//! stream) and changes only which thread computes a subject and when.

use super::{compute_subject, SubjectComputation};
use crate::algorithm::blastx::{
    args::ResolvedOptions,
    parameters::ContextParameters,
    query_setup::PreparedQueryBatch,
    search::PreliminarySink,
    seed::{Diagonals, Lookup},
    subjects::Subjects,
};
use anyhow::Result;
use std::sync::OnceLock;

/// No NCBI counterpart: the subjects of one task, the oid of the first subject `compute_subject`
/// refused and its error. A task stops at its first error, so the owner hands the subjects before
/// it to the sink and then returns the error, as the waves do.
type Task = (
    Vec<(usize, SubjectComputation)>,
    Option<(usize, anyhow::Error)>,
);

/// EXPERIMENT (LOSAT_X_BXSUBJFLOW[=k]): subjects per task; `None` when the switch is off.
pub(super) fn subjects_per_task() -> Option<usize> {
    static K: OnceLock<Option<usize>> = OnceLock::new();
    *K.get_or_init(|| {
        crate::utils::xswitch::value("LOSAT_X_BXSUBJFLOW", 2..=65_536, super::X_BX_CHUNK)
    })
}

/// EXPERIMENT (LOSAT_X_BXTASKRES, in the strict set: on unless `=0`): whether the flow's tasks are
/// cut by residues as well as by subjects (`XTasks`).
fn task_residues() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXTASKRES"))
}

/// No NCBI counterpart: the subjects of each task of a flow (scheduling only). `Grid`: task `t` is
/// the subjects `t * k..(t + 1) * k`. `Cut` (`LOSAT_X_BXTASKRES`): the first subject of each task,
/// then the number of subjects; a task ends after `k` subjects or before the subject that would
/// take its residues past `k` times the mean subject length, so a run of long subjects (the
/// isoforms of one long protein, stored next to each other) is spread over several tasks instead
/// of making one task many times the others, which the ordered flow then waits for at its end.
/// With `-db` read in blocks the bounds for each `k` are made once per run and shared
/// (`XStreamDb::x_task_starts`): a database of tens of millions of sequences would otherwise read
/// its index and hold a copy of the bounds for every query batch.
pub(super) enum XTasks {
    Grid { k: usize, subjects: usize },
    Cut(std::sync::Arc<[u32]>),
}

impl XTasks {
    /// The tasks of `k` subjects (at most) of `subjects`. `Grid` without the switch, and where the
    /// lengths cannot be read (an index read error: the task that reads the subject then meets the
    /// error where it did before) or the subjects do not fit `u32`.
    pub(super) fn new(k: usize, subjects: Subjects<'_>) -> Self {
        let n = subjects.len();
        let grid = XTasks::Grid { k, subjects: n };
        if !task_residues() {
            return grid;
        }
        let tasks = Self::cut(k, subjects);
        crate::utils::xreach::x_reach!("LOSAT_X_BXTASKRES", tasks.is_some());
        tasks.unwrap_or(grid)
    }

    /// `Cut`, or `None` where `new` keeps the grid.
    fn cut(k: usize, subjects: Subjects<'_>) -> Option<Self> {
        match subjects {
            Subjects::Stream(stream) => stream.db.x_task_starts(k, || Self::starts(k, subjects)),
            _ => Self::starts(k, subjects),
        }
        .map(XTasks::Cut)
    }

    /// `Cut`'s bounds for tasks of `k` subjects of `subjects`, or `None` where `new` keeps the grid.
    fn starts(k: usize, subjects: Subjects<'_>) -> Option<std::sync::Arc<[u32]>> {
        let n = subjects.len();
        if n == 0 || u32::try_from(n).is_err() {
            return None;
        }
        let total = match subjects.total_length() {
            Ok(total) if total > 0 => total as u128,
            _ => return None,
        };
        // `k` times the mean length, rounded up.
        let budget = (k as u128 * total).div_ceil(n as u128) as u64;
        let mut starts = vec![0u32];
        let (mut count, mut residues, mut oid) = (0usize, 0u64, 0usize);
        while oid < n {
            let lengths = subjects.x_lengths_from(oid).ok()?;
            if lengths.is_empty() {
                return None;
            }
            for length in lengths {
                let length = u64::from(length);
                if count == k || (count > 0 && residues + length > budget) {
                    starts.push(oid as u32);
                    (count, residues) = (0, 0);
                }
                count += 1;
                residues += length;
                oid += 1;
            }
        }
        starts.push(n as u32);
        Some(starts.into())
    }

    /// The number of tasks.
    pub(super) fn count(&self) -> usize {
        match self {
            XTasks::Grid { k, subjects } => subjects.div_ceil(*k),
            XTasks::Cut(starts) => starts.len() - 1,
        }
    }

    /// The subjects of task `task`.
    pub(super) fn range(&self, task: usize) -> std::ops::Range<usize> {
        match self {
            XTasks::Grid { k, subjects } => task * k..(*subjects).min((task + 1) * k),
            XTasks::Cut(starts) => starts[task] as usize..starts[task + 1] as usize,
        }
    }
}

/// No NCBI counterpart: reads the shadow switch; it does not change any value NCBI computes.
/// `LOSAT_X_BXSUBJFLOWSHADOW`: the owner searches every task again on its own table and compares.
fn shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXSUBJFLOWSHADOW").is_some())
}

/// No NCBI counterpart: whether the subject loop runs as a flow, and with how many subjects per
/// task (fewer inside the capped batch flow, `x_bx_mem::flow_subjects_per_task`). Not with
/// diagnostics (they keep every subject's seeds), without a pool, with one task only, or when a
/// diagonal offset could pass 2^30 (module comment).
pub(super) fn chosen(
    pool: &crate::utils::threading::SearchPool<'_>,
    diagnostics: bool,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<usize> {
    // Dispatch point (LOSAT_X_BXCAPFILL, in the strict set): inside the capped batch flow, tasks of
    // fewer subjects (x_bx_mem::flow_subjects_per_task).
    let k = super::x_bx_mem::flow_subjects_per_task(pool, subjects.len(), subjects_per_task()?);
    let own = fits(k, pool, diagnostics, options, subjects);
    crate::utils::xreach::x_reach!("LOSAT_X_BXSUBJFLOW", own);
    own.then_some(k)
}

/// No NCBI counterpart: `chosen`'s rule for tasks of `k` subjects, without recording a side.
pub(super) fn fits(
    k: usize,
    pool: &crate::utils::threading::SearchPool<'_>,
    diagnostics: bool,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> bool {
    !(diagnostics || !pool.enabled() || pool.threads() < 2 || subjects.len() <= k)
        && offsets_stay_below_2_30(options.window_size, subjects.longest_length())
}

/// No NCBI counterpart: the fewest subjects per pool thread with which `LOSAT_X_BXCHUNKFLOW`
/// searches the chunks of a split query one after another (`chunks_in_order`): four tasks of 64.
const CHUNK_FLOW_SUBJECTS_PER_THREAD: usize = 256;
/// No NCBI counterpart: the fewest subjects in all for `chunks_in_order`: 32 tasks of 64.
const CHUNK_FLOW_MIN_SUBJECTS: usize = 2_048;

/// EXPERIMENT (LOSAT_X_BXCHUNKFLOW, in the strict set: on unless `=0`): whether
/// the chunks of a split query are searched one after another, each chunk's subjects as the flow
/// on the whole pool, instead of side by side, one chunk per thread. It holds
/// where the flow runs in each chunk (`chosen`) and there are at least
/// `CHUNK_FLOW_SUBJECTS_PER_THREAD` subjects per thread and `CHUNK_FLOW_MIN_SUBJECTS` subjects.
/// With fewer subjects per thread the flow leaves threads idle in every chunk; with fewer subjects
/// in all, each chunk's setup on the loop's thread (its lookup table, the helpers' start)
/// outweighs the flow at two to four threads. There the chunks side by side are faster. The rule
/// counts subjects, not tasks, so that `LOSAT_X_BXSUBJFLOW=k` does not move it. No NCBI
/// counterpart: it chooses between two schedules of the same chunk loop and changes no value NCBI
/// computes (the sink sees the chunks in order either way, `search::search_core`).
pub(super) fn chunks_in_order(
    pool: &crate::utils::threading::SearchPool<'_>,
    diagnostics: bool,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    if !*ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCHUNKFLOW")) {
        return false;
    }
    let own = subjects_per_task().is_some_and(|k| {
        fits(k, pool, diagnostics, options, subjects)
            && subjects.len()
                >= CHUNK_FLOW_MIN_SUBJECTS.max(CHUNK_FLOW_SUBJECTS_PER_THREAD * pool.threads())
    });
    crate::utils::xreach::x_reach!("LOSAT_X_BXCHUNKFLOW", own);
    own
}

/// EXPERIMENT (LOSAT_X_BXCHUNKCROSS, in the strict set: on unless `=0`): `Some(k)` where the chunk
/// loop of `search::search_core` would search `chunks` (two or more) one after another, each with
/// the subject flow in tasks of `k` (`chosen`): they then run as one flow (`run_chunks`). No
/// condition of its own.
pub(super) fn crosses_chunks(
    pool: &crate::utils::threading::SearchPool<'_>,
    diagnostics: bool,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
    chunks: usize,
) -> Option<usize> {
    static ON: OnceLock<bool> = OnceLock::new();
    if chunks < 2 || !*ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCHUNKCROSS")) {
        return None;
    }
    let k = chosen(pool, diagnostics, options, subjects);
    crate::utils::xreach::x_reach!("LOSAT_X_BXCHUNKCROSS", k.is_some());
    k
}

/// No NCBI counterpart: whether `-window_size` plus twice the longest subject is at most 2^28,
/// so that no diagonal offset of any table history passes 2^30 (module comment).
pub(super) fn offsets_stay_below_2_30(window_size: i32, longest: usize) -> bool {
    u64::try_from(window_size)
        .ok()
        .and_then(|window| window.checked_add(2 * longest as u64))
        .is_some_and(|reach| reach <= 1 << 28)
}

#[cfg(test)]
mod tests {
    use super::{offsets_stay_below_2_30, XTasks};
    use crate::algorithm::blastx::{input::FastaRecord, subjects::Subjects};

    fn records(lengths: &[usize]) -> Vec<FastaRecord> {
        lengths
            .iter()
            .enumerate()
            .map(|(i, &length)| FastaRecord {
                internal_id: format!("s{i}"),
                title: String::new(),
                sequence: vec![b'A'; length],
                lowercase_masks: Vec::new(),
                warnings: Vec::new(),
                source: None,
            })
            .collect()
    }

    // The cut covers the subjects in order, no task holds more than `k` subjects, and a task of
    // two or more subjects holds at most `k` times the mean length; equal lengths keep the grid.
    #[test]
    fn the_tasks_cut_by_residues_cover_the_subjects_in_order() {
        let mut lengths = vec![300; 1_000];
        for length in &mut lengths[400..420] {
            *length = 20_000;
        }
        lengths[999] = 7;
        for k in [1, 2, 7, 64, 2_000] {
            let set = records(&lengths);
            let tasks = XTasks::cut(k, Subjects::Fasta(&set)).unwrap();
            let budget = (k * lengths.iter().sum::<usize>()).div_ceil(lengths.len());
            let mut next = 0;
            for task in 0..tasks.count() {
                let range = tasks.range(task);
                assert_eq!(range.start, next, "k {k} task {task}");
                assert!(!range.is_empty() && range.len() <= k, "k {k} task {task}");
                let residues: usize = lengths[range.clone()].iter().sum();
                assert!(range.len() == 1 || residues <= budget, "k {k} task {task}");
                next = range.end;
            }
            assert_eq!(next, lengths.len(), "k {k}");
            if k == 64 {
                // The twenty long subjects (subjects 400-419 of the grid's task 384..448): two per
                // task under the budget of 64 x the mean length (693.707), 26 tasks for the grid's
                // 16.
                assert_eq!(budget, 44_398);
                assert_eq!(tasks.count(), 26);
                assert_eq!(tasks.range(6), 384..401);
                assert_eq!(tasks.range(7), 401..403);
                assert_eq!(tasks.range(16), 419..483);
            }
        }
        let set = records(&[300; 1_000]);
        let tasks = XTasks::cut(64, Subjects::Fasta(&set)).unwrap();
        assert_eq!(tasks.count(), 16);
        assert_eq!(tasks.range(15), 960..1_000);
        assert!(XTasks::cut(64, Subjects::Fasta(&[])).is_none());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-175
    // ```c++
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    // ```
    // The offset a subject is scanned at is below INT4_MAX / 4 plus the previous subject and a
    // window; its entries add at most the subject's own length.
    #[test]
    fn the_flow_is_kept_only_where_no_offset_passes_2_30() {
        assert!(offsets_stay_below_2_30(40, 5_000_000));
        assert!(offsets_stay_below_2_30((1 << 28) - 2 * 35_000, 35_000));
        assert!(!offsets_stay_below_2_30((1 << 28) - 2 * 35_000 + 1, 35_000));
        assert!(!offsets_stay_below_2_30(536_870_000, 4_910));
        assert!(!offsets_stay_below_2_30(-1, 10));
        // The bound holds for the worst history: INT4_MAX / 4 - 1 plus two subjects and a
        // window stays below 2^30.
        let worst = i64::from(i32::MAX / 4) - 1 + (1 << 28);
        assert!(worst < 1 << 30);
    }
}

/// No NCBI counterpart: the subjects `range` of one task, searched on `diagonals`; see `Task`.
#[allow(clippy::too_many_arguments)]
fn search_task(
    range: std::ops::Range<usize>,
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    db_length: i64,
    x_containment: bool,
) -> Task {
    let mut kept = Vec::new();
    let (first, end) = (range.start, range.end);
    // EXPERIMENT (LOSAT_X_BXDBSTREAM): the task's subjects, read at once when the subjects are read
    // in blocks; an error stops the flow at the task's first subject.
    let x_block = match subjects.x_read_range(first, end) {
        Ok(block) => block,
        Err(error) => return (kept, Some((first, error))),
    };
    let subjects = subjects.x_with_block(x_block.as_ref());
    for oid in first..end {
        match compute_subject(
            subjects,
            batch,
            params,
            options,
            lookup,
            diagonals,
            xdrop,
            false,
            db_length,
            oid,
            x_containment,
        ) {
            Ok(mut computed) => {
                // Only the diagnostic stages read these (not with the flow).
                computed.seeds = Vec::new();
                computed.initial = Vec::new();
                if !computed.is_empty() {
                    kept.push((oid, computed));
                }
            }
            Err(error) => return (kept, Some((oid, error))),
        }
    }
    (kept, None)
}

/// EXPERIMENT (LOSAT_X_BXSUBJFLOWSHADOW): tasks compared so far.
static X_SHADOW_TASKS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

// No NCBI counterpart: prints the shadow-mode summary (from `x_seed_scan::print_shadow_stats`);
// it does not change any value NCBI computes.
pub(crate) fn x_print_shadow_summary() {
    if subjects_per_task().is_some() && shadow() {
        eprintln!(
            "[X_SHADOW] BXSUBJFLOW tasks_compared={} (subject results identical)",
            X_SHADOW_TASKS.load(std::sync::atomic::Ordering::Relaxed)
        );
    }
}

/// No NCBI counterpart: the shadow check of one task: `again` (the owner's second search of it)
/// must give the subjects and the error position of the value the loop used.
fn x_shadow_check(kept: &[(usize, SubjectComputation)], failed_at: Option<usize>, again: &Task) {
    assert_eq!(
        failed_at,
        again.1.as_ref().map(|(oid, _)| *oid),
        "LOSAT_X_BXSUBJFLOWSHADOW: the first refused subject differs"
    );
    assert_eq!(
        kept.len(),
        again.0.len(),
        "LOSAT_X_BXSUBJFLOWSHADOW: the subjects with results differ"
    );
    for ((oa, ca), (ob, cb)) in kept.iter().zip(&again.0) {
        assert!(
            oa == ob
                && ca.raw == cb.raw
                && ca.purged == cb.purged
                && ca.containment == cb.containment
                && ca.prepared_link.map(|(c, f)| (c, f.to_bits()))
                    == cb.prepared_link.map(|(c, f)| (c, f.to_bits())),
            "LOSAT_X_BXSUBJFLOWSHADOW: subject {oa} differs"
        );
    }
    X_SHADOW_TASKS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
/// ```c
///     /* iterate over all subject sequences */
///     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
///            != BLAST_SEQSRC_EOF) {
/// ```
/// The subject loop of one query chunk: every subject of `subjects` in oid order, with the sink
/// calls of the waves' replay loop (`subject_start`, `containment`, `preliminary_purge`,
/// `subject`); `diagonals` is the owner's table, `helper_tables` one table per helper (its own, or
/// with `LOSAT_X_BXDIAGKEEP` the one its thread keeps, `x_bx_mem::XScratch`); `capped` = the chunk's
/// batch is one of the capped batch flow (`x_bx_mem::in_capped_batch`). Returns the number of
/// subjects handed to the sink.
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    k: usize,
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    helper_tables: Vec<super::x_bx_mem::XScratch<'_>>,
    xdrop: i32,
    db_length: i64,
    x_containment: bool,
    capped: bool,
    sink: &mut impl PreliminarySink,
) -> Result<usize> {
    let bounds = XTasks::new(k, subjects);
    let tasks = bounds.count();
    // The fused loop's choice for the lookup table, worked out here once for every helper (a
    // helper's thread may have searched another table's tasks in between).
    let x_ratio =
        crate::algorithm::blastx::seed::x_seed_diag::x_shared_table_ratio(lookup, options);
    let compute = |task: usize, scratch: &mut super::x_bx_mem::XScratch<'_>| {
        crate::algorithm::blastx::seed::x_seed_diag::x_take_table_ratio(lookup, x_ratio);
        let range = bounds.range(task);
        let table = match scratch.get() {
            Ok(table) => table,
            Err(error) => return (Vec::new(), Some((range.start, error))),
        };
        search_task(
            range,
            subjects,
            batch,
            params,
            options,
            lookup,
            table,
            xdrop,
            db_length,
            x_containment,
        )
    };
    // Four tasks per thread ahead of the owner: a helper always finds one while the owner hands a
    // task to the sink, and the tasks held at once stay few (only subjects with hits are kept).
    // Inside the capped batch flow (LOSAT_X_BXUNITCAP), the helpers leave once every task is
    // claimed: the other batches searched at the time want the threads (x_bx_mem.rs).
    let ahead = if capped {
        crate::utils::xahead::Ahead::x_new_leaving(4 * (helper_tables.len() + 1))
    } else {
        crate::utils::xahead::Ahead::new(4 * (helper_tables.len() + 1))
    };
    let x_t = super::super::x_timing::lap();
    let mut call = 0usize;
    let mut own = super::x_bx_mem::XScratch::Borrowed(diagonals);
    let result = crate::utils::xahead::with_helpers(Some(&ahead), helper_tables, &compute, || {
        for task in 0..tasks {
            ahead.announce(task, tasks, |_| true);
            let (kept, error) = match ahead.take(task, &mut own, &compute) {
                Some(value) => value,
                None => compute(task, &mut own),
            };
            if shadow() {
                // The owner searches the task again on its own table (one that has seen other
                // tasks than the helper's), also when it searched the task itself.
                x_shadow_check(
                    &kept,
                    error.as_ref().map(|(oid, _)| *oid),
                    &compute(task, &mut own),
                );
            }
            let range = bounds.range(task);
            let end = error.as_ref().map_or(range.end, |(oid, _)| *oid);
            let mut kept = kept.into_iter().peekable();
            for oid in range.start..end {
                let computed = match kept.next_if(|(kept_oid, _)| *kept_oid == oid) {
                    Some((_, computed)) => computed,
                    None => SubjectComputation::empty(),
                };
                sink.subject_start(oid, batch, params, computed.prepared_link)?;
                for (h, contained) in &computed.containment {
                    sink.containment(h, *contained);
                }
                let SubjectComputation { raw, purged, .. } = computed;
                if options.gapped && !raw.is_empty() {
                    sink.preliminary_purge(oid, &raw, &purged)?;
                }
                sink.subject(oid, &purged, batch, params)?;
                call += 1;
            }
            if let Some((_, error)) = error {
                return Err(error);
            }
        }
        Ok(())
    });
    super::super::x_timing::subject_flow(x_t, tasks as u64, ahead.counts());
    result.map(|()| call)
}

/// No NCBI counterpart: tasks per thread that the flow across chunks searches ahead of the loop's
/// thread (`run_chunks`; scheduling only).
const X_CROSS_AHEAD: usize = 16;

/// No NCBI counterpart: what the subject loop of one query chunk needs, made by the statements of
/// `search::search_core` before that loop (`LOSAT_X_BXCHUNKCROSS`): the chunk's Karlin blocks with
/// its subject parameters, its lookup table, its X-drop, and a diagonal table of its length for
/// the loop's thread.
pub(super) struct XChunkSetup {
    pub(super) params: Vec<ContextParameters>,
    pub(super) lookup: Lookup,
    pub(super) xdrop: i32,
    pub(super) table: Diagonals,
    pub(super) length: usize,
}

// No NCBI counterpart: a chunk's setup in the flow across chunks, with the fused loop's choice for
// its lookup table (worked out once, as `run` does) and the length of its diagonal tables.
struct XChunkReady {
    params: Vec<ContextParameters>,
    lookup: Lookup,
    xdrop: i32,
    ratio: Option<u64>,
    length: usize,
    table: std::sync::Mutex<Option<Diagonals>>,
}

// No NCBI counterpart: the state of one chunk's setup in the flow across chunks. `Failed` keeps the
// error until the loop's thread reaches the chunk and returns it; `Ended`: the loop has ended the
// chunk and dropped its setup.
enum XSetupState {
    Empty,
    Building,
    Ready(std::sync::Arc<XChunkReady>),
    Skipped,
    Failed(Option<anyhow::Error>),
    Ended,
}

// No NCBI counterpart: what a thread finds for a chunk.
enum XGot {
    Ready(std::sync::Arc<XChunkReady>),
    Skipped,
    Failed,
}

// No NCBI counterpart: a mutex lock that ignores poisoning (no lock is held while a setup is made).
fn x_lock<T>(mutex: &std::sync::Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

// No NCBI counterpart: the setups of the chunks of one flow across chunks, each made once by the
// first thread that needs it (`get`) or that starts it ahead (`start`).
struct XSetups<'s> {
    slots: Vec<std::sync::Mutex<XSetupState>>,
    make: &'s (dyn Fn(usize, usize) -> Result<Option<XChunkSetup>> + Sync),
    options: &'s ResolvedOptions,
}

// No NCBI counterpart: puts a chunk's setup back to `Empty` if making it unwinds, so that the next
// thread that needs it makes it again (and meets the same panic) instead of waiting forever.
struct XUnbuilt<'a>(&'a std::sync::Mutex<XSetupState>, bool);

impl Drop for XUnbuilt<'_> {
    fn drop(&mut self) {
        if self.1 {
            *x_lock(self.0) = XSetupState::Empty;
        }
    }
}

impl XSetups<'_> {
    // No NCBI counterpart: makes chunk `c`'s setup on this thread (its slot is `Building`), its
    // lookup table in `parts`.
    fn make(&self, c: usize, parts: usize) {
        let slot = &self.slots[c];
        let mut unbuilt = XUnbuilt(slot, true);
        let state = match (self.make)(c, parts) {
            Ok(Some(setup)) => {
                let ratio = crate::algorithm::blastx::seed::x_seed_diag::x_shared_table_ratio(
                    &setup.lookup,
                    self.options,
                );
                XSetupState::Ready(std::sync::Arc::new(XChunkReady {
                    params: setup.params,
                    lookup: setup.lookup,
                    xdrop: setup.xdrop,
                    ratio,
                    length: setup.length,
                    table: std::sync::Mutex::new(Some(setup.table)),
                }))
            }
            Ok(None) => XSetupState::Skipped,
            Err(error) => XSetupState::Failed(Some(error)),
        };
        unbuilt.1 = false;
        *x_lock(slot) = state;
    }

    // No NCBI counterpart: starts chunk `c`'s setup on this thread if nobody has (a helper ahead of
    // the chunk), its lookup table in `parts`; returns at once if it is being made or made.
    fn start(&self, c: usize, parts: usize) {
        {
            let mut state = x_lock(&self.slots[c]);
            if !matches!(*state, XSetupState::Empty) {
                return;
            }
            *state = XSetupState::Building;
        }
        self.make(c, parts);
    }

    // No NCBI counterpart: chunk `c`'s setup: made here if nobody has started it, else the one
    // another thread makes (waited for: a setup takes a few milliseconds).
    fn get(&self, c: usize) -> XGot {
        let mut idle = 0u32;
        loop {
            let mine = {
                let mut state = x_lock(&self.slots[c]);
                match &*state {
                    XSetupState::Ready(ready) => return XGot::Ready(std::sync::Arc::clone(ready)),
                    XSetupState::Skipped => return XGot::Skipped,
                    XSetupState::Failed(_) | XSetupState::Ended => return XGot::Failed,
                    XSetupState::Building => false,
                    XSetupState::Empty => {
                        *state = XSetupState::Building;
                        true
                    }
                }
            };
            if mine {
                self.make(c, 1);
                continue;
            }
            idle = idle.saturating_add(1);
            if idle < 128 {
                std::hint::spin_loop();
            } else {
                std::thread::yield_now();
            }
        }
    }

    // No NCBI counterpart: the error of chunk `c`'s setup, once.
    fn error(&self, c: usize) -> Option<anyhow::Error> {
        match &mut *x_lock(&self.slots[c]) {
            XSetupState::Failed(error) => error.take(),
            _ => None,
        }
    }

    // No NCBI counterpart: the loop has ended chunk `c`: its setup is dropped (once the helpers
    // holding it let go).
    fn end(&self, c: usize) {
        *x_lock(&self.slots[c]) = XSetupState::Ended;
    }
}

// No NCBI counterpart (buffer reuse): a thread's diagonal table in the flow across chunks, made again
// when a task's chunk has another length (the last chunk of a split query is shorter): the table
// its thread keeps (`LOSAT_X_BXDIAGKEEP`, `x_bx_mem::XDiagSlot`) or one of its own.
enum XCrossScratch {
    Slot(usize, Option<super::x_bx_mem::XDiagSlot>),
    Owned(usize, Option<Diagonals>),
}

impl XCrossScratch {
    // No NCBI counterpart: whether the table serves a chunk of `length` letters (`Diagonals::x_fits`).
    fn fits(&self, length: usize, window: i32) -> bool {
        match self {
            Self::Slot(have, slot) => *have == length && slot.is_some(),
            Self::Owned(_, table) => table
                .as_ref()
                .is_some_and(|table| table.x_fits(length, window)),
        }
    }

    // No NCBI counterpart: buffer reuse only (a new table is `Diagonals::new`, the NCBI port).
    fn get(&mut self, length: usize, window: i32) -> Result<&mut Diagonals> {
        let fits = self.fits(length, window);
        match self {
            Self::Slot(have, slot) => {
                if !fits {
                    // The slot left gives its table back to this thread first.
                    *slot = None;
                    *slot = Some(super::x_bx_mem::XDiagSlot::new(length, window));
                    *have = length;
                }
                slot.as_mut().expect("set above").get()
            }
            Self::Owned(have, table) => {
                if !fits {
                    *table = Some(Diagonals::new(length, window)?);
                    *have = length;
                }
                Ok(table.as_mut().expect("set above"))
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-271
/// ```c
///         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
/// ...
///                     retval =
///                         CPrelimSearchRunner(*chunk_data, opts_memento.get())();
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
/// ```c
///     /* iterate over all subject sequences */
///     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
///            != BLAST_SEQSRC_EOF) {
/// ```
/// EXPERIMENT (LOSAT_X_BXCHUNKCROSS, in the strict set): the chunk loop of a split query whose chunks
/// `search::search_core` searches one after another, each with the subject flow
/// (`LOSAT_X_BXCHUNKFLOW`), as one flow of tasks across the chunks. Task `c * tasks + t` is task `t`
/// (at most `k` subjects, `XTasks`) of chunk `c`; the helpers search up to a window ahead of this thread across a
/// chunk's end, so they do not wait while this thread hands the last tasks of a chunk to the sink
/// and the next chunk is set up; the tasks offered are those of this thread's chunk and the next.
/// A chunk's setup (`make(c, parts)`: `search_core`'s statements before the chunk's subject loop, a
/// function of the chunk alone) is made once, by the first thread that needs it or ahead by a
/// helper near the end of this thread's chunk, on one thread (so at most two chunks' setups are
/// alive at a time), and the first chunk's on the pool before the helpers start (its lookup table
/// in `first_parts`, `LOSAT_X_BXBUILDPAR`). This thread makes the sink calls of the loop in
/// its order: per chunk `chunk_start`, the subjects as `run` hands them, `chunk_end`; a chunk
/// without a valid context gets `chunk_start` only; a setup that failed returns its error after
/// that chunk's `chunk_start`. A task's value depends on its chunk and subjects alone (module
/// comment), so only which thread computes what, and when, changes. Returns the number of subjects
/// handed to the sink.
#[allow(clippy::too_many_arguments)]
pub(super) fn run_chunks(
    k: usize,
    chunks: &[crate::algorithm::blastx::split::QueryChunk],
    make: &(dyn Fn(usize, usize) -> Result<Option<XChunkSetup>> + Sync),
    first_parts: usize,
    subjects: Subjects<'_>,
    options: &ResolvedOptions,
    helpers: usize,
    keep: bool,
    db_length: i64,
    x_containment: bool,
    capped: bool,
    is_split: bool,
    sink: &mut impl PreliminarySink,
) -> Result<usize> {
    let bounds = XTasks::new(k, subjects);
    let tasks = bounds.count();
    let total = chunks.len() * tasks;
    let setups = XSetups {
        slots: (0..chunks.len())
            .map(|_| std::sync::Mutex::new(XSetupState::Empty))
            .collect(),
        make,
        options,
    };
    // The first chunk's setup while the pool is idle, before the helpers start: the parts of its
    // lookup table's build (`first_parts`) go to the other threads. A setup made later is made on
    // one thread: inside a helper, a build in parts could, waiting at its end for a part another
    // thread took, run a helper's loop that is not started yet and wait there for this very setup.
    setups.start(0, first_parts);
    // A helper this many tasks before the end of the loop's chunk makes the next chunk's setup.
    let lead = 2 * (helpers + 1);
    let window = options.window_size;
    // The chunk the loop's thread is at: the tasks offered are those of this chunk and the next
    // (`x_announce_until`), and only the next chunk's setup is made ahead, so at most two chunks'
    // setups are alive at a time.
    let at_chunk = std::sync::atomic::AtomicUsize::new(0);
    let x_probe = std::env::var_os("LOSAT_X_S17STAMP").is_some();
    let x_tasks = std::sync::Mutex::new(Vec::<(usize, usize, u64, u64, usize, usize)>::new());
    let compute_inner = |index: usize, scratch: &mut XCrossScratch| -> Task {
        let (c, task) = (index / tasks, index % tasks);
        // A chunk without a valid context (or whose setup failed) has no task to search; the loop
        // does not take its values.
        let XGot::Ready(ready) = setups.get(c) else {
            return (Vec::new(), None);
        };
        if task + lead >= tasks
            && c + 1 < chunks.len()
            && c == at_chunk.load(std::sync::atomic::Ordering::Relaxed)
        {
            setups.start(c + 1, 1);
        }
        crate::algorithm::blastx::seed::x_seed_diag::x_take_table_ratio(&ready.lookup, ready.ratio);
        let range = bounds.range(task);
        let table = match scratch.get(ready.length, window) {
            Ok(table) => table,
            Err(error) => return (Vec::new(), Some((range.start, error))),
        };
        search_task(
            range,
            subjects,
            &chunks[c].prepared,
            &ready.params,
            options,
            &ready.lookup,
            table,
            ready.xdrop,
            db_length,
            x_containment,
        )
    };
    let compute = |index: usize, scratch: &mut XCrossScratch| -> Task {
        if !x_probe {
            return compute_inner(index, scratch);
        }
        let t0 = crate::utils::xtiming::x_since_us();
        let value = compute_inner(index, scratch);
        let t1 = crate::utils::xtiming::x_since_us();
        let hsps: usize = value.0.iter().map(|(_, c)| c.purged.len()).sum();
        x_lock(&x_tasks).push((index, rayon::current_thread_index().unwrap_or(999), t0, t1, value.0.len(), hsps));
        value
    };
    crate::utils::xtiming::x_stamp("flow_start");
    // Sixteen tasks per thread ahead of this thread (`run` keeps four): a task much slower than
    // the others (a subject with many hits) holds this thread at its place in the order while the
    // helpers search the rest of the window, and with four per thread they ran out of tasks and
    // spun (16 threads, fly 2L 50 kb x the fly proteins: 1.6 s of `sched_yield` in 0.69 s); with
    // sixteen they seldom do (0.08 s), and the tasks held at once stay few next to the lookup
    // tables (only subjects with hits are kept). Inside the capped batch flow the helpers leave
    // once every task is claimed.
    let ahead = if capped {
        crate::utils::xahead::Ahead::x_new_leaving(X_CROSS_AHEAD * (helpers + 1))
    } else {
        crate::utils::xahead::Ahead::new(X_CROSS_AHEAD * (helpers + 1))
    };
    let scratches: Vec<_> = (0..helpers)
        .map(|_| {
            if keep {
                XCrossScratch::Slot(0, None)
            } else {
                XCrossScratch::Owned(0, None)
            }
        })
        .collect();
    let x_t = super::super::x_timing::lap();
    let mut call = 0usize;
    let mut own = XCrossScratch::Owned(0, None);
    let result = crate::utils::xahead::with_helpers(Some(&ahead), scratches, &compute, || {
        for (c, chunk) in chunks.iter().enumerate() {
            at_chunk.store(c, std::sync::atomic::Ordering::Relaxed);
            sink.chunk_start(chunk)?;
            let ready = match setups.get(c) {
                XGot::Ready(ready) => ready,
                // No valid context: the serial loop goes on to the next chunk.
                XGot::Skipped => continue,
                XGot::Failed => {
                    return Err(setups
                        .error(c)
                        .expect("a failed setup keeps its error until the loop reaches it"))
                }
            };
            if !own.fits(ready.length, window) {
                if let Some(table) = x_lock(&ready.table).take() {
                    own = XCrossScratch::Owned(ready.length, Some(table));
                }
            }
            let batch = &chunk.prepared;
            let params = &ready.params;
            for task in 0..tasks {
                let index = c * tasks + task;
                ahead.x_announce_until(index, ((c + 2) * tasks).min(total), total, |_| true);
                let (kept, error) = match ahead.take(index, &mut own, &compute) {
                    Some(value) => value,
                    None => compute(index, &mut own),
                };
                if shadow() {
                    x_shadow_check(
                        &kept,
                        error.as_ref().map(|(oid, _)| *oid),
                        &compute(index, &mut own),
                    );
                }
                let range = bounds.range(task);
                let end = error.as_ref().map_or(range.end, |(oid, _)| *oid);
                let mut kept = kept.into_iter().peekable();
                for oid in range.start..end {
                    let computed = match kept.next_if(|(kept_oid, _)| *kept_oid == oid) {
                        Some((_, computed)) => computed,
                        None => SubjectComputation::empty(),
                    };
                    sink.subject_start(oid, batch, params, computed.prepared_link)?;
                    for (h, contained) in &computed.containment {
                        sink.containment(h, *contained);
                    }
                    let SubjectComputation { raw, purged, .. } = computed;
                    if options.gapped && !raw.is_empty() {
                        sink.preliminary_purge(oid, &raw, &purged)?;
                    }
                    sink.subject(oid, &purged, batch, params)?;
                    call += 1;
                }
                if let Some((_, error)) = error {
                    return Err(error);
                }
            }
            drop(ready);
            setups.end(c);
            sink.chunk_end(chunk, is_split)?;
        }
        Ok(())
    });
    super::super::x_timing::subject_flow(x_t, total as u64, ahead.counts());
    if x_probe {
        use std::fmt::Write as _;
        let mut text = crate::utils::xtiming::X_DEFERRED
            .lock()
            .unwrap_or_else(|e| e.into_inner());
        for (index, thread, t0, t1, kept, hsps) in x_lock(&x_tasks).iter() {
            let _ = writeln!(text, "XTASK {index} {thread} {t0} {t1} {kept} {hsps}");
        }
    }
    result.map(|()| call)
}
