//! EXPERIMENT (LOSAT_X_BXSUBJFLOW[=k], in the strict set: on unless `=0` or `LOSAT_X_DEFAULT=0`):
//! the subject loop of a query chunk as a flow of tasks of `k` consecutive subjects
//! (`search::X_BX_CHUNK` = 64, the waves' subjects per worker, when the variable is `1` or not a
//! number from 2 to 65,536), searched ahead by the other threads of the pool and handed to the sink
//! in subject order by the thread that runs the loop. Shadow: `LOSAT_X_BXSUBJFLOWSHADOW`.
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

/// No NCBI counterpart: reads the shadow switch; it does not change any value NCBI computes.
/// `LOSAT_X_BXSUBJFLOWSHADOW`: the owner searches every task again on its own table and compares.
fn shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXSUBJFLOWSHADOW").is_some())
}

/// No NCBI counterpart: whether the subject loop runs as a flow, and with how many subjects per
/// task. Not with diagnostics (they keep every subject's seeds), without a pool, with one task
/// only, or when a diagonal offset could pass 2^30 (module comment).
pub(super) fn chosen(
    pool: &crate::utils::threading::SearchPool<'_>,
    diagnostics: bool,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<usize> {
    let k = subjects_per_task()?;
    let own = fits(k, pool, diagnostics, options, subjects);
    crate::utils::xreach::x_reach!("LOSAT_X_BXSUBJFLOW", own);
    own.then_some(k)
}

/// No NCBI counterpart: `chosen`'s rule for tasks of `k` subjects, without recording a side.
fn fits(
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

/// EXPERIMENT (LOSAT_X_BXCHUNKFLOW, in the strict set: on unless `=0` or `LOSAT_X_DEFAULT=0`): whether
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
    use super::offsets_stay_below_2_30;

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

/// No NCBI counterpart: the subjects of task `task`, searched on `diagonals`; see `Task`.
#[allow(clippy::too_many_arguments)]
fn search_task(
    task: usize,
    k: usize,
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
    let (first, end) = (task * k, subjects.len().min((task + 1) * k));
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
/// `subject`); `diagonals` is the owner's table, `helper_tables` one table per helper. Returns the
/// number of subjects handed to the sink.
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
    sink: &mut impl PreliminarySink,
) -> Result<usize> {
    let tasks = subjects.len().div_ceil(k);
    let compute = |task: usize, scratch: &mut super::x_bx_mem::XScratch<'_>| {
        let table = match scratch.get() {
            Ok(table) => table,
            Err(error) => return (Vec::new(), Some((task * k, error))),
        };
        search_task(
            task,
            k,
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
    let ahead = crate::utils::xahead::Ahead::new(4 * (helper_tables.len() + 1));
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
            let end = error
                .as_ref()
                .map_or(subjects.len().min((task + 1) * k), |(oid, _)| *oid);
            let mut kept = kept.into_iter().peekable();
            for oid in task * k..end {
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
