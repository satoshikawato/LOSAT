//! EXPERIMENT (LOSAT_X_BPSUBJRUN / LOSAT_X_BPSUBJRUNSHADOW): the parallel preliminary loop of a
//! `-subject` search in runs of consecutive subjects, one diagonal table fill per run.
//!
//! Without the switch, the parallel loop of `run_resolved_in_pool` searches each subject on its own:
//! `for_each_init` makes a new `BlastpSubjectScratch` for every piece rayon splits the subjects into
//! (about 2,300 per query block of the human x yeast blastp-fast job at 16 threads), and each subject
//! starts with `prepare_independent_subject`, which fills the whole diagonal table (65,536 entries of
//! 4 bytes for a 50,000-letter query block). In that job at 16 threads (workstation,
//! `LOSAT_X_BPQBLOCK=50000`) the scratches took 3.7 s and the fills 7.4 s of the loop's 39 s of CPU.
//!
//! With the switch on, the subjects are cut into contiguous runs (`x_cut_runs`: guided
//! self-scheduling, each run the next `1 / (X_PARTS_PER_THREAD x threads)` of the letters left, so
//! the runs shrink to single subjects at the end), and the search runs one task per thread of the
//! pool in one `rayon::scope`; each task makes one scratch when it takes its first run and takes
//! the runs, in order, from one counter. A run starts as a subject of the loop without the switch
//! does: `prepare_independent_subject` at the diagonal offset that the one-thread order gives its
//! first subject (`blastp_subject_diag_offsets`). The later subjects of the run are searched as the
//! one-thread loop searches them: the table is not filled again, and the offset is the one that
//! the previous subject's `finish_subject` left (`Blast_ExtendWordExit`, below), which is the
//! one-thread order's offset of the subject, because the run holds consecutive subjects and
//! `blastp_subject_diag_offsets` is the same recurrence (a subject without letters is not searched
//! and moves neither). Every other per-subject reset of `prepare_independent_subject`
//! (`BlastpSubjectScratch::reset_subject_state`: the initial HSPs, the chaining nodes, the
//! preliminary HSP lists, the restricted-alignment arrays, the interval tree) still runs before
//! each subject. The results do not need it (`process_subject` takes, clears or resets each of
//! these fields itself, which the database path relies on for the later subjects of a block); it
//! is kept so that a later subject of a run starts from the state of the loop without the switch,
//! field by field, except for the diagonal table. Each subject's result goes to its own slot, and the slots are merged in subject
//! order afterwards as without the switch (`merge_blastp_subject_results_in_order`).
//!
//! Why the reports do not change. A subject's result depends on the diagonal table only through
//! the test of each word hit (`aa_ungapped.c`, below; the fused loops of `LOSAT_X_BPSEED` make
//! the same decisions). For an entry with flag 0, a cleared entry (last_hit 0 at an offset of at
//! least `window`, or `-window` after `s_BlastDiagClear` at the offset `window`) gives
//! `diff >= window` for every subject offset (`diff >= 0` for one hit, where `window` is 0), so the
//! hit only sets `last_hit` to `subject_offset + diag_offset`; an entry with flag 0 and
//! `last_hit <= diag_offset - window` does the same. An entry with flag 1 and
//! `last_hit <= diag_offset` takes the "start a new hit" branch for every subject offset, which
//! also sets `last_hit` to `subject_offset + diag_offset` and the flag to 0, as the cleared entry
//! does. Every entry a subject of length L writes at offset `d` has `last_hit <= d + L - 1`
//! (`subject_offset + diag_offset` with `subject_offset <= L - wordsize`, or
//! `s_last_off - (wordsize - 1) + diag_offset` with `s_last_off <= L` and `wordsize >= 2`), and
//! the next subject's offset is `d + L + window`, so at the start of every
//! later subject of a run each entry reads as a cleared one: the subject's decisions, and so its
//! result, are those of the loop without the switch, which starts it from a cleared table at the
//! same offset. When the offset reaches INT4_MAX / 4, `finish_subject` puts it back to `window`
//! and clears the table (`s_BlastDiagClear`) inside the run, as the one-thread loop does, and
//! `blastp_subject_diag_offsets` gives the next subject the offset `window` too. This is the
//! property the database path uses for its blocks (`database::search_blocks`).
//!
//! The argument is in integers; `last_hit` is a 31-bit field (`DiagStruct`), whose stored values
//! are exact while `d + L < 2^30`. The offsets stay below INT4_MAX / 4 plus one subject and one
//! window, so this holds for every subject shorter than 2^29 letters. Only one subject of more than
//! 2^29 letters in one-hit mode could leave a wrapped entry that a later subject of its run reads
//! differently from a cleared one; there the run path keeps the leftovers as the one-thread loop
//! and NCBI do, and it is the loop without the switch that departs from them. The shadow's
//! recompute would show such a difference; its table check, which reads the sign-extended field,
//! would not.
//!
//! Which thread
//! searches a run and in what order the runs end change nothing: the results are merged by
//! subject index, and no floating-point operation changes order.
//!
//! A strict-set switch (`utils/xswitch.rs`): on unless `LOSAT_X_BPSUBJRUN=0`.
//! It acts only where
//! the parallel loop runs: a `-subject` search with more than one thread whose subjects
//! `blastp_should_parallelize_stage` gives to the pool. The one-thread loop and the database path
//! (`-db`) already search with one table and are not changed.
//!
//! `LOSAT_X_BPSUBJRUNSHADOW=1` (shadow mode, implies the switch): before each later subject of a
//! run, the run's offset is compared with `blastp_subject_diag_offsets` (a wrong offset) and every
//! entry of the table is checked to read as a cleared one at that offset (`x_first_uncleared`: an
//! entry left above the offset, as a reset that did not clear the table would leave, fails); then
//! every subject is searched again with a second scratch of the task, as the loop without the
//! switch searches it (`prepare_independent_subject` at the same offset), and the two results are
//! compared (every HSP list of every query, every field). The first difference stops the run
//! (panic) after printing the line that is otherwise printed at exit:
//! `[X_SHADOW] BPSUBJRUN searches= runs= subjects= fills_saved= resets_in_runs= subjects_compared=
//! tables_checked= identical` (`fills_saved`: subjects searched without a fill of their own;
//! `resets_in_runs`: INT4_MAX / 4 resets that happened inside a run; `tables_checked`: tables found
//! to read as cleared before a later subject of a run).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:160-172
//! ```c
//! Int2
//! Blast_ExtendWordExit(Blast_ExtendWord * ewp, Int4 subject_length)
//! {
//!     if (!ewp)
//!         return -1;
//!
//!     if (ewp->diag_table) {
//!         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//!             ewp->diag_table->offset = ewp->diag_table->window;
//!             s_BlastDiagClear(ewp->diag_table);
//!         } else {
//!             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//!         }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-103
//! ```c
//! static Int4 s_BlastDiagClear(BLAST_DiagTable * diag)
//! {
//!     Int4 i, n;
//!     DiagStruct *diag_struct_array;
//!
//!     if (diag == NULL)
//!         return 0;
//!
//!     n = diag->diag_array_length;
//!
//!     diag->offset = diag->window;
//!
//!     diag_struct_array = diag->hit_level_array;
//!
//!     for (i = 0; i < n; i++) {
//!         diag_struct_array[i].flag = 0;
//!         diag_struct_array[i].last_hit = -diag->window;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-544
//! ```c
//!             diag_coord = (query_offset - subject_offset) & diag_mask;
//!
//!             /* If the reset bit is set, an extension just happened. */
//!             if (diag_array[diag_coord].flag) {
//!                 /* If we've already extended past this hit, skip it. */
//!                 if ((Int4) (subject_offset + diag_offset) <
//!                     diag_array[diag_coord].last_hit) {
//!                     continue;
//!                 }
//!                 /* Otherwise, start a new hit. */
//!                 else {
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     diag_array[diag_coord].flag = 0;
//!                 }
//!             }
//!             /* If the reset bit is cleared, try to start an extension. */
//!             else {
//!                 /* find the distance to the last hit on this diagonal */
//!                 last_hit = diag_array[diag_coord].last_hit - diag_offset;
//!                 diff = subject_offset - last_hit;
//!
//!                 if (diff >= window) {
//!                     /* We are beyond the window for this diagonal; start a
//!                        new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//!                 }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413,1553-1554
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ...
//!          /* Save the results. */
//!          status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
//! ```

use std::ops::Range;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::utils::threading::SearchPool;

use super::{BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, EncodedProtein};

/// EXPERIMENT (LOSAT_X_BPSUBJRUN): the mode the variables select.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XSubjRun {
    Off,
    On,
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BPSUBJRUN / LOSAT_X_BPSUBJRUNSHADOW switches once; it does
// not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPSUBJRUN): a strict-set switch (on unless it is `0`);
/// `LOSAT_X_BPSUBJRUNSHADOW` (any value, not part of the set) selects the
/// shadow mode.
pub(super) fn x_bp_subj_run() -> XSubjRun {
    static MODE: OnceLock<XSubjRun> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPSUBJRUNSHADOW").is_some() {
            XSubjRun::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BPSUBJRUN") {
            XSubjRun::On
        } else {
            XSubjRun::Off
        }
    })
}

// No NCBI counterpart: how the subjects are cut into runs; it changes no value NCBI computes (any
// cut gives the same results).
/// The parts of the remaining work per thread that the next run takes (`x_cut_runs`). A run's fill
/// costs a few microseconds, but the runs that end last decide when the loop ends: runs of equal
/// size cost exact blastp up to 2.5 % of its time at 16 threads on the workstation (a run of 1/256
/// of a block's subjects is about 13 ms of its search there), while runs that shrink with the work
/// left, down to single subjects, end together (perf-mb-20261014 `logs/subjrun_report.md`).
const X_PARTS_PER_THREAD: usize = 4;

// No NCBI counterpart: the cut of the subjects into runs; it changes no value NCBI computes.
/// Contiguous runs of the subject indices `0..lengths.len()`, cut as guided self-scheduling cuts
/// work: each run is the shortest run from where the previous one ended whose weight is at least
/// the weight left divided by `parts` (rounded up), so the runs shrink as the work left shrinks,
/// down to single subjects. A subject weighs its letters and one more, so that subjects without
/// letters count too. Every run holds at least one subject, and the runs cover every index once,
/// in order. The cut depends only on the lengths and `parts`, not on which thread takes a run.
fn x_cut_runs(lengths: &[usize], parts: usize) -> Vec<Range<usize>> {
    let n = lengths.len();
    let parts = parts.max(1) as u64;
    // prefix[i]: the weight of the subjects before i.
    let mut prefix = Vec::with_capacity(n + 1);
    let mut weight = 0u64;
    prefix.push(weight);
    for &length in lengths {
        weight += length as u64 + 1;
        prefix.push(weight);
    }
    let mut runs = Vec::new();
    let mut start = 0;
    while start < n {
        let wanted = prefix[start] + (prefix[n] - prefix[start]).div_ceil(parts);
        // The first end whose run reaches the wanted weight (there is one: prefix[n] does).
        let end = prefix
            .partition_point(|&before| before < wanted)
            .max(start + 1);
        runs.push(start..end);
        start = end;
    }
    runs
}

static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_RUNS: AtomicU64 = AtomicU64::new(0);
static X_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static X_FILLS_SAVED: AtomicU64 = AtomicU64::new(0);
static X_RESETS_IN_RUNS: AtomicU64 = AtomicU64::new(0);
static X_SUBJECTS_COMPARED: AtomicU64 = AtomicU64::new(0);
static X_TABLES_CHECKED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BPSUBJRUN searches={} runs={} subjects={} fills_saved={} resets_in_runs={} subjects_compared={} tables_checked={}{}",
        X_SEARCHES.load(Ordering::Relaxed),
        X_RUNS.load(Ordering::Relaxed),
        X_SUBJECTS.load(Ordering::Relaxed),
        X_FILLS_SAVED.load(Ordering::Relaxed),
        X_RESETS_IN_RUNS.load(Ordering::Relaxed),
        X_SUBJECTS_COMPARED.load(Ordering::Relaxed),
        X_TABLES_CHECKED.load(Ordering::Relaxed),
        if identical { " identical" } else { " DIFFERENT" },
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPSUBJRUNSHADOW): the line `main` prints at exit.
pub(crate) fn x_bp_subj_run_print_summary() {
    if x_bp_subj_run() == XSubjRun::Shadow {
        x_shadow_print_line(true);
    }
}

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    // The panic ends the run before `main` prints the line at exit: it is printed here.
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BPSUBJRUN {message}");
}

// No NCBI counterpart: a check of the shadow mode (the property in the module's header); it
// changes no value NCBI computes.
/// The first entry of `table` that does not read as a cleared entry at the diagonal offset
/// `diag_offset` for every subject offset: an entry with flag 0 must have
/// `last_hit <= diag_offset - window`, an entry with flag 1 `last_hit <= diag_offset`.
fn x_first_uncleared(table: &[DiagStruct], diag_offset: i32, window: i32) -> Option<usize> {
    let flag0_bound = i64::from(diag_offset) - i64::from(window);
    table.iter().position(|entry| {
        let last_hit = i64::from(entry.last_hit());
        if entry.flag() == 0 {
            last_hit > flag0_bound
        } else {
            last_hit > i64::from(diag_offset)
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
// ```c
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-172
// ```c
//         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//             ewp->diag_table->offset = ewp->diag_table->window;
//             s_BlastDiagClear(ewp->diag_table);
//         } else {
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//         }
// ```
/// The subjects of one run, `run`, searched with `scratch` into `slots` (the run's slots of the
/// subject results), and with `reference` (shadow mode) searched again as the loop without the
/// switch searches them and compared.
#[allow(clippy::too_many_arguments)]
fn x_search_run<K>(
    scratch: &mut BlastpSubjectScratch,
    mut reference: Option<&mut BlastpSubjectScratch>,
    run: Range<usize>,
    slots: &mut [Option<BlastpSubjectResult>],
    subjects: &[EncodedProtein],
    subject_diag_offsets: &[i32],
    window: i32,
    process_subject: &K,
) where
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let shadow = reference.is_some();
    let (mut resets, mut checked) = (0u64, 0u64);
    for (s_idx, slot) in run.clone().zip(slots.iter_mut()) {
        let diag_offset = subject_diag_offsets[s_idx];
        if s_idx == run.start {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-103
            // (`s_BlastDiagClear`, module header): the run starts from a cleared table at the
            // one-thread order's offset of its first subject, as each subject does without the
            // switch.
            scratch.prepare_independent_subject(diag_offset);
        } else {
            // The previous subject's `finish_subject` (`Blast_ExtendWordExit`) left the
            // one-thread order's offset, and its entries read as cleared ones (module header).
            if shadow {
                if scratch.diag_offset != diag_offset {
                    x_shadow_difference(format!(
                        "subject {s_idx}: the run's diagonal offset is {}, the one-thread order's {diag_offset}",
                        scratch.diag_offset
                    ));
                }
                if let Some(entry) = x_first_uncleared(&scratch.diag_array, diag_offset, window) {
                    x_shadow_difference(format!(
                        "subject {s_idx}: diagonal entry {entry} (last_hit {}, flag {}) does not read as cleared at offset {diag_offset}",
                        scratch.diag_array[entry].last_hit(),
                        scratch.diag_array[entry].flag()
                    ));
                }
                checked += 1;
            } else {
                debug_assert_eq!(scratch.diag_offset, diag_offset, "subject {s_idx}");
            }
            if diag_offset < subject_diag_offsets[s_idx - 1] {
                resets += 1;
            }
            // Every per-subject reset of `prepare_independent_subject` but the table's fill.
            scratch.reset_subject_state();
        }
        let mut subject_result = BlastpSubjectResult::new(s_idx);
        process_subject(
            scratch,
            s_idx,
            BlastpSubjectRef::of(&subjects[s_idx]),
            &mut subject_result,
        );
        if let Some(reference) = reference.as_deref_mut() {
            reference.prepare_independent_subject(diag_offset);
            let mut expected = BlastpSubjectResult::new(s_idx);
            process_subject(
                reference,
                s_idx,
                BlastpSubjectRef::of(&subjects[s_idx]),
                &mut expected,
            );
            let lists = format!("{:?}", subject_result.query_hsp_lists);
            let expected_lists = format!("{:?}", expected.query_hsp_lists);
            let compact = format!("{:?}", subject_result.x_compact_hsp_lists);
            let expected_compact = format!("{:?}", expected.x_compact_hsp_lists);
            if lists != expected_lists || compact != expected_compact {
                x_shadow_difference(format!(
                    "subject {s_idx}: the run's result differs from the result of a cleared table: {lists} {compact} / {expected_lists} {expected_compact}"
                ));
            }
            X_SUBJECTS_COMPARED.fetch_add(1, Ordering::Relaxed);
        }
        *slot = Some(subject_result);
    }
    if shadow {
        X_RUNS.fetch_add(1, Ordering::Relaxed);
        X_SUBJECTS.fetch_add(run.len() as u64, Ordering::Relaxed);
        X_FILLS_SAVED.fetch_add(run.len() as u64 - 1, Ordering::Relaxed);
        X_RESETS_IN_RUNS.fetch_add(resets, Ordering::Relaxed);
        X_TABLES_CHECKED.fetch_add(checked, Ordering::Relaxed);
    }
}

/// The runs `runs` searched by `workers` tasks of the current rayon pool (one `rayon::scope`),
/// each with one scratch made by `new_scratch` when it takes its first run (two in shadow mode),
/// into `subject_results` (one slot per subject, every slot of a run written by its task).
#[allow(clippy::too_many_arguments)]
fn x_search_runs<S, K>(
    workers: usize,
    runs: &[Range<usize>],
    subjects: &[EncodedProtein],
    subject_diag_offsets: &[i32],
    window: i32,
    new_scratch: &S,
    process_subject: &K,
    shadow: bool,
    subject_results: &mut [Option<BlastpSubjectResult>],
) where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    assert_eq!(subject_results.len(), subjects.len());
    // The slots of each run, taken by the task that takes the run.
    let mut rest = subject_results;
    let mut run_slots = Vec::with_capacity(runs.len());
    for run in runs {
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(run.len());
        run_slots.push(Mutex::new(head));
        rest = tail;
    }
    debug_assert!(rest.is_empty());
    let next_run = AtomicUsize::new(0);
    // No NCBI counterpart: scheduling of the preliminary loop over the subjects (NCBI's threads
    // take subjects from one iterator, blast_engine.c above); it changes no value NCBI computes.
    rayon::scope(|scope| {
        for _ in 0..workers.clamp(1, runs.len().max(1)) {
            scope.spawn(|_| {
                let mut scratch: Option<BlastpSubjectScratch> = None;
                let mut reference: Option<BlastpSubjectScratch> = None;
                loop {
                    let index = next_run.fetch_add(1, Ordering::Relaxed);
                    let Some(run) = runs.get(index) else {
                        break;
                    };
                    let scratch = scratch.get_or_insert_with(new_scratch);
                    let reference = shadow.then(|| reference.get_or_insert_with(new_scratch));
                    let mut slots = run_slots[index]
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner);
                    x_search_run(
                        scratch,
                        reference,
                        run.clone(),
                        &mut slots,
                        subjects,
                        subject_diag_offsets,
                        window,
                        process_subject,
                    );
                }
            });
        }
    });
}

// No NCBI counterpart: scheduling of the preliminary loop over the subjects (the cut into runs and
// their tasks); it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPSUBJRUN): the parallel loop of a `-subject` search in runs of
/// subjects (module header) on `pool`, into `subject_results`.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_subject_runs<S, K>(
    pool: &SearchPool<'_>,
    subjects: &[EncodedProtein],
    subject_diag_offsets: &[i32],
    window: i32,
    new_scratch: &S,
    process_subject: &K,
    shadow: bool,
    subject_results: &mut [Option<BlastpSubjectResult>],
) where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let workers = pool.threads().max(1);
    let lengths: Vec<usize> = subjects.iter().map(|subject| subject.aa_len).collect();
    let runs = x_cut_runs(&lengths, workers.saturating_mul(X_PARTS_PER_THREAD));
    if shadow {
        X_SEARCHES.fetch_add(1, Ordering::Relaxed);
    }
    pool.install(|| {
        x_search_runs(
            workers,
            &runs,
            subjects,
            subject_diag_offsets,
            window,
            new_scratch,
            process_subject,
            shadow,
            subject_results,
        )
    });
}

#[cfg(test)]
mod tests {
    use super::super::blastp_subject_diag_offsets;
    use super::*;

    const WINDOW: i32 = 40;
    const TABLE: i32 = 64;

    /// Lengths from a fixed linear congruential sequence: some subjects without letters, the
    /// others of 1-599 letters.
    fn lengths(n: usize, seed: u64) -> Vec<usize> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let value = (state >> 33) as usize;
                if value % 11 == 0 {
                    0
                } else {
                    1 + value % 599
                }
            })
            .collect()
    }

    fn proteins(lengths: &[usize]) -> Vec<EncodedProtein> {
        lengths
            .iter()
            .map(|&aa_len| EncodedProtein {
                // Only the length is read here; the letters are a sentinel pair.
                aa_seq: vec![0, 0],
                aa_len,
            })
            .collect()
    }

    #[test]
    fn runs_cover_every_subject_once_in_order_and_shrink_with_the_work_left() {
        let mut cases = vec![
            Vec::new(),
            vec![5],
            vec![0, 0, 0],
            vec![1_000_000, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 1_000_000],
        ];
        for seed in 1..6 {
            cases.push(lengths(97, seed));
        }
        for lengths in &cases {
            let weight = |range: Range<usize>| -> u64 {
                lengths[range].iter().map(|&length| length as u64 + 1).sum()
            };
            for parts in [0, 1, 2, 3, 7, 16, 64, 97, 1000] {
                let runs = x_cut_runs(lengths, parts);
                let what = format!("{lengths:?} {parts}: {runs:?}");
                let mut next = 0;
                for run in &runs {
                    assert_eq!(run.start, next, "{what}");
                    assert!(run.end > run.start, "{what}");
                    // The shortest run that reaches its share of the weight left.
                    let wanted = weight(run.start..lengths.len()).div_ceil(parts.max(1) as u64);
                    assert!(weight(run.clone()) >= wanted, "{what}");
                    assert!(
                        run.len() == 1 || weight(run.start..run.end - 1) < wanted,
                        "{what}"
                    );
                    next = run.end;
                }
                assert_eq!(next, lengths.len(), "{what}");
                if parts <= 1 && !lengths.is_empty() {
                    assert_eq!(runs.len(), 1, "{what}");
                }
            }
        }
        // Equal weights (4 each, 48 in all) in 4 parts: 12, then 9, 6, 4, 3, 2, 1 of the weight left.
        assert_eq!(
            x_cut_runs(&[3; 12], 4),
            vec![0..3, 3..6, 6..8, 8..9, 9..10, 10..11, 11..12]
        );
    }

    #[test]
    fn an_entry_reads_as_cleared_only_below_the_offset() {
        let offset = 1_000;
        let mut table = vec![DiagStruct::default(); 8];
        assert_eq!(x_first_uncleared(&table, WINDOW, WINDOW), None);
        assert_eq!(x_first_uncleared(&table, offset, WINDOW), None);
        let cleared = vec![DiagStruct::clear(WINDOW); 8];
        assert_eq!(x_first_uncleared(&cleared, WINDOW, WINDOW), None);
        table[3].set_last_hit(offset - WINDOW);
        assert_eq!(x_first_uncleared(&table, offset, WINDOW), None);
        table[3].set_last_hit(offset - WINDOW + 1);
        assert_eq!(x_first_uncleared(&table, offset, WINDOW), Some(3));
        // One hit (`window` 0): an entry up to the offset reads as cleared.
        assert_eq!(x_first_uncleared(&table, offset, 0), None);
        table[3].set_flag(1);
        table[3].set_last_hit(offset);
        assert_eq!(x_first_uncleared(&table, offset, WINDOW), None);
        table[3].set_last_hit(offset + 1);
        assert_eq!(x_first_uncleared(&table, offset, WINDOW), Some(3));
    }

    /// What the test kernel saw at the start of a subject's search: the subject, the table's
    /// offset, and the first entry that did not read as cleared.
    type Seen = (usize, i32, Option<usize>);

    /// How the test kernel ends a subject.
    #[derive(Clone, Copy, PartialEq)]
    enum Kernel {
        /// As `process_subject`: entries written below the next offset, then `finish_subject`.
        Correct,
        /// The offset is not advanced.
        NoFinish,
        /// At INT4_MAX / 4 the offset goes back to `window` but the table is not cleared.
        MissedReset,
        /// An entry is left where the next subject's test reads it (`last_hit` above that
        /// subject's offset less `window`).
        EntryTooHigh,
    }

    /// The subjects searched in runs by `threads` tasks with a kernel that writes entries as
    /// `process_subject` writes them (or as `kernel` breaks that), records what it saw, and
    /// moves the offset with `finish_subject`.
    fn search(
        lengths: &[usize],
        threads: usize,
        parts: usize,
        shadow: bool,
        kernel: Kernel,
    ) -> (Vec<Seen>, Vec<Option<BlastpSubjectResult>>) {
        let subjects = proteins(lengths);
        let offsets = blastp_subject_diag_offsets(&subjects, WINDOW);
        let runs = x_cut_runs(lengths, parts);
        let seen = Mutex::new(Vec::new());
        let new_scratch = || BlastpSubjectScratch::new(16, TABLE, WINDOW, 100);
        let process_subject = |scratch: &mut BlastpSubjectScratch,
                               s_idx: usize,
                               subject: BlastpSubjectRef<'_>,
                               result: &mut BlastpSubjectResult| {
            result.reset(s_idx);
            if subject.aa_len == 0 {
                return;
            }
            let offset = scratch.diag_offset;
            seen.lock().unwrap().push((
                s_idx,
                offset,
                x_first_uncleared(&scratch.diag_array, offset, WINDOW),
            ));
            let length = subject.aa_len as i32;
            for k in 0..5 {
                let entry = &mut scratch.diag_array[(s_idx * 7 + k) % TABLE as usize];
                entry.set_last_hit(offset + (k as i32 * 13) % length);
                entry.set_flag((k % 2) as u32);
            }
            match kernel {
                Kernel::Correct => scratch.finish_subject(subject.aa_len, WINDOW),
                Kernel::NoFinish => {}
                Kernel::MissedReset => {
                    if scratch.diag_offset >= i32::MAX / 4 {
                        scratch.diag_offset = WINDOW;
                    } else {
                        scratch.finish_subject(subject.aa_len, WINDOW);
                    }
                }
                Kernel::EntryTooHigh => {
                    // One above the highest `last_hit` the next subject reads as cleared.
                    scratch.diag_array[s_idx % TABLE as usize].set_last_hit(offset + length + 1);
                    scratch.finish_subject(subject.aa_len, WINDOW);
                }
            }
        };
        let mut results: Vec<Option<BlastpSubjectResult>> = std::iter::repeat_with(|| None)
            .take(lengths.len())
            .collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            x_search_runs(
                threads,
                &runs,
                &subjects,
                &offsets,
                WINDOW,
                &new_scratch,
                &process_subject,
                shadow,
                &mut results,
            )
        });
        let mut seen = seen.into_inner().unwrap();
        seen.sort_unstable();
        (seen, results)
    }

    /// Every subject with letters searched at the one-thread order's offset from a table that
    /// reads as cleared, `times` times, and every slot filled with its subject's result.
    fn assert_one_thread_offsets(
        lengths: &[usize],
        seen: &[Seen],
        results: &[Option<BlastpSubjectResult>],
        times: usize,
        what: &str,
    ) {
        let offsets = blastp_subject_diag_offsets(&proteins(lengths), WINDOW);
        let expected: Vec<Seen> = lengths
            .iter()
            .enumerate()
            .filter(|&(_, &length)| length > 0)
            .flat_map(|(s_idx, _)| std::iter::repeat_n((s_idx, offsets[s_idx], None), times))
            .collect();
        assert_eq!(seen, expected.as_slice(), "{what}");
        for (s_idx, result) in results.iter().enumerate() {
            assert_eq!(
                result.as_ref().map(|result| result.subject_index),
                Some(s_idx),
                "{what}"
            );
        }
    }

    #[test]
    fn runs_search_each_subject_at_its_one_thread_offset() {
        for seed in 1..4 {
            let lengths = lengths(150, seed);
            for threads in [1, 2, 3, 7] {
                for parts in [1, 2, 5, threads * 4, 150, 400] {
                    for shadow in [false, true] {
                        let (seen, results) =
                            search(&lengths, threads, parts, shadow, Kernel::Correct);
                        let what = format!("seed {seed} threads {threads} parts {parts}");
                        assert_one_thread_offsets(
                            &lengths,
                            &seen,
                            &results,
                            1 + usize::from(shadow),
                            &what,
                        );
                    }
                }
            }
        }
    }

    // Subjects long enough for the offset to reach INT4_MAX / 4 several times: the reset
    // (`finish_subject`, `Blast_ExtendWordExit`) happens inside the runs, clears the table, and the
    // next subject starts at `window` as in the one-thread order.
    #[test]
    fn the_reset_at_int4_max_quarter_happens_inside_a_run() {
        let mut lengths = vec![200_000_000; 12];
        lengths[5] = 0;
        lengths[6] = 3;
        let offsets = blastp_subject_diag_offsets(&proteins(&lengths), WINDOW);
        let resets = offsets.windows(2).filter(|pair| pair[1] < pair[0]).count();
        assert!(resets >= 2, "{offsets:?}");
        for threads in [1, 2, 3] {
            for parts in [1, 2, 3] {
                for shadow in [false, true] {
                    let (seen, results) = search(&lengths, threads, parts, shadow, Kernel::Correct);
                    let what = format!("threads {threads} parts {parts} shadow {shadow}");
                    assert_one_thread_offsets(
                        &lengths,
                        &seen,
                        &results,
                        1 + usize::from(shadow),
                        &what,
                    );
                }
            }
        }
    }

    /// The panic message of a shadow search with `kernel`, which must stop it.
    fn shadow_stops(lengths: &[usize], kernel: Kernel) -> String {
        let lengths = lengths.to_vec();
        let panic = std::panic::catch_unwind(move || {
            search(&lengths, 2, 1, true, kernel);
        })
        .expect_err("the shadow must stop the run");
        panic
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| panic.downcast_ref::<&str>().map(|text| text.to_string()))
            .unwrap_or_default()
    }

    // The shadow's checks catch a run that searches a subject at another offset, a reset that does
    // not clear the table, and an entry left where the next subject reads it.
    #[test]
    fn the_shadow_stops_a_wrong_offset_a_missed_reset_and_a_stale_entry() {
        let short = [50, 60, 70, 80];
        let message = shadow_stops(&short, Kernel::NoFinish);
        assert!(
            message.contains("[X_SHADOW] BPSUBJRUN") && message.contains("offset"),
            "{message}"
        );
        let long = [200_000_000, 200_000_000, 200_000_000, 200_000_000, 50];
        let message = shadow_stops(&long, Kernel::MissedReset);
        assert!(message.contains("does not read as cleared"), "{message}");
        let message = shadow_stops(&short, Kernel::EntryTooHigh);
        assert!(message.contains("does not read as cleared"), "{message}");
    }
}
