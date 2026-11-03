//! EXPERIMENT (LOSAT_X_BPSUBJRUN): the parallel preliminary loop of a `-subject` search in runs of
//! consecutive subjects, one diagonal table per task, filled only where it would not read as
//! cleared.
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
//! the runs, in order, from one counter. A run starts at the diagonal offset that the one-thread
//! order gives its first subject (`blastp_subject_diag_offsets`), from a table that reads as
//! cleared there (`x_start_run`): the task fills its table (`prepare_independent_subject`) only
//! when that offset is below the one its table reads as cleared from, and otherwise sets the offset
//! and runs the other per-subject resets (below). The later subjects of the run are searched as the
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
//! the test of each word hit (`aa_ungapped.c`, below; the fused loops make
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
//! Why a run needs no fill. The same holds for a task's table from one of its runs to the next:
//! a new scratch's table (every entry 0 with flag 0, at the offset `window`) reads as cleared at
//! every offset from `window` on, and after each `finish_subject` every entry reads as cleared at
//! every offset from the one it leaves on (the bounds above only get looser as the offset grows;
//! after the reset the entries of `s_BlastDiagClear` read so from `window` on). Only the task's own
//! subjects write its table, so the subjects other tasks search in between, and the resets of the
//! one-thread order among them, do not matter: a run whose first offset is at least the table's
//! offset starts from a table that reads as cleared, as a filled one does, and only a run whose
//! first offset is below it needs the fill. A task takes its runs in increasing order, so that
//! happens only after the one-thread order's offset has gone back to `window` since the task's
//! previous run, and never in a search whose offsets stay below INT4_MAX / 4. A run also fills
//! the table when the table's offset has reached 2^30, from where an entry's `last_hit` may have
//! wrapped in its 31 bits (below). Until 2026-11-03 every run filled the table (Species4 x Species6
//! blastp-fast at 16 threads: 68,627 fills, 17.7 GB of writes in all; structural line S1).
//!
//! The argument is in integers; `last_hit` is a 31-bit field (`DiagStruct`), whose stored values
//! are exact while `d + L < 2^30`. An offset below INT4_MAX / 4 moves by one subject and one
//! window, so a value can wrap only where consecutive subjects have about 2^29 letters or more
//! together. A later subject of the same run reads such an entry as the one-thread loop and NCBI
//! do (it may read differently from a cleared one), and it is the loop without the switch that
//! departs from them. Every entry is at most the offset that `finish_subject` left, so none has
//! wrapped while the table's offset is below 2^30; from there a run starts from a filled table, as
//! every run did before 2026-11-03 (`x_start_run`), so which task takes the run does not matter.
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
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

use crate::utils::threading::SearchPool;

use super::{BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, EncodedProtein};

// No NCBI counterpart: reads the LOSAT_X_BPSUBJRUN switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSUBJRUN): a strict-set switch (on unless it is `0`).
pub(super) fn x_bp_subj_run() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPSUBJRUN"))
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
// No NCBI counterpart: when a run's start may skip the fill (`x_start_run`); it changes no value
// NCBI computes.
/// The first offset at which a table's entries may hold `last_hit` values that wrapped in the
/// 31-bit field (`DiagStruct`, 2^30): every entry is at most the offset `finish_subject` left, so
/// below it none has (module header).
const X_EXACT_LAST_HIT: i32 = 1 << 30;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-103
// (`s_BlastDiagClear`, module header)
/// The start of a run whose first subject's offset in the one-thread order is `diag_offset`: the
/// task's table at that offset, reading as cleared there. `scratch.diag_offset` is the offset from
/// which every entry of the table reads as cleared (module header, "Why a run needs no fill"), so
/// the table is filled (`prepare_independent_subject`) only when `diag_offset` is below it, or when
/// the table's offset has reached `X_EXACT_LAST_HIT` (an entry may have wrapped); otherwise the
/// offset is set and the other per-subject resets run, as before every later subject of a run.
fn x_start_run(scratch: &mut BlastpSubjectScratch, diag_offset: i32) {
    if diag_offset < scratch.diag_offset || scratch.diag_offset >= X_EXACT_LAST_HIT {
        scratch.prepare_independent_subject(diag_offset);
    } else {
        scratch.diag_offset = diag_offset;
        scratch.reset_subject_state();
    }
}

/// The subjects of one run, `run`, searched with `scratch` into `slots` (the run's slots of the
/// subject results).
fn x_search_run<K>(
    scratch: &mut BlastpSubjectScratch,
    run: Range<usize>,
    slots: &mut [Option<BlastpSubjectResult>],
    subjects: &[EncodedProtein],
    subject_diag_offsets: &[i32],
    process_subject: &K,
) where
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    for (s_idx, slot) in run.clone().zip(slots.iter_mut()) {
        let diag_offset = subject_diag_offsets[s_idx];
        if s_idx == run.start {
            x_start_run(scratch, diag_offset);
        } else {
            // The previous subject's `finish_subject` (`Blast_ExtendWordExit`) left the
            // one-thread order's offset, and its entries read as cleared ones (module header).
            debug_assert_eq!(scratch.diag_offset, diag_offset, "subject {s_idx}");
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
        *slot = Some(subject_result);
    }
}

/// The runs `runs` searched by `workers` tasks of the current rayon pool (one `rayon::scope`),
/// each with one scratch made by `new_scratch` when it takes its first run, into `subject_results`
/// (one slot per subject, every slot of a run written by its task).
fn x_search_runs<S, K>(
    workers: usize,
    runs: &[Range<usize>],
    subjects: &[EncodedProtein],
    subject_diag_offsets: &[i32],
    new_scratch: &S,
    process_subject: &K,
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
                loop {
                    let index = next_run.fetch_add(1, Ordering::Relaxed);
                    let Some(run) = runs.get(index) else {
                        break;
                    };
                    let scratch = scratch.get_or_insert_with(new_scratch);
                    let mut slots = run_slots[index]
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner);
                    x_search_run(
                        scratch,
                        run.clone(),
                        &mut slots,
                        subjects,
                        subject_diag_offsets,
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
pub(super) fn x_search_subject_runs<S, K>(
    pool: &SearchPool<'_>,
    subjects: &[EncodedProtein],
    subject_diag_offsets: &[i32],
    new_scratch: &S,
    process_subject: &K,
    subject_results: &mut [Option<BlastpSubjectResult>],
) where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let workers = pool.threads().max(1);
    let lengths: Vec<usize> = subjects.iter().map(|subject| subject.aa_len).collect();
    let runs = x_cut_runs(&lengths, workers.saturating_mul(X_PARTS_PER_THREAD));
    pool.install(|| {
        x_search_runs(
            workers,
            &runs,
            subjects,
            subject_diag_offsets,
            new_scratch,
            process_subject,
            subject_results,
        )
    });
}

#[cfg(test)]
mod tests {
    use super::super::blastp_subject_diag_offsets;
    use super::*;
    use crate::algorithm::tblastx::blast_extend::DiagStruct;

    const WINDOW: i32 = 40;
    const TABLE: i32 = 64;

    // No NCBI counterpart: a check of the unit kernel (the property in the module's header); it
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

    // A run keeps the task's table at or above the offset it reads as cleared from, and fills it
    // below that offset (`x_start_run`).
    #[test]
    fn a_run_fills_the_table_only_below_its_offset() {
        let mut scratch = BlastpSubjectScratch::new(16, TABLE, WINDOW, 100);
        let fresh = scratch.diag_array.clone();
        // A new table reads as cleared from `window` on.
        x_start_run(&mut scratch, WINDOW + 500);
        assert_eq!(scratch.diag_offset, WINDOW + 500);
        // A subject of 40 letters writes entries below its end.
        let d = scratch.diag_offset;
        scratch.diag_array[3].set_last_hit(d + 39);
        scratch.diag_array[3].set_flag(1);
        scratch.diag_array[5].set_last_hit(d + 10);
        scratch.finish_subject(40, WINDOW);
        let left = scratch.diag_offset;
        assert_eq!(left, d + 40 + WINDOW);
        let written = scratch.diag_array.clone();
        // A run at or beyond that offset keeps the entries, which read as cleared there.
        for start in [left, left + 1_000] {
            x_start_run(&mut scratch, start);
            assert_eq!(scratch.diag_offset, start);
            assert_eq!(scratch.diag_array, written);
            assert_eq!(x_first_uncleared(&scratch.diag_array, start, WINDOW), None);
        }
        // Below it (the one-thread order went back to `window`) the entries need not read as
        // cleared, so the run fills the table.
        assert_eq!(x_first_uncleared(&written, WINDOW, WINDOW), Some(3));
        x_start_run(&mut scratch, WINDOW);
        assert_eq!(scratch.diag_offset, WINDOW);
        assert_eq!(scratch.diag_array, fresh);
        // From 2^30 on an entry may have wrapped in its 31 bits: the run fills the table at any
        // offset.
        scratch.diag_array[3].set_last_hit(5);
        scratch.diag_offset = X_EXACT_LAST_HIT;
        x_start_run(&mut scratch, X_EXACT_LAST_HIT + 100);
        assert_eq!(scratch.diag_offset, X_EXACT_LAST_HIT + 100);
        assert_eq!(scratch.diag_array, fresh);
        scratch.diag_array[3].set_last_hit(5);
        scratch.diag_offset = X_EXACT_LAST_HIT - 1;
        x_start_run(&mut scratch, X_EXACT_LAST_HIT + 100);
        assert_eq!(scratch.diag_array[3].last_hit(), 5);
    }

    /// What the test kernel saw at the start of a subject's search: the subject, the table's
    /// offset, and the first entry that did not read as cleared.
    type Seen = (usize, i32, Option<usize>);

    /// The subjects searched in runs by `threads` tasks with a kernel that writes entries as
    /// `process_subject` writes them, records what it saw, and moves the offset with
    /// `finish_subject`.
    fn search(
        lengths: &[usize],
        threads: usize,
        parts: usize,
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
            scratch.finish_subject(subject.aa_len, WINDOW);
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
                &new_scratch,
                &process_subject,
                &mut results,
            )
        });
        let mut seen = seen.into_inner().unwrap();
        seen.sort_unstable();
        (seen, results)
    }

    /// Every subject with letters searched at the one-thread order's offset from a table that
    /// reads as cleared, and every slot filled with its subject's result.
    fn assert_one_thread_offsets(
        lengths: &[usize],
        seen: &[Seen],
        results: &[Option<BlastpSubjectResult>],
        what: &str,
    ) {
        let offsets = blastp_subject_diag_offsets(&proteins(lengths), WINDOW);
        let expected: Vec<Seen> = lengths
            .iter()
            .enumerate()
            .filter(|&(_, &length)| length > 0)
            .map(|(s_idx, _)| (s_idx, offsets[s_idx], None))
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
                    let (seen, results) = search(&lengths, threads, parts);
                    let what = format!("seed {seed} threads {threads} parts {parts}");
                    assert_one_thread_offsets(&lengths, &seen, &results, &what);
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
                let (seen, results) = search(&lengths, threads, parts);
                let what = format!("threads {threads} parts {parts}");
                assert_one_thread_offsets(&lengths, &seen, &results, &what);
            }
        }
    }

    // One task's table across runs that skip the subjects of other tasks, with and without a
    // reset of the one-thread order among the skipped ones, and runs that start without a fill
    // after a reset (by another task's subject, or by the task's own last one) (the order of
    // `x_search_runs` when other tasks take the runs in between): every subject still starts at
    // its one-thread offset from a table that reads as cleared.
    #[test]
    fn a_task_skipping_runs_starts_each_run_from_a_cleared_table() {
        let seen = Mutex::new(Vec::new());
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
            for k in 0..TABLE as usize {
                let entry = &mut scratch.diag_array[k];
                entry.set_last_hit(offset + (k as i32 * 13) % length);
                entry.set_flag((k % 2) as u32);
            }
            scratch.finish_subject(subject.aa_len, WINDOW);
        };
        // Short subjects (no reset), and subjects long enough to reset between the task's runs.
        let short: Vec<usize> = (0..12).map(|i| 30 + 7 * i).collect();
        let long = vec![200_000_000; 12];
        let cases: [(&Vec<usize>, Vec<Range<usize>>, bool); 4] = [
            (&short, vec![0..2, 5..7, 10..12], false),
            (&long, vec![0..2, 4..6, 9..11], true),
            (&long, vec![0..1, 6..7], false),
            (&long, vec![3..4, 6..7], true),
        ];
        for (lengths, task_runs, goes_back) in cases {
            let subjects = proteins(lengths);
            let offsets = blastp_subject_diag_offsets(&subjects, WINDOW);
            let back = task_runs
                .windows(2)
                .filter(|pair| offsets[pair[1].start] < offsets[pair[0].end - 1])
                .count();
            assert_eq!(back >= 1, goes_back, "{task_runs:?} {offsets:?}");
            let mut scratch = BlastpSubjectScratch::new(16, TABLE, WINDOW, 100);
            seen.lock().unwrap().clear();
            for run in &task_runs {
                let mut slots: Vec<Option<BlastpSubjectResult>> =
                    std::iter::repeat_with(|| None).take(run.len()).collect();
                x_search_run(
                    &mut scratch,
                    run.clone(),
                    &mut slots,
                    &subjects,
                    &offsets,
                    &process_subject,
                );
            }
            let expected: Vec<Seen> = task_runs
                .iter()
                .flat_map(|run| run.clone())
                .map(|s_idx| (s_idx, offsets[s_idx], None))
                .collect();
            assert_eq!(*seen.lock().unwrap(), expected, "{lengths:?}");
        }
    }
}
