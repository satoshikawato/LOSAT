//! EXPERIMENT (LOSAT_X_BPSUBJRUN): the parallel preliminary loop of a `-subject` search in runs of
//! consecutive subjects, one diagonal table fill per run.
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
//! The argument is in integers; `last_hit` is a 31-bit field (`DiagStruct`), whose stored values
//! are exact while `d + L < 2^30`. The offsets stay below INT4_MAX / 4 plus one subject and one
//! window, so this holds for every subject shorter than 2^29 letters. Only one subject of more than
//! 2^29 letters in one-hit mode could leave a wrapped entry that a later subject of its run reads
//! differently from a cleared one; there the run path keeps the leftovers as the one-thread loop
//! and NCBI do, and it is the loop without the switch that departs from them.
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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-103
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-544
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413,1553-1554

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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-172
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-103
            // (`s_BlastDiagClear`, module header): the run starts from a cleared table at the
            // one-thread order's offset of its first subject, as each subject does without the
            // switch.
            scratch.prepare_independent_subject(diag_offset);
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
#[path = "x_subject_runs_tests.rs"]
mod tests;
