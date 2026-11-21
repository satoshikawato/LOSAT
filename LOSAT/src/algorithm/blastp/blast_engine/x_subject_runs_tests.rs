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
            if value.is_multiple_of(11) {
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
