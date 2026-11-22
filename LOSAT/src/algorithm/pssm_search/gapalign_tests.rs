use super::*;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::blast_util::BLASTAA_SIZE;

/// A PSSM built from BLOSUM62's rows of `query` (NCBIstdaa): the PSSM kernels must then give
/// the alignments of blastp's matrix kernels.
fn blosum62_pssm(query: &[u8]) -> QueryPssm {
    QueryPssm {
        rows: query
            .iter()
            .map(|&q| {
                let mut row = [0i32; BLASTAA_SIZE];
                for (s, cell) in row.iter_mut().enumerate() {
                    *cell = crate::utils::matrix::blosum62_score(q, s as u8);
                }
                row
            })
            .collect(),
    }
}

fn encode(text: &str) -> Vec<u8> {
    text.bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect()
}

// The PSSM kernels on BLOSUM62 rows agree with blastp's matrix kernels (gapalign.rs), which
// the BLASTP fixtures check against NCBI: same ends, score and edit script, on a pair with a
// gap in each sequence.
#[test]
fn blosum62_rows_give_blastps_alignments() {
    // Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
    let query = encode("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ");
    let subject_letters =
        encode("MKTAYIAKQRQISFVKSHFSRQDLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVKVKVKALPDAQ");
    let mut subject = subject_letters.clone();
    subject.push(0);
    let rows = blosum62_pssm(&query);
    let mut scratch = PssmGapScratch::new();
    let (q_off, s_off) = (30usize, 31usize);
    let ours = pssm_score_only_gapped_alignment(
        &rows,
        &subject,
        subject_letters.len(),
        q_off,
        s_off,
        11,
        1,
        38,
        &mut scratch,
    );
    let theirs = crate::algorithm::blastp::gapalign::blastp_score_only_gapped_alignment(
        &query,
        &subject_letters,
        q_off,
        s_off,
        crate::config::ScoringMatrix::Blosum62,
        11,
        1,
        38,
        crate::algorithm::blastp::gapalign::BlastpGappedAlignmentMode::Exact,
    );
    assert_eq!(
        (
            ours.query_start,
            ours.query_stop,
            ours.subject_start,
            ours.subject_stop,
            ours.score
        ),
        (
            theirs.query_start,
            theirs.query_stop,
            theirs.subject_start,
            theirs.subject_stop,
            theirs.score
        )
    );
    let traced = pssm_gapped_alignment_with_traceback(
        &rows,
        &subject,
        subject_letters.len(),
        q_off,
        s_off,
        11,
        1,
        64,
        &mut scratch,
    );
    let reference = crate::algorithm::blastp::gapalign::blast_gapped_alignment_with_traceback(
        &query,
        &subject_letters,
        q_off,
        s_off,
        crate::config::ScoringMatrix::Blosum62,
        None,
        11,
        1,
        64,
    )
    .expect("an alignment");
    assert_eq!(traced.query_start, reference.query_start);
    assert_eq!(traced.query_stop, reference.query_stop);
    assert_eq!(traced.subject_start, reference.subject_start);
    assert_eq!(traced.subject_stop, reference.subject_stop);
    assert_eq!(traced.score, reference.score);
    assert_eq!(traced.edit_script, reference.edit_script);
    assert!(traced
        .edit_script
        .iter()
        .any(|op| matches!(op, GapEditOp::Del(_))));
}

// NCBI's `ALIGN_EX` (cited at `semi_gapped_align_score_only`) raises an X-drop below the
// gap costs to `gap_open_extend`: X-drops of 0 and below (RPS-BLAST's `-xdrop_gap 0`, `-1`)
// align as the cost itself does, in `Int4` arithmetic.
#[test]
fn x_drops_at_or_below_zero_align_as_the_gap_costs() {
    let query = encode("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ");
    let subject_letters =
        encode("MKTAYIAKQRQISFVKSHFSRQDLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVKVKVKALPDAQ");
    let mut subject = subject_letters.clone();
    subject.push(0);
    let rows = blosum62_pssm(&query);
    let mut scratch = PssmGapScratch::new();
    let (q_off, s_off) = (30usize, 31usize);
    let score_only = |x_drop: i32, scratch: &mut PssmGapScratch| {
        let found = pssm_score_only_gapped_alignment(
            &rows,
            &subject,
            subject_letters.len(),
            q_off,
            s_off,
            11,
            1,
            x_drop,
            scratch,
        );
        (
            found.query_start,
            found.query_stop,
            found.subject_start,
            found.subject_stop,
            found.score,
        )
    };
    let traced = |x_drop: i32, scratch: &mut PssmGapScratch| {
        let found = pssm_gapped_alignment_with_traceback(
            &rows,
            &subject,
            subject_letters.len(),
            q_off,
            s_off,
            11,
            1,
            x_drop,
            scratch,
        );
        (
            found.query_start,
            found.query_stop,
            found.subject_start,
            found.subject_stop,
            found.score,
            found.edit_script,
        )
    };
    let at_cost = score_only(12, &mut scratch);
    let traced_at_cost = traced(12, &mut scratch);
    for x_drop in [0, -1, -200, i32::MIN] {
        assert_eq!(score_only(x_drop, &mut scratch), at_cost, "{x_drop}");
        assert_eq!(traced(x_drop, &mut scratch), traced_at_cost, "{x_drop}");
    }
    // The cost is a real X-drop: a larger one can find more.
    assert!(score_only(38, &mut scratch).4 >= at_cost.4);
}

#[test]
fn gapped_start_is_the_middle_of_the_best_window() {
    let query = encode("AAAAAAAAAAWWWWWWWWWWWWAAAA");
    let rows = blosum62_pssm(&query);
    // The subject equals the query: the windows grow to the first all-W one (ending at 20).
    let start = pssm_get_start_for_gapped_alignment(&rows, &query, 0, 26, 0);
    assert_eq!(start, 20 - HSP_MAX_WINDOW / 2);
    // A short HSP starts at its middle.
    assert_eq!(
        pssm_get_start_for_gapped_alignment(&rows, &query, 4, 10, 4),
        9
    );
}

/// The score-only DP before D5 of PS12s (`reverse` a value, each column's letter read by
/// `subject_letter`): the reference that `semi_gapped_align_score_only` must equal.
#[allow(clippy::too_many_arguments)]
fn score_only_reference<R: PssmRows>(
    rows: &R,
    m: usize,
    query_offset: usize,
    subject: &[u8],
    s_base: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
) -> (usize, usize, i32) {
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = x_drop.max(gap_open_extend);
    if n == 0 || m == 0 {
        return (0, 0, 0);
    }
    let mut dp = vec![
        GapDp {
            best: MININT,
            best_gap: MININT,
        };
        n + 2
    ];
    let mut score = -gap_open_extend;
    dp[0].best = 0;
    dp[0].best_gap = -gap_open_extend;
    let mut i = 1usize;
    while i <= n {
        if score < -x_dropoff {
            break;
        }
        dp[i].best = score;
        dp[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        i += 1;
    }
    let mut b_size = i;
    let mut best_score = 0i32;
    let mut a_offset = 0usize;
    let mut b_offset = 0usize;
    let mut first_b_index = 0usize;
    for a_index in 1..=m {
        let matrix_row = pssm_row(rows, m, a_index, query_offset, reverse);
        let mut score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;
        while b_index < b_size {
            let letter = subject_letter(subject, s_base, n, b_index, reverse);
            let mut score_gap_col = dp[b_index].best_gap;
            let next_score = dp[b_index].best + matrix_row[letter as usize];
            if score < score_gap_col {
                score = score_gap_col;
            }
            if score < score_gap_row {
                score = score_gap_row;
            }
            if best_score - score > x_dropoff {
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    dp[b_index].best = MININT;
                }
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score = score;
                    a_offset = a_index;
                    b_offset = b_index;
                }
                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                dp[b_index].best_gap = (score - gap_open_extend).max(score_gap_col);
                score_gap_row = (score - gap_open_extend).max(score_gap_row);
                dp[b_index].best = score;
            }
            score = next_score;
            b_index += 1;
        }
        if first_b_index == b_size {
            break;
        }
        if last_b_index + 1 < b_size {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= (best_score - x_dropoff) && b_size <= n {
                dp[b_size].best = score_gap_row;
                dp[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n {
            dp[b_size].best = MININT;
            dp[b_size].best_gap = MININT;
            b_size += 1;
        }
    }
    (a_offset, b_offset, best_score)
}

// D5 of PS12s: the score-only DP with `reverse` a constant and the subject's letters taken
// once per call gives the reference's ends and score, on random PSSMs and subjects, both
// directions, scratch kept from call to call, and calls whose columns read the sentinel or
// past the subject (letter 0): extensions that reach the subject's end, a subject cut short,
// a base past it.
#[test]
fn the_score_only_dp_equals_its_reference() {
    let mut state = 7u64;
    let mut next = |bound: u64| -> u64 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % bound
    };
    let mut scratch = PssmGapScratch::new();
    for case in 0..3000 {
        let length = 1 + next(120) as usize;
        let rows = QueryPssm {
            rows: (0..length)
                .map(|_| {
                    let mut row = [0i32; BLASTAA_SIZE];
                    for cell in row.iter_mut() {
                        *cell = next(16) as i32 - 9;
                    }
                    // Matches often score well, so that alignments grow.
                    row[(next(20) + 1) as usize] = 5 + next(6) as i32;
                    row
                })
                .collect(),
        };
        let subject_length = next(150) as usize;
        let mut subject: Vec<u8> = (0..subject_length).map(|_| 1 + next(24) as u8).collect();
        if next(4) != 0 {
            subject.push(0);
        }
        let m = 1 + next(length as u64) as usize;
        let query_offset = next((length - m + 1) as u64) as usize;
        let s_base = next(subject.len() as u64 + 3) as usize;
        let n = next(subject.len() as u64 + 4) as usize;
        let (gap_open, gap_extend) = [(11, 1), (1100, 100), (0, 1), (9, 2)][next(4) as usize];
        let x_drop = [0, 7, 38, 160, 3000][next(5) as usize];
        for reverse in [false, true] {
            // The reverse DP reads rows `m - a_index` (0 to m - 1), the forward one rows
            // `a_index + query_offset` (up to m + query_offset): both inside `length`.
            if !reverse && m + query_offset >= length {
                continue;
            }
            let expected = score_only_reference(
                &rows,
                m,
                query_offset,
                &subject,
                s_base,
                n,
                gap_open,
                gap_extend,
                x_drop,
                reverse,
            );
            let got = if reverse {
                semi_gapped_align_score_only::<_, true>(
                    &rows,
                    m,
                    query_offset,
                    &subject,
                    s_base,
                    n,
                    gap_open,
                    gap_extend,
                    x_drop,
                    &mut scratch,
                )
            } else {
                semi_gapped_align_score_only::<_, false>(
                    &rows,
                    m,
                    query_offset,
                    &subject,
                    s_base,
                    n,
                    gap_open,
                    gap_extend,
                    x_drop,
                    &mut scratch,
                )
            };
            assert_eq!(got, expected, "case {case} reverse {reverse}");
        }
    }
}
