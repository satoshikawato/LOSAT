//! Hand-made matches for BLASTP's Smith-Waterman redo (`redo_one_match_smith_waterman`): the
//! window's alignments in NCBI's list order, the forbidden ranges between them, the heap's test
//! before the first alignment, the cutoff score of sum statistics and the regime of a significant
//! score of 0; and random matches for the redo resumed at a heap test left open.

use super::*;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{build_matrix_info, read_aa_composition};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, build_query_word_hashes, redo_calc_lambda, BlastCompoAdjustMode,
    BlastCompoAlignmentContext, BlastCompoGappingParams, BlastCompoSequenceData,
    BlastCompoSequenceRange, BlastRedoRangeResult, EMatrixAdjustRule,
};
use crate::core::composition_adjustment::smith_waterman_match::{SwRedoEndless, SwRedoFailed};
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::aa_char_to_ncbistdaa;

const SCALE: f64 = 32.0;

fn encode(text: &str) -> Vec<u8> {
    text.bytes().map(aa_char_to_ncbistdaa).collect()
}

/// BLASTP's redo parameters of composition mode 1 with BLOSUM62 11/1, scaled by 32 as
/// `RedoParamsBuilder::build` scales them.
fn params(do_link_hsps: bool, cutoff_score: i32) -> BlastRedoAlignParams {
    params_in(
        BlastCompoAdjustMode::CompositionBasedStats,
        do_link_hsps,
        cutoff_score,
    )
}

/// `params` in composition mode `mode`.
fn params_in(
    mode: BlastCompoAdjustMode,
    do_link_hsps: bool,
    cutoff_score: i32,
) -> BlastRedoAlignParams {
    BlastRedoAlignParams {
        matrix_info: build_matrix_info(ScoringMatrix::Blosum62, 0.3176 / SCALE).unwrap(),
        gapping_params: BlastCompoGappingParams {
            gap_open: (11.0 * SCALE) as i32,
            gap_extend: SCALE as i32,
            x_dropoff: 2000,
            context: std::cell::Cell::new(None),
        },
        compo_adjust_mode: mode,
        alphsize: crate::utils::matrix::BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: 0.041f64.glibc_ln(),
        score_divisor: SCALE,
        smith_waterman: true,
        near_identical_cutoff: 0.0,
        position_based: false,
        ccat_query_length: 0,
        query_is_translated: false,
        subject_is_translated: false,
        cutoff_score,
        do_link_hsps,
    }
}

/// The whole query and subject, unmasked (`s_SequenceGetProteinRange` without SEG).
fn whole_range<'seq>(
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    _query_words: Option<&[u64]>,
    _align: &BlastCompoAlignment,
    _should_test_identical: bool,
    _subject_maybe_biased: bool,
    _params: &BlastRedoAlignParams,
) -> Result<BlastRedoRangeResult> {
    let mut subject = BlastCompoSequenceData::from_ncbistdaa(matching_seq.data);
    subject.fit_protein_range_in_place(subject_range);
    Ok(BlastRedoRangeResult {
        query: BlastCompoSequenceData::copy_query_range_with_selenocysteine_fix(
            orig_query,
            query_range,
        ),
        subject,
        subject_maybe_biased: true,
        x_subject_composition: None,
    })
}

fn no_x_drop_redo(
    _: &BlastCompoAlignment,
    _: EMatrixAdjustRule,
    _: Option<&AdjustedProteinMatrix>,
    _: &BlastCompoSequenceData,
    _: &BlastCompoSequenceRange,
    _: i32,
    _: &BlastCompoSequenceData,
    _: &BlastCompoSequenceRange,
    _: i32,
    _: &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    panic!("the Smith-Waterman redo makes no X-drop redo of an HSP")
}

const CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range: whole_range,
    redo_one_alignment: no_x_drop_redo,
};

/// The incoming alignments of query 0: `(query_start, query_end, subject_start, subject_end)`.
fn incoming(hsps: &[(i32, i32, i32, i32)]) -> Option<Box<BlastCompoAlignment>> {
    let mut head = None;
    for &(qs, qe, ss, se) in hsps.iter().rev() {
        let mut align = blast_compo_alignment_new(
            ((qe - qs) as f64 * 4.0 * SCALE) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            qs,
            qe,
            0,
            ss,
            se,
            0,
            None,
        );
        align.next = head;
        head = Some(align);
    }
    head
}

fn query_info(query: &[u8], searchsp: f64) -> BlastCompoQueryInfo {
    BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(query),
        composition: read_aa_composition(query),
        eff_search_space: searchsp,
        words: Some(build_query_word_hashes(query)),
    }
}

struct Run {
    result: Result<BlastRedoOneMatchResult>,
    matrix: Option<AdjustedProteinMatrix>,
    heap_asks: usize,
}

fn run(
    query: &[u8],
    subject: &[u8],
    hsps: &[(i32, i32, i32, i32)],
    params: &BlastRedoAlignParams,
    searchsp: f64,
    cutoff_e: f64,
    heap_answer: bool,
) -> Run {
    run_from(
        query,
        subject,
        hsps,
        params,
        (searchsp, cutoff_e),
        None,
        &mut || SwHeapAnswer::Known(heap_answer),
    )
    .0
}

/// `run` with its heap test answered by `heap`, from `resume` where given; with the state that the
/// call keeps where `heap` leaves the test open.
fn run_from(
    query: &[u8],
    subject: &[u8],
    hsps: &[(i32, i32, i32, i32)],
    params: &BlastRedoAlignParams,
    (searchsp, cutoff_e): (f64, f64),
    resume: Option<&SwMatchResume>,
    heap: &mut dyn FnMut() -> SwHeapAnswer,
) -> (Run, Option<Box<SwMatchResume>>) {
    let info = query_info(query, searchsp);
    let matching_seq = BlastCompoMatchingSequence::new(7, subject);
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    let mut matrix = None;
    let mut heap_asks = 0;
    let mut scratch = SwMatchScratch::new();
    let mut kept = None;
    let result = redo_one_match_smith_waterman(
        &incoming(hsps),
        params,
        &matching_seq,
        std::slice::from_ref(&info),
        0.267 / SCALE,
        &CALLBACKS,
        &mut workspace,
        &mut matrix,
        cutoff_e,
        resume,
        &mut |query_index, _evalue, _score| {
            assert_eq!(query_index, 0);
            heap_asks += 1;
            heap()
        },
        &mut kept,
        &mut scratch,
    );
    (
        Run {
            result,
            matrix,
            heap_asks,
        },
        kept,
    )
}

/// The alignments of query 0, the list's head (the last one made) first.
fn alignments(result: &BlastRedoOneMatchResult) -> Vec<&BlastCompoAlignment> {
    let mut out = Vec::new();
    let mut current = result.alignments_by_query[0].as_deref();
    while let Some(align) = current {
        out.push(align);
        current = align.next.as_deref();
    }
    out
}

/// The best local score of the whole rectangle (Gotoh, NCBI's gap cost `open + k * extend`).
fn whole_matrix_best(rows: &[[i32; 28]], subject: &[u8], open: i32, extend: i32) -> i32 {
    let n = subject.len();
    let neg = i32::MIN / 4;
    let mut h_prev = vec![0i32; n + 1];
    let mut f_prev = vec![neg; n + 1];
    let mut best = 0;
    for row in rows {
        let mut h = vec![0i32; n + 1];
        let mut f = vec![neg; n + 1];
        let mut e = neg;
        for j in 1..=n {
            e = (e - extend).max(h[j - 1] - open - extend);
            f[j] = (f_prev[j] - extend).max(h_prev[j] - open - extend);
            let diagonal = h_prev[j - 1] + row[usize::from(subject[j - 1])];
            h[j] = 0.max(diagonal).max(e).max(f[j]);
            best = best.max(h[j]);
        }
        h_prev = h;
        f_prev = f;
    }
    best
}

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const QUERY: &str = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV";
const SPACER: &str = "PGSNPGSNPGSNPGSNPGSN";

// Two copies of the query's part in the subject and two HSPs: NCBI's loop makes the best local
// alignment, forbids its range and makes the second copy's; the list's head is the second.
#[test]
fn two_copies_give_two_alignments_in_ncbis_list_order() {
    let query = encode(QUERY);
    let subject_text = format!("{}{}{}", &QUERY[..40], SPACER, &QUERY[8..50]);
    let subject = encode(&subject_text);
    let p = params(false, 1);
    let run = run(
        &query,
        &subject,
        &[(0, 40, 0, 40), (8, 50, 60, 102)],
        &p,
        1e6,
        10.0,
        true,
    );
    let result = run.result.unwrap();
    let list = alignments(&result);
    assert_eq!(list.len(), 2);
    let (second, first) = (list[0], list[1]);
    assert!(first.score >= second.score);
    // The first alignment is the rectangle's best local alignment, on the adjusted matrix.
    let matrix = run.matrix.as_ref().unwrap();
    let rows: Vec<[i32; 28]> = query
        .iter()
        .map(|&q| matrix.scores[usize::from(q)])
        .collect();
    assert_eq!(
        first.score,
        whole_matrix_best(
            &rows,
            &subject,
            p.gapping_params.gap_open,
            p.gapping_params.gap_extend
        )
    );
    // The second copy's subject range does not meet the first's.
    assert!(second.match_end <= first.match_start || first.match_end <= second.match_start);
    for align in [first, second] {
        assert_eq!(
            align.matrix_adjust_rule,
            EMatrixAdjustRule::CompoScaleOldMatrix
        );
        assert_eq!(align.query_index, 0);
        assert!(matches!(
            align.context,
            Some(BlastCompoAlignmentContext::EditScript(ref script)) if !script.is_empty()
        ));
    }
    // The heap is asked once, before the first alignment.
    assert_eq!(run.heap_asks, 1);
    assert!(result.lambda_ratio.is_some());
}

// One HSP: one round, one alignment.
#[test]
fn one_hsp_gives_one_alignment() {
    let query = encode(QUERY);
    let subject_text = format!("{}{}{}", &QUERY[..40], SPACER, &QUERY[8..50]);
    let subject = encode(&subject_text);
    let run = run(
        &query,
        &subject,
        &[(0, 40, 0, 40)],
        &params(false, 1),
        1e6,
        10.0,
        true,
    );
    assert_eq!(alignments(&run.result.unwrap()).len(), 1);
}

// A heap that would not take the first alignment: no alignment from this match.
#[test]
fn a_heap_that_would_not_insert_leaves_no_alignment() {
    let query = encode(QUERY);
    let subject = encode(&format!("{}{}{}", &QUERY[..40], SPACER, &QUERY[8..50]));
    let run = run(
        &query,
        &subject,
        &[(0, 40, 0, 40), (8, 50, 60, 102)],
        &params(false, 1),
        1e6,
        10.0,
        false,
    );
    assert!(alignments(&run.result.unwrap()).is_empty());
    assert_eq!(run.heap_asks, 1);
}

// Sum statistics: the cutoff score decides, the heap is not asked.
#[test]
fn sum_statistics_read_the_cutoff_score() {
    let query = encode(QUERY);
    let subject = encode(&format!("{}{}{}", &QUERY[..40], SPACER, &QUERY[8..50]));
    let hsps = [(0, 40, 0, 40), (8, 50, 60, 102)];
    let none = run(
        &query,
        &subject,
        &hsps,
        &params(true, i32::MAX),
        1e6,
        10.0,
        true,
    );
    assert!(alignments(&none.result.unwrap()).is_empty());
    assert_eq!(none.heap_asks, 0);
    let both = run(
        &query,
        &subject,
        &hsps,
        &params(true, 2000),
        1e6,
        10.0,
        true,
    );
    assert_eq!(alignments(&both.result.unwrap()).len(), 2);
    assert_eq!(both.heap_asks, 0);
}

// A significant score of 0 (`-evalue` above searchsp x K) with two HSPs: NCBI's loop comes to a
// round of score 0 at the first cell; the first query letter W against the subject's first letter
// W scores above 0 (endless), against D 0 or less (the redo fails).
#[test]
fn a_significant_score_zero_is_decided_by_the_first_cell() {
    let query = encode(&format!("W{QUERY}"));
    let hsps = [(1, 41, 1, 41), (9, 51, 61, 103)];
    for (first, endless) in [("W", true), ("D", false)] {
        let subject = encode(&format!(
            "{first}{}{}{}",
            &QUERY[..40],
            SPACER,
            &QUERY[8..50]
        ));
        let error = run(&query, &subject, &hsps, &params(false, 1), 1e6, 1e6, true)
            .result
            .unwrap_err();
        assert_eq!(error.is::<SwRedoEndless>(), endless, "{first}");
        assert_eq!(error.is::<SwRedoFailed>(), !endless, "{first}");
    }
    // With -evalue below e(0) = searchsp x K the same match ends.
    let subject = encode(&format!("W{}{}{}", &QUERY[..40], SPACER, &QUERY[8..50]));
    let run = run(&query, &subject, &hsps, &params(false, 1), 1e6, 10.0, true);
    assert_eq!(alignments(&run.result.unwrap()).len(), 2);
}

/// A deterministic generator of the random matches (xorshift64).
struct Rng(u64);

impl Rng {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % n as u64) as usize
    }
}

/// A call's outcome, an error by its text.
fn outcome(run: &Run) -> std::result::Result<&BlastRedoOneMatchResult, String> {
    run.result.as_ref().map_err(|error| format!("{error:#}"))
}

// Random matches of the query's parts: one to three copies of a part in the subject, each with
// letters changed at random and an HSP of its own, after a first letter W or D or none; the query
// with a first W or without; composition modes 1, 2 and 3; `-evalue` 10; 1e6, at which a score
// of 0 is significant, so that a match of two or more HSPs fails or never ends (by the first
// cell); or 1e-100, at which no first round is significant, so the heap is never asked (the
// call keeps no state). A heap test left open ends the call as its FALSE does, and the call resumed from the state
// that it kept ends as the TRUE answer does (alignments, p-value, lambda ratio, the matrix left,
// or the failed or endless redo), without asking the heap again.
#[test]
fn a_redo_resumed_at_its_open_heap_test_equals_the_redo_answered_true() {
    const LETTERS: &[u8] = b"ARNDCQEGHILKMFPSTWYV";
    let modes = [
        BlastCompoAdjustMode::CompositionBasedStats,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        BlastCompoAdjustMode::ForceFullMatrixAdjust,
    ];
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let (mut aligned, mut failed, mut endless, mut not_asked) = (0, 0, 0, 0);
    for case in 0..400 {
        let w_first = rng.below(2) == 0;
        let query = encode(&format!("{}{QUERY}", if w_first { "W" } else { "" }));
        let offset = i32::from(w_first);
        let mut subject_text = match rng.below(3) {
            0 => "W".to_string(),
            1 => "D".to_string(),
            _ => String::new(),
        };
        let mut hsps = Vec::new();
        for _ in 0..1 + rng.below(3) {
            let length = 12 + rng.below(30);
            let start = rng.below(QUERY.len() - length);
            let part: String = QUERY[start..start + length]
                .bytes()
                .map(|letter| {
                    if rng.below(6) == 0 {
                        char::from(LETTERS[rng.below(LETTERS.len())])
                    } else {
                        char::from(letter)
                    }
                })
                .collect();
            let at = subject_text.len() as i32;
            subject_text.push_str(&part);
            let (start, length) = (start as i32 + offset, length as i32);
            hsps.push((start, start + length, at, at + length));
            subject_text.push_str(&SPACER[..1 + rng.below(SPACER.len())]);
        }
        let subject = encode(&subject_text);
        let p = params_in(modes[rng.below(modes.len())], false, 1);
        let evalues = (
            1e6,
            match rng.below(4) {
                0 => 1e6,
                1 => 1e-100,
                _ => 10.0,
            },
        );
        let call = |resume: Option<&SwMatchResume>, heap: &mut dyn FnMut() -> SwHeapAnswer| {
            run_from(&query, &subject, &hsps, &p, evalues, resume, heap)
        };
        let (with_true, _) = call(None, &mut || SwHeapAnswer::Known(true));
        let (with_false, _) = call(None, &mut || SwHeapAnswer::Known(false));
        let (open, kept) = call(None, &mut || SwHeapAnswer::Deferred);
        assert_eq!(outcome(&open), outcome(&with_false), "case {case}");
        assert_eq!(open.matrix, with_false.matrix, "case {case}");
        assert_eq!(kept.is_some(), with_true.heap_asks == 1, "case {case}");
        let Some(kept) = kept else {
            not_asked += 1;
            continue;
        };
        let (resumed, again) = call(Some(&kept), &mut || {
            panic!("case {case}: the resumed redo asks the heap again")
        });
        assert!(again.is_none(), "case {case}");
        assert_eq!(outcome(&resumed), outcome(&with_true), "case {case}");
        assert_eq!(resumed.matrix, with_true.matrix, "case {case}");
        match &with_true.result {
            Ok(result) if !alignments(result).is_empty() => aligned += 1,
            Err(error) if error.is::<SwRedoFailed>() => failed += 1,
            Err(error) if error.is::<SwRedoEndless>() => endless += 1,
            _ => {}
        }
    }
    assert!(
        aligned > 100 && failed > 10 && endless > 10 && not_asked > 0,
        "aligned {aligned} failed {failed} endless {endless} not asked {not_asked}"
    );
}
