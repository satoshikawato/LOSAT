//! Hand-made matches for the Smith-Waterman match body on TBLASTN's translated subject windows
//! (`redo_one_match_smith_waterman` with `subject_is_translated`): the heap test before the
//! query's first alignment of the call (not of each window), the redo without composition-based
//! statistics (the unadjusted matrix, no adjustment rule, no SEG) and an error of the heap test.
//! BLASTP's protein windows are tested through BLASTP's call
//! (`algorithm::blastp::kappa::smith_waterman` tests) and its manifest rows.

use super::*;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::{build_matrix_info, read_aa_composition};
use crate::core::composition_adjustment::redo_alignment::{
    build_query_word_hashes, translated_subject_get_range, BlastCompoAdjustMode,
    BlastCompoGappingParams, BlastCompoSequenceData, BlastRedoAlignCallbacks,
};
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::aa_char_to_ncbistdaa;

const QUERY: &str = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV";
const SPACER: &str = "PGSNPGSNPGSNPGSNPGSN";

fn encode(text: &str) -> Vec<u8> {
    text.bytes().map(aa_char_to_ncbistdaa).collect()
}

/// One codon of the standard code per letter.
fn back_translate(protein: &str) -> String {
    protein
        .chars()
        .map(|aa| match aa {
            'A' => "GCT",
            'R' => "CGT",
            'N' => "AAT",
            'D' => "GAT",
            'C' => "TGT",
            'Q' => "CAA",
            'E' => "GAA",
            'G' => "GGT",
            'H' => "CAT",
            'I' => "ATT",
            'L' => "CTT",
            'K' => "AAA",
            'M' => "ATG",
            'F' => "TTT",
            'P' => "CCT",
            'S' => "TCT",
            'T' => "ACT",
            'W' => "TGG",
            'Y' => "TAT",
            'V' => "GTT",
            other => panic!("no codon for {other}"),
        })
        .collect()
}

fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// TBLASTN's redo parameters with BLOSUM62 11/1: composition mode 2 scaled by 32, or none
/// (scale 1, the unadjusted matrix).
fn params(
    mode: BlastCompoAdjustMode,
    do_link_hsps: bool,
    cutoff_score: i32,
) -> BlastRedoAlignParams {
    let scale = if mode == BlastCompoAdjustMode::NoCompositionBasedStats {
        1.0
    } else {
        32.0
    };
    BlastRedoAlignParams {
        matrix_info: build_matrix_info(ScoringMatrix::Blosum62, 0.3176 / scale).unwrap(),
        gapping_params: BlastCompoGappingParams {
            gap_open: (11.0 * scale) as i32,
            gap_extend: scale as i32,
            x_dropoff: (2000.0 * scale) as i32,
            context: std::cell::Cell::new(None),
        },
        compo_adjust_mode: mode,
        alphsize: crate::utils::matrix::BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: 0.0,
        score_divisor: scale,
        smith_waterman: true,
        near_identical_cutoff: 0.0,
        position_based: false,
        ccat_query_length: 0,
        query_is_translated: false,
        subject_is_translated: true,
        cutoff_score,
        do_link_hsps,
    }
}

#[allow(clippy::too_many_arguments)]
fn no_x_drop_redo(
    _: &BlastCompoAlignment,
    _: EMatrixAdjustRule,
    _: Option<&AdjustedProteinMatrix>,
    _: &BlastCompoSequenceData,
    _: &crate::core::composition_adjustment::redo_alignment::BlastCompoSequenceRange,
    _: i32,
    _: &BlastCompoSequenceData,
    _: &crate::core::composition_adjustment::redo_alignment::BlastCompoSequenceRange,
    _: i32,
    _: &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    panic!("the Smith-Waterman redo makes no X-drop redo of an HSP")
}

const CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range: translated_subject_get_range,
    redo_one_alignment: no_x_drop_redo,
};

/// The incoming alignments of query 0: `(query_start, query_end, frame, subject_start,
/// subject_end)` in the frame's letters.
fn incoming(hsps: &[(i32, i32, i32, i32, i32)]) -> Option<Box<BlastCompoAlignment>> {
    let mut head = None;
    for &(qs, qe, frame, ss, se) in hsps.iter().rev() {
        let mut align = blast_compo_alignment_new(
            (qe - qs) * 4,
            EMatrixAdjustRule::DontAdjustMatrix,
            qs,
            qe,
            0,
            ss,
            se,
            frame,
            None,
        );
        align.next = head;
        head = Some(align);
    }
    head
}

struct Run {
    result: Result<BlastRedoOneMatchResult>,
    heap_asks: usize,
}

/// The subject: the query's letters 0-40 in frame +1, a spacer, and the query's letters 8-50 on the
/// minus strand (frame -1, letters 0-42 of that frame); one HSP in each frame, two windows of
/// query 0.
fn run(
    params: &BlastRedoAlignParams,
    searchsp: f64,
    cutoff_e: f64,
    heap: &mut dyn FnMut() -> Result<bool>,
) -> Run {
    let query = encode(QUERY);
    let subject = format!(
        "{}{}{}",
        back_translate(&QUERY[..40]),
        back_translate(SPACER),
        reverse_complement(&back_translate(&QUERY[8..50]))
    );
    let info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(&query),
        composition: read_aa_composition(&query),
        eff_search_space: searchsp,
        words: Some(build_query_word_hashes(&query)),
    };
    let matching_seq =
        BlastCompoMatchingSequence::new_translated(3, subject.as_bytes().into(), 1).unwrap();
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    let mut matrix = None;
    let unadjusted = unadjusted_square_matrix(ScoringMatrix::Blosum62);
    let mut heap_asks = 0;
    let mut scratch = SwMatchScratch::new();
    let scale = params.score_divisor;
    let result = redo_one_match_smith_waterman(
        &incoming(&[(0, 40, 1, 0, 40), (8, 50, -1, 0, 42)]),
        params,
        0.267 / scale,
        0.041f64.glibc_ln(),
        &matching_seq,
        std::slice::from_ref(&info),
        &CALLBACKS,
        &mut workspace,
        &mut matrix,
        Some(&unadjusted),
        cutoff_e,
        &mut |query_index, _evalue, _score| {
            assert_eq!(query_index, 0);
            heap_asks += 1;
            heap()
        },
        &mut scratch,
    );
    Run { result, heap_asks }
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

// Two windows of one query (frames +1 and -1): NCBI asks the heap before the query's first
// alignment of the call (while the query's list of alignments is empty), so once, not once per
// window.
#[test]
fn the_heap_is_asked_before_the_querys_first_alignment_of_the_call() {
    for mode in [
        BlastCompoAdjustMode::NoCompositionBasedStats,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
    ] {
        let run = run(&params(mode, false, 1), 1e6, 10.0, &mut || Ok(true));
        let result = run.result.unwrap();
        let list = alignments(&result);
        assert_eq!(list.len(), 2, "{mode:?}");
        let mut frames: Vec<_> = list.iter().map(|align| align.frame).collect();
        frames.sort_unstable();
        assert_eq!(frames, [-1, 1], "{mode:?}");
        assert_eq!(run.heap_asks, 1, "{mode:?}");
    }
}

// A heap that would not take the query's first alignment: the first window makes none, so the
// second window asks again (its query still has none), and the match has no alignment.
#[test]
fn a_heap_that_refuses_is_asked_again_by_the_next_window() {
    let run = run(
        &params(BlastCompoAdjustMode::NoCompositionBasedStats, false, 1),
        1e6,
        10.0,
        &mut || Ok(false),
    );
    assert!(alignments(&run.result.unwrap()).is_empty());
    assert_eq!(run.heap_asks, 2);
}

// Without composition-based statistics the rows are the unadjusted matrix's: an ungapped copy of
// the query scores its BLOSUM62 sum, the rule is `eDontAdjustMatrix`, no p-value or lambda ratio.
#[test]
fn without_composition_statistics_the_matrix_is_not_adjusted() {
    let run = run(
        &params(BlastCompoAdjustMode::NoCompositionBasedStats, false, 1),
        1e6,
        10.0,
        &mut || Ok(true),
    );
    let result = run.result.unwrap();
    let query = encode(QUERY);
    let unadjusted = unadjusted_square_matrix(ScoringMatrix::Blosum62);
    let copy_score = |range: std::ops::Range<usize>| -> i32 {
        query[range]
            .iter()
            .map(|&letter| unadjusted.scores[usize::from(letter)][usize::from(letter)])
            .sum()
    };
    for align in alignments(&result) {
        assert_eq!(
            align.matrix_adjust_rule,
            EMatrixAdjustRule::DontAdjustMatrix
        );
        let expected = if align.frame == 1 {
            copy_score(0..40)
        } else {
            copy_score(8..50)
        };
        assert_eq!(align.score, expected, "frame {}", align.frame);
    }
    assert_eq!(result.pvalue_for_this_pair, None);
    assert_eq!(result.lambda_ratio, None);
}

// Sum statistics: the cutoff score decides, the heap is not asked.
#[test]
fn sum_statistics_read_the_cutoff_score_in_every_window() {
    let run = run(
        &params(BlastCompoAdjustMode::NoCompositionBasedStats, true, 1),
        1e6,
        10.0,
        &mut || panic!("no heap test with sum statistics"),
    );
    assert_eq!(alignments(&run.result.unwrap()).len(), 2);
    assert_eq!(run.heap_asks, 0);
}

// An error of the heap test ends the call with it.
#[test]
fn an_error_of_the_heap_test_ends_the_call() {
    let run = run(
        &params(BlastCompoAdjustMode::NoCompositionBasedStats, false, 1),
        1e6,
        10.0,
        &mut || Err(anyhow::anyhow!("heap past the array")),
    );
    let error = run.result.unwrap_err();
    assert_eq!(error.to_string(), "heap past the array");
    assert_eq!(run.heap_asks, 1);
}

// The refusal and NCBI's error of a failed batch name the program.
#[test]
fn the_batch_errors_name_the_program() {
    let crash = sw_batch_failed_error("tblastn", true).to_string();
    assert!(
        crash.contains("NCBI BLAST+ 2.17.0's tblastn crashes"),
        "{crash}"
    );
    assert!(
        crash.ends_with("is not supported by LOSAT's TBLASTN"),
        "{crash}"
    );
    let failed = sw_batch_failed_error("tblastn", false);
    assert!(failed
        .to_string()
        .contains("LOSAT engine error: Traceback failed"));
    let endless = sw_endless_error("blastp", "why").to_string();
    assert_eq!(
        endless,
        "a match whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 (why) is not supported by LOSAT's BLASTP"
    );
}

/// A batch status whose last redone match failed.
fn failed_status() -> SwBatchStatus {
    let mut status = SwBatchStatus::default();
    let redone = status.redone::<()>(Err(SwRedoFailed.into()), || anyhow::anyhow!("not endless"));
    assert!(matches!(redone, Ok(None)));
    assert!(status.failed());
    status
}

// Line RG item X4b-1 (review 2, S3 a): NCBI's order at each match of the redo's stream, early
// termination first, then the fetch (`SwBatchStatus::reaches_redo`). After a failed redo, a match
// whose subject the fetch excludes clears the status unless early termination stops it first; a
// match that the test stops keeps the status whether or not its subject is excluded, so the failed
// redo before it still fails the batch. Testing the exclusion first, or not testing early
// termination, fails this test.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3575
#[test]
fn an_excluded_list_clears_the_status_unless_early_termination_stops_it() {
    // (terminates early, excluded at the fetch, redone, the batch still failed)
    for (terminates_early, excluded, redone, still_failed) in [
        (true, true, false, true),
        (true, false, false, true),
        (false, true, false, false),
        (false, false, true, true),
    ] {
        let mut status = failed_status();
        let case = format!("terminates early {terminates_early}, excluded {excluded}");
        assert_eq!(
            status.reaches_redo(terminates_early, excluded),
            redone,
            "{case}"
        );
        assert_eq!(status.failed(), still_failed, "{case}");
    }
    // A match that is redone sets the status by its own redo.
    let mut status = failed_status();
    assert!(status.reaches_redo(false, false));
    let redone = status.redone(Ok(7), || anyhow::anyhow!("not endless"));
    assert!(matches!(redone, Ok(Some(7))));
    assert!(!status.failed());
    // Without a negative seqid list no list of the stream is marked.
    let lists = [3u32, 1, 2];
    let (marks, _) = crate::algorithm::common::database::FetchMarks::of_noting(
        None,
        &lists,
        |&oid| Some(oid),
        |_| (),
    );
    assert!((0..=lists.len()).all(|index| !marks.excluded(index)));
}
