use super::*;
use crate::core::composition_adjustment::adjust_scores::{
    blast_adjust_scores, build_matrix_info, read_aa_composition, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::redo_calc_lambda;
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;
use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa};
use crate::utils::seg::{SegMasker, SegParams};
use std::fs;

/// A file of `tests/fasta`, from the crate directory and not the working directory: a WASI
/// module has no working directory (the runner preopens "/").
fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

fn fasta_sequence_by_id(contents: &str, wanted_id: &str) -> Vec<u8> {
    let mut current_id = None::<&str>;
    let mut sequence = Vec::new();
    for line in contents.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let id = header.split_whitespace().next().unwrap_or_default();
            if current_id == Some(wanted_id) {
                break;
            }
            current_id = Some(id);
            continue;
        }
        if current_id == Some(wanted_id) {
            sequence.extend(line.trim().bytes().map(aa_char_to_ncbistdaa));
        }
    }
    sequence
}

fn mask_query_sequence(sequence: &[u8]) -> Vec<u8> {
    let mut masked = sequence.to_vec();
    let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));
    for interval in masker.mask_sequence(sequence) {
        for pos in interval.start..interval.end {
            masked[pos] = ncbistdaa::X;
        }
    }
    masked
}

fn mask_subject_sequence(sequence: &[u8]) -> Vec<u8> {
    let mut masked = sequence.to_vec();
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));
    for interval in masker.mask_sequence(sequence) {
        for pos in interval.start..interval.end {
            masked[pos] = ncbistdaa::X;
        }
    }
    masked
}

#[test]
fn test_adjust_subject_range_keeps_short_subject() {
    let (offset, length, shift) = adjust_subject_range(100, 1000, 50, 400);
    assert_eq!(offset, 100);
    assert_eq!(length, 1000);
    assert_eq!(shift, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4176-4190
#[test]
fn test_adjust_subject_range_shifts_long_subject_ncbi_style() {
    let subject_offset = 7_000usize;
    let subject_length = 100_000usize;
    let query_offset = 5usize;
    let query_length = 8usize;
    let (adjusted_offset, adjusted_length, start_shift) =
        adjust_subject_range(subject_offset, subject_length, query_offset, query_length);

    assert_eq!(adjusted_offset, query_offset + MAX_TOTAL_GAPS);
    assert_eq!(
        start_shift,
        subject_offset - (query_offset + MAX_TOTAL_GAPS)
    );
    assert_eq!(
        adjusted_length,
        subject_length.min(subject_offset + 3_003) - start_shift
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2488-2533
#[test]
fn test_prelim_edit_ops_to_gap_edit_script_reverses_forward_ops_ncbi_style() {
    let merged = prelim_edit_ops_to_gap_edit_script(
        vec![GapEditOp::Sub(2), GapEditOp::Del(1)],
        vec![GapEditOp::Ins(3), GapEditOp::Del(4), GapEditOp::Sub(5)],
    );
    assert_eq!(
        merged,
        vec![
            GapEditOp::Sub(2),
            GapEditOp::Del(1),
            GapEditOp::Sub(5),
            GapEditOp::Del(4),
            GapEditOp::Ins(3),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2496-2533
#[test]
fn test_prelim_edit_ops_to_gap_edit_script_merges_last_forward_op_ncbi_style() {
    let merged = prelim_edit_ops_to_gap_edit_script(
        vec![GapEditOp::Sub(2), GapEditOp::Del(1)],
        vec![GapEditOp::Ins(3), GapEditOp::Del(4), GapEditOp::Del(5)],
    );
    assert_eq!(
        merged,
        vec![
            GapEditOp::Sub(2),
            GapEditOp::Del(6),
            GapEditOp::Del(4),
            GapEditOp::Ins(3),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4681-4707
#[test]
fn test_prune_terminal_gap_ops_updates_score_and_boundaries_ncbi_style() {
    let mut edit_script = vec![GapEditOp::Ins(2), GapEditOp::Sub(3), GapEditOp::Del(1)];
    let mut score_left = 30;
    let mut score_right = 40;
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    let mut query_start = 10i32;
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    let mut query_stop = 15i32;
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    let mut subject_start = 20i32;
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    let mut subject_stop = 24i32;

    prune_terminal_gap_ops(
        &mut edit_script,
        11,
        1,
        &mut score_left,
        &mut score_right,
        &mut query_start,
        &mut query_stop,
        &mut subject_start,
        &mut subject_stop,
    );

    assert_eq!(edit_script, vec![GapEditOp::Sub(3)]);
    assert_eq!(score_left, 43);
    assert_eq!(score_right, 52);
    assert_eq!(query_start, 12);
    assert_eq!(query_stop, 15);
    assert_eq!(subject_start, 20);
    assert_eq!(subject_stop, 23);
}

#[test]
fn test_blastp_traceback_builds_gap_script() {
    let query = [1u8, 2, 3, 4, 5];
    let subject = [1u8, 2, 9, 4, 5];
    let result = blast_gapped_alignment_with_traceback(
        &query,
        &subject,
        1,
        1,
        ScoringMatrix::Blosum62,
        None,
        11,
        1,
        25,
    )
    .expect("traceback alignment");
    assert!(!result.edit_script.is_empty());
    assert!(result.score > 0);
    assert!(result.query_stop > result.query_start);
}

#[test]
fn test_reevaluate_blastp_hit_keeps_best_subalignment() {
    let query = [1u8, 2, 3, 4];
    let subject = [1u8, 2, 21, 4];
    let mut hit = Hit {
        identity: 0.0,
        length: 4,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 4,
        s_start: 1,
        s_end: 4,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 4,
        q_idx: 0,
        s_idx: 0,
        raw_score: 0,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![GapEditOp::Sub(4)]),
        num_positives: 0,
    };
    let delete = reevaluate_blastp_hit_with_traceback(
        &mut hit,
        &query,
        &query,
        &subject,
        ScoringMatrix::Blosum62,
        11,
        1,
        1,
    );
    assert!(!delete);
    assert!(hit.raw_score > 0);
    assert!(hit.q_end >= hit.q_start);
}

#[test]
fn test_reevaluate_blastp_hit_trims_partial_leading_sub_run() {
    let matrix = ScoringMatrix::Blosum62;
    let mismatch = (1u8..=u8::MAX)
        .find(|&residue| protein_score(matrix, 1, residue) < 0)
        .expect("negative mismatch for residue 1");
    let match_score = protein_score(matrix, 1, 1);
    let query = [1u8, 1, 1, 1];
    let subject = [mismatch, 1, 1, 1];
    let mut hit = Hit {
        identity: 0.0,
        length: 4,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 4,
        s_start: 1,
        s_end: 4,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 4,
        q_idx: 0,
        s_idx: 0,
        raw_score: 0,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![GapEditOp::Sub(4)]),
        num_positives: 0,
    };
    let delete =
        reevaluate_blastp_hit_with_traceback(&mut hit, &query, &query, &subject, matrix, 11, 1, 1);

    assert!(!delete);
    assert_eq!(hit.q_start, 2);
    assert_eq!(hit.q_end, 4);
    assert_eq!(hit.s_start, 2);
    assert_eq!(hit.s_end, 4);
    assert_eq!(hit.raw_score, match_score * 3);
    assert_eq!(hit.gap_info, Some(vec![GapEditOp::Sub(3)]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:618-637
#[test]
fn test_reevaluate_blastp_hit_extends_left_over_exact_low_residue_match() {
    let matrix = ScoringMatrix::Blosum62;
    let mismatch = (1u8..=u8::MAX)
        .find(|&residue| protein_score(matrix, 5, residue) < 0)
        .expect("negative mismatch for residue 5");
    let match_score = protein_score(matrix, 1, 1);
    let query = [9u8, 9, 1, 1, 5];
    let subject = [1u8, 1, mismatch];
    let mut hit = Hit {
        identity: 0.0,
        length: 3,
        mismatch: 0,
        gapopen: 1,
        q_start: 3,
        q_end: 5,
        s_start: 2,
        s_end: 3,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 0,
        query_length: query.len(),
        q_idx: 0,
        s_idx: 0,
        raw_score: 0,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![GapEditOp::Ins(1), GapEditOp::Sub(2)]),
        num_positives: 0,
    };

    let delete =
        reevaluate_blastp_hit_with_traceback(&mut hit, &query, &query, &subject, matrix, 11, 1, 1);

    assert!(!delete);
    assert_eq!(hit.q_start, 3);
    assert_eq!(hit.q_end, 4);
    assert_eq!(hit.s_start, 1);
    assert_eq!(hit.s_end, 2);
    assert_eq!(hit.raw_score, match_score);
    assert_eq!(hit.gap_info, Some(vec![GapEditOp::Sub(2)]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:547-549
#[test]
fn test_reevaluate_blastp_hit_masks_query_residue_with_ncbi_k_residue_mask() {
    use crate::utils::matrix::ncbistdaa;

    let matrix = ScoringMatrix::Blosum62;
    let query = [ncbistdaa::S];
    let subject = [ncbistdaa::S];
    let mut hit = Hit {
        identity: 0.0,
        length: 1,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 1,
        s_start: 1,
        s_end: 1,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 1,
        q_idx: 0,
        s_idx: 0,
        raw_score: 0,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(vec![GapEditOp::Sub(1)]),
        num_positives: 0,
    };

    let delete = reevaluate_blastp_hit_with_traceback(
        &mut hit, &query, &query, &subject, matrix, 11, 1, -10,
    );

    assert!(!delete);
    assert_eq!(
        hit.raw_score,
        protein_score(matrix, ncbistdaa::S & 0x0f, ncbistdaa::S)
    );
    assert_eq!(hit.num_ident, 1);
}

#[test]
fn test_traceback_keeps_seed_when_start_is_last_residue() {
    let query = [1u8, 2, 3];
    let subject = [8u8, 9, 3];
    let result = blast_gapped_alignment_with_traceback(
        &query,
        &subject,
        2,
        2,
        ScoringMatrix::Blosum62,
        None,
        11,
        1,
        25,
    )
    .expect("traceback keeps terminal seed");
    assert_eq!(result.query_start, 2);
    assert_eq!(result.query_stop, 3);
    assert_eq!(result.subject_start, 2);
    assert_eq!(result.subject_stop, 3);
}

#[test]
fn test_score_only_alignment_keeps_seed_when_start_is_last_residue() {
    let query = [1u8, 2, 3];
    let subject = [8u8, 9, 3];
    let result = blastp_score_only_gapped_alignment(
        &query,
        &subject,
        2,
        2,
        ScoringMatrix::Blosum62,
        11,
        1,
        25,
        BlastpGappedAlignmentMode::Exact,
    );
    assert_eq!(result.query_start, 2);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4338-4344
    assert_eq!(result.query_stop, 2);
    assert_eq!(result.subject_start, 2);
    assert_eq!(result.subject_stop, 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4257-4291
#[test]
fn test_score_only_alignment_exact_preserves_original_subject_coordinates_with_shift() {
    let query = vec![1u8; 8];
    let subject = vec![1u8; 100_000];
    let q_start = 5usize;
    let s_start = 7_000usize;

    let result = blastp_score_only_gapped_alignment(
        &query,
        &subject,
        q_start,
        s_start,
        ScoringMatrix::Blosum62,
        11,
        1,
        25,
        BlastpGappedAlignmentMode::Exact,
    );

    assert_eq!(result.query_start, 0);
    assert_eq!(result.query_stop, 8);
    assert_eq!(result.subject_start, 6_995);
    assert_eq!(result.subject_stop, 7_003);
    assert!(result.score > 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4273-4291
#[test]
fn test_score_only_alignment_restricted_preserves_original_subject_coordinates_with_shift() {
    let query = vec![1u8; 8];
    let subject = vec![1u8; 100_000];
    let q_start = 5usize;
    let s_start = 7_000usize;

    let result = blastp_score_only_gapped_alignment(
        &query,
        &subject,
        q_start,
        s_start,
        ScoringMatrix::Blosum62,
        11,
        1,
        25,
        BlastpGappedAlignmentMode::Restricted,
    );

    assert_eq!(result.query_start, 0);
    assert_eq!(result.query_stop, 8);
    assert_eq!(result.subject_start, 6_995);
    assert_eq!(result.subject_stop, 7_003);
    assert!(result.score > 0);
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:797-808
#[test]
fn test_score_only_alignment_reused_scratch_matches_fresh_scratch() {
    let query = [
        ncbistdaa::A,
        ncbistdaa::C,
        ncbistdaa::D,
        ncbistdaa::E,
        ncbistdaa::F,
        ncbistdaa::G,
        ncbistdaa::H,
        ncbistdaa::I,
    ];
    let subject = [
        ncbistdaa::A,
        ncbistdaa::C,
        ncbistdaa::X,
        ncbistdaa::E,
        ncbistdaa::F,
        ncbistdaa::G,
        ncbistdaa::K,
        ncbistdaa::I,
    ];
    let cases = [
        (2usize, 2usize, BlastpGappedAlignmentMode::Exact),
        (5usize, 5usize, BlastpGappedAlignmentMode::Exact),
        (5usize, 5usize, BlastpGappedAlignmentMode::Restricted),
    ];
    let mut scratch = GapAlignScratch::new();

    for (q_start, s_start, mode) in cases {
        let fresh = blastp_score_only_gapped_alignment(
            &query,
            &subject,
            q_start,
            s_start,
            ScoringMatrix::Blosum62,
            11,
            1,
            25,
            mode,
        );
        let reused = blastp_score_only_gapped_alignment_with_scratch(
            &query,
            &subject,
            q_start,
            s_start,
            ScoringMatrix::Blosum62,
            11,
            1,
            25,
            mode,
            &mut scratch,
        );
        assert_eq!(reused, fresh);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4715-4716
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:523-556
#[test]
fn test_traceback_score_matches_rescored_edit_script_wssv_537() {
    let query_fasta = fs::read_to_string(fixture("WSSV.faa")).expect("query fasta");
    let subject_fasta = fs::read_to_string(fixture("PajaWSV.faa")).expect("subject fasta");
    let query = mask_query_sequence(&fasta_sequence_by_id(&query_fasta, "YP_009220537.1"));
    let subject = mask_subject_sequence(&fasta_sequence_by_id(&subject_fasta, "BDT63528.1"));

    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176 / 32.0).unwrap();
    let query_composition = read_aa_composition(&query);
    let subject_composition = read_aa_composition(&subject);
    let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
    let adjusted = blast_adjust_scores(
        &matrix_info,
        &query_composition,
        query.len() as i32,
        &subject_composition,
        subject.len() as i32,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        0,
        &mut composition_workspace,
        redo_calc_lambda,
    )
    .unwrap()
    .expect("adjusted matrix");

    let aligned = blast_gapped_alignment_with_traceback(
        &query,
        &subject,
        8,
        37,
        ScoringMatrix::Blosum62,
        Some(&adjusted.adjusted_matrix),
        11 * 32,
        32,
        100_000,
    )
    .expect("traceback alignment");
    let rescored = score_edit_script_local_alignment(
        &query,
        &subject,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:523-550: query/subject pointers walk the valid edit script.
        usize::try_from(aligned.query_start).unwrap(),
        usize::try_from(aligned.subject_start).unwrap(),
        &aligned.edit_script,
        BlastpScoreMatrix::Adjusted(&adjusted.adjusted_matrix),
        11 * 32,
        32,
    );

    assert_eq!(aligned.score, rescored);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1902-1945
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4549-4715
#[test]
#[ignore = "diagnostic start scan for the residual WSSV/PajaWSV CBS=2 traceback diff"]
fn test_debug_traceback_start_scan_wssv_488_vs_bdt63510() {
    fn push_aligned_op(edit_script: &mut Vec<GapEditOp>, op: GapEditOp) {
        match (edit_script.last_mut(), op) {
            (Some(GapEditOp::Sub(lhs)), GapEditOp::Sub(rhs))
            | (Some(GapEditOp::Del(lhs)), GapEditOp::Del(rhs))
            | (Some(GapEditOp::Ins(lhs)), GapEditOp::Ins(rhs)) => *lhs += rhs,
            (_, op) => edit_script.push(op),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:523-556
    fn edit_script_from_aligned_strings(query_aln: &str, subject_aln: &str) -> Vec<GapEditOp> {
        let mut edit_script = Vec::new();
        for (q, s) in query_aln.bytes().zip(subject_aln.bytes()) {
            match (q == b'-', s == b'-') {
                (false, false) => push_aligned_op(&mut edit_script, GapEditOp::Sub(1)),
                (true, false) => push_aligned_op(&mut edit_script, GapEditOp::Del(1)),
                (false, true) => push_aligned_op(&mut edit_script, GapEditOp::Ins(1)),
                (true, true) => unreachable!("pairwise alignment cannot gap both rows"),
            }
        }
        edit_script
    }

    let query_fasta = fs::read_to_string(fixture("WSSV.faa")).expect("query fasta");
    let subject_fasta = fs::read_to_string(fixture("PajaWSV.faa")).expect("subject fasta");
    let query = fasta_sequence_by_id(&query_fasta, "YP_009220488.1");
    let subject = mask_subject_sequence(&fasta_sequence_by_id(&subject_fasta, "BDT63510.1"));

    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176 / 32.0).unwrap();
    let query_composition = read_aa_composition(&query);
    let subject_composition = read_aa_composition(&subject);
    let mut composition_workspace = BlastCompositionWorkspace::new_blosum62();
    let adjusted = blast_adjust_scores(
        &matrix_info,
        &query_composition,
        query.len() as i32,
        &subject_composition,
        subject.len() as i32,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        0,
        &mut composition_workspace,
        redo_calc_lambda,
    )
    .unwrap()
    .expect("adjusted matrix");

    let pairwise = fs::read_to_string("tests/blast_out/blastp_pairwise/WSSV.PajaWSV.blastp.out")
        .expect("pairwise reference");
    let mut in_target_query = false;
    let mut in_target_subject = false;
    let mut query_aln = String::new();
    let mut subject_aln = String::new();
    for line in pairwise.lines() {
        if line.starts_with("Query= ") {
            in_target_query = line.contains("YP_009220488.1");
            in_target_subject = false;
            continue;
        }
        if !in_target_query {
            continue;
        }
        if line.starts_with("> ") {
            in_target_subject = line.contains("BDT63510.1");
            continue;
        }
        if !in_target_subject {
            continue;
        }
        if line.starts_with("Query") {
            query_aln.push_str(
                line.split_whitespace()
                    .nth(2)
                    .expect("pairwise query aligned chunk"),
            );
        } else if line.starts_with("Sbjct") {
            subject_aln.push_str(
                line.split_whitespace()
                    .nth(2)
                    .expect("pairwise subject aligned chunk"),
            );
        } else if !query_aln.is_empty() && line.starts_with("Lambda") {
            break;
        }
    }

    let ncbi_edit_script = edit_script_from_aligned_strings(&query_aln, &subject_aln);
    let ncbi_alignment_score = score_edit_script_local_alignment(
        &query,
        &subject,
        78,
        0,
        &ncbi_edit_script,
        BlastpScoreMatrix::Adjusted(&adjusted.adjusted_matrix),
        11 * 32,
        32,
    );
    println!("ncbi_alignment_score={ncbi_alignment_score}");

    let mut best = None;
    for s_start in 0..665usize {
        let q_start = s_start + 78;
        if q_start >= query.len() || s_start >= subject.len() {
            continue;
        }
        let Some(aligned) = blast_gapped_alignment_with_traceback(
            &query,
            &subject,
            q_start,
            s_start,
            ScoringMatrix::Blosum62,
            Some(&adjusted.adjusted_matrix),
            11 * 32,
            32,
            100_000,
        ) else {
            continue;
        };
        if best
            .as_ref()
            .is_none_or(|(score, _, _, _, _, _, _)| aligned.score > *score)
        {
            best = Some((
                aligned.score,
                q_start,
                s_start,
                aligned.query_start,
                aligned.query_stop,
                aligned.subject_start,
                aligned.subject_stop,
            ));
        }
        if aligned.score == 1532 {
            println!(
                "score=1532 q_start={} s_start={} q={}..{} s={}..{}",
                q_start,
                s_start,
                aligned.query_start,
                aligned.query_stop,
                aligned.subject_start,
                aligned.subject_stop
            );
        }
    }

    println!("best={best:?}");
    panic!("debug scan complete");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:569-575
#[test]
fn test_align_ex_protein_sets_fence_hit_on_ncbi_fence_sentry() {
    let query = [1u8, 1, 1, 1];
    let subject = [1u8, FENCE_SENTRY, 1, 1];
    let mut scratch = GapAlignScratch::new();
    let mut fence_hit = false;

    let (_a_offset, _b_offset, score, edit_ops) = align_ex_protein(
        &query,
        &subject,
        3,
        3,
        BlastpScoreMatrix::Standard(ScoringMatrix::Blosum62),
        11,
        1,
        18,
        false,
        &mut scratch,
        &mut fence_hit,
        false,
    );

    assert!(fence_hit);
    assert_eq!(score, 0);
    assert!(edit_ops.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4620-4652,4714-4716
// Independent pinned C boundary oracle returns status=0, fence=1 and
// endpoints 3..1 with no edit script. The caller owns the fence retry.
#[test]
fn test_blast_gapped_alignment_with_traceback_skips_right_extension_after_left_fence_hit() {
    let query = [1u8, 1, 1, 1];
    let subject = [1u8, 1, FENCE_SENTRY, 1];
    let mut scratch = GapAlignScratch::new();
    let mut fence_hit = false;
    let result = blast_gapped_alignment_with_traceback_with_scratch(
        &query,
        &subject,
        2,
        2,
        ScoringMatrix::Blosum62,
        None,
        11,
        1,
        18,
        &mut scratch,
        Some(&mut fence_hit),
    )
    .expect("NCBI returns successful status even after hitting the fence");
    let expected = include_str!("../../../tests/unit/blastp_fence_traceback_expected.tsv");
    let mut lines = expected.lines();
    let fields: Vec<i32> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(|f| f.parse().unwrap())
        .collect();
    assert_eq!(fields[0], 0); // successful C status maps to Some
    assert_eq!(i32::from(fence_hit), fields[1]);
    assert_eq!(result.score, fields[2]);
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    assert_eq!(result.query_start, fields[3]);
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    assert_eq!(result.query_stop, fields[4]);
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    assert_eq!(result.subject_start, fields[5]);
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:85-88
    assert_eq!(result.subject_stop, fields[6]);
    assert_eq!(lines.next(), Some("SCRIPT\t-1"));
    assert!(result.edit_script.is_empty());
}
// NCBI 2.17.0: c++/src/algo/blast/core/blast_gapalign.c:431-432,531-669,689-726
// A changing band and scratch reuse must preserve every used row byte,
// not only the resulting score or final tabular output.
#[test]
fn test_traceback_contiguous_rows_reused_scratch_matches_fresh() {
    // No NCBI counterpart: the test reads scalar-kernel scratch rows, so it calls the scalar
    // port (`align_ex_protein_scalar`) directly; the SIMD kernel writes no such rows.
    let mut adjusted = AdjustedProteinMatrix {
        scores: [[-3; 28]; 28],
    };
    for i in 1..28 {
        adjusted.scores[i][i] = 2;
    }
    let mut scratch = GapAlignScratch::new();
    let mut saw_offset = false;
    let mut saw_gap = false;
    let mut saw_fence = false;
    for score_matrix in [
        BlastpScoreMatrix::standard(ScoringMatrix::Blosum62),
        BlastpScoreMatrix::Standard(ScoringMatrix::Blosum50),
        BlastpScoreMatrix::Adjusted(&adjusted),
    ] {
        for reverse in [false, true] {
            for (gap_open, gap_extend, x_drop) in [(11, 1, 25), (2, 0, 4), (0, 0, 2)] {
                for len in [128usize, 1, 32, 0, 7, 2] {
                    for variant in 0..4 {
                        let query: Vec<u8> = (0..=len)
                            .map(|i| [ncbistdaa::A, ncbistdaa::C, ncbistdaa::G][i % 3])
                            .collect();
                        let mut subject = query.clone();
                        if variant == 1 {
                            subject.splice(len / 2..len / 2, [ncbistdaa::W; 3]);
                        } else if variant == 2 && len > 3 {
                            subject.drain(len / 2..len / 2 + 2);
                        } else if variant == 3 && len > 3 {
                            subject[len / 2] = FENCE_SENTRY;
                        }
                        let run = |scratch: &mut GapAlignScratch| {
                            let mut fence = false;
                            let result = align_ex_protein_scalar(
                                &query,
                                &subject,
                                len,
                                subject.len() - 1,
                                score_matrix,
                                gap_open,
                                gap_extend,
                                x_drop,
                                reverse,
                                scratch,
                                &mut fence,
                                false,
                            );
                            // Empty inputs return before resetting existing scratch.
                            let used = if len == 0 { 0 } else { scratch.trace_rows_used };
                            (
                                result,
                                fence,
                                scratch.trace_offsets[..used].to_vec(),
                                scratch.trace_rows[..used].to_vec(),
                            )
                        };
                        let fresh = run(&mut GapAlignScratch::new());
                        let reused = run(&mut scratch);
                        assert_eq!(reused, fresh);
                        saw_offset |= reused.2.iter().any(|&offset| offset > 0);
                        saw_gap |= reused
                            .0
                             .3
                            .iter()
                            .any(|op| !matches!(op, GapEditOp::Sub(_)));
                        saw_fence |= reused.1;
                    }
                }
            }
        }
    }
    assert!(saw_offset && saw_gap && saw_fence);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
// The SIMD kernel with each matrix's tables (`x_standard_tables`) gives the scalar port's
// result, score-only and with traceback, in both directions: on 8-bit lanes (the largest
// X-drop they take, x + goe + 2M + 1 = 127 in `fits_8bit`), on 16-bit lanes (one more), and
// with a sentinel (code 0) in either sequence or NCBI's fence in the subject.
#[test]
fn x_standard_tables_give_the_scalar_port_for_every_matrix() {
    use crate::utils::xdrop_simd::{xdrop_align, ColSeq, RowSeq};
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut below = move |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    let mut scratch = GapAlignScratch::new();
    for matrix in ScoringMatrix::ALL {
        assert!(x_standard_tables(matrix).is_some(), "{matrix:?}");
        let top = (1..28u8)
            .flat_map(|q| (1..28u8).map(move |s| standard_score(matrix, q, s)))
            .max()
            .unwrap();
        let score_matrix = BlastpScoreMatrix::Standard(matrix);
        let mut fast_runs = 0;
        for (gap_open, gap_extend) in [(11, 1), (14, 2), (5, 2)] {
            let goe = gap_open + gap_extend;
            let x8 = 126 - goe - 2 * top;
            assert!(x8 >= goe, "{matrix:?}: x8 {x8}");
            for x_drop in [10, 38, x8, x8 + 1] {
                for variant in 0..4 {
                    let len1 = 20 + below(120) as usize;
                    let mut query: Vec<u8> = (0..=len1).map(|_| 1 + below(27) as u8).collect();
                    // the subject: the query with substitutions, insertions and deletions
                    let mut subject = vec![1 + below(27) as u8];
                    for &r in &query[1..] {
                        match below(20) {
                            0 => {}
                            1 => subject.extend([r, 1 + below(27) as u8]),
                            2..=4 => subject.push(1 + below(27) as u8),
                            _ => subject.push(r),
                        }
                    }
                    let len2 = subject.len() - 1;
                    match variant {
                        1 => subject[1 + below(len2 as u64) as usize] = 0,
                        2 => query[1 + below(len1 as u64) as usize] = 0,
                        3 => subject[1 + below(len2 as u64) as usize] = FENCE_SENTRY,
                        _ => {}
                    }
                    for reverse in [false, true] {
                        let at = format!(
                            "{matrix:?} {gap_open}/{gap_extend} x {x_drop} variant {variant} \
                                 reverse {reverse}"
                        );
                        let mut fence = (false, false);
                        let traceback = align_ex_protein(
                            &query,
                            &subject,
                            len1,
                            len2,
                            score_matrix,
                            gap_open,
                            gap_extend,
                            x_drop,
                            reverse,
                            &mut scratch,
                            &mut fence.0,
                            false,
                        );
                        let expect = align_ex_protein_scalar(
                            &query,
                            &subject,
                            len1,
                            len2,
                            score_matrix,
                            gap_open,
                            gap_extend,
                            x_drop,
                            reverse,
                            &mut GapAlignScratch::new(),
                            &mut fence.1,
                            false,
                        );
                        assert_eq!(traceback, expect, "{at}: traceback");
                        assert_eq!(fence.0, fence.1, "{at}: fence");
                        if variant == 3 {
                            // the fence is read by the traceback only
                            continue;
                        }
                        let score_only = align_ex_protein_score_only(
                            &query,
                            0,
                            &subject,
                            0,
                            len1,
                            len2,
                            score_matrix,
                            gap_open,
                            gap_extend,
                            x_drop,
                            reverse,
                            &mut scratch,
                        );
                        let expect = align_ex_protein_score_only_scalar(
                            &query,
                            0,
                            &subject,
                            0,
                            len1,
                            len2,
                            score_matrix,
                            gap_open,
                            gap_extend,
                            x_drop,
                            reverse,
                            &mut GapAlignScratch::new(),
                        );
                        assert_eq!(score_only, expect, "{at}: score only");
                    }
                    // the kernel itself takes the problem at the widest X-drops, so the
                    // calls above compared it with the scalar port
                    if crate::utils::xdrop_simd::KERNEL_COMPILED && variant == 0 && x_drop >= x8 {
                        let q = RowSeq::Bytes {
                            data: &query,
                            base: 0,
                            step: 1,
                        };
                        let s = ColSeq {
                            data: &subject,
                            base: 1,
                            step: 1,
                            zero_from: usize::MAX,
                        };
                        let taken = with_x_scores(score_matrix, |scores| {
                            xdrop_align(
                                &q,
                                &s,
                                scores,
                                len1,
                                len2,
                                gap_open,
                                gap_extend,
                                x_drop,
                                false,
                                false,
                                &mut scratch.fast,
                            )
                        });
                        assert!(taken.is_some(), "{matrix:?} x {x_drop}: kernel declined");
                        fast_runs += 1;
                    }
                }
            }
        }
        assert!(
            fast_runs == 6 || !crate::utils::xdrop_simd::KERNEL_COMPILED,
            "{matrix:?}: {fast_runs}"
        );
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:563-575,634-639
// Stop scanning the row immediately upon encountering `FENCE_SENTRY`.
// The fence cell is never written. Only the preceding initialized bytes
// may become visible when leaving the first main-band row.
#[test]
fn test_traceback_fence_exposes_only_initialized_prefix() {
    // No NCBI counterpart: the test reads scalar-kernel scratch rows, so it calls the scalar
    // port (`align_ex_protein_scalar`) directly; the SIMD kernel writes no such rows.
    for (reverse, capacity, gap_extend) in
        [(false, 48, 1), (true, 48, 1), (false, 8, 0), (true, 8, 0)]
    {
        for fence_b_index in [0usize, 1, 7] {
            let query = [ncbistdaa::A; 9];
            let mut subject = query;
            let fence_index = if reverse {
                7 - fence_b_index
            } else {
                fence_b_index + 1
            };
            subject[fence_index] = FENCE_SENTRY;
            let mut scratch = GapAlignScratch::new();
            scratch.trace_rows = vec![vec![u8::MAX; capacity]; 2];
            scratch.trace_offsets = vec![0; 2];
            let mut fence = false;
            let result = align_ex_protein_scalar(
                &query,
                &subject,
                8,
                8,
                BlastpScoreMatrix::standard(ScoringMatrix::Blosum62),
                11,
                gap_extend,
                100,
                reverse,
                &mut scratch,
                &mut fence,
                false,
            );
            assert!(fence);
            assert!(result.3.is_empty());
            assert_eq!(scratch.trace_rows_used, 2);
            assert_eq!(scratch.trace_rows[1].len(), fence_b_index);
            assert!(scratch.trace_rows[1].iter().all(|&byte| byte != u8::MAX));
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4620-4653,4670-4716
// The independent direct C oracle supplies all expected status, flag,
// score, signed endpoints and complete edit operations. Input starts vary
// independently, including 0, 1, the prior start2 boundary and the last
// residue; the subject fence visits every position in both directions.
#[test]
fn test_fence_signed_traceback_grid_matches_pinned_ncbi() {
    let expected = include_str!("../../../tests/unit/blastp_fence_signed_traceback_expected.tsv");
    let mut actual = String::from(expected.lines().next().unwrap());
    actual.push('\n');
    for line in expected.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        let qstart: usize = fields[0].parse().unwrap();
        let sstart: usize = fields[1].parse().unwrap();
        let position: i32 = fields[2].parse().unwrap();
        let query = [1u8; 6];
        let mut subject = [1u8; 6];
        if position >= 0 {
            subject[usize::try_from(position).unwrap()] = FENCE_SENTRY;
        }
        let mut scratch = GapAlignScratch::new();
        let mut fence = false;
        let result = blast_gapped_alignment_with_traceback_with_scratch(
            &query,
            &subject,
            qstart,
            sstart,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            18,
            &mut scratch,
            Some(&mut fence),
        )
        .expect("NCBI status is zero for every valid fence-grid input");
        let edits = if result.edit_script.is_empty() {
            "-".to_owned()
        } else {
            result
                .edit_script
                .iter()
                .map(|op| {
                    let kind = match op {
                        GapEditOp::Sub(_) => 3,
                        GapEditOp::Del(_) => 0,
                        GapEditOp::Ins(_) => 6,
                    };
                    format!("{kind}:{}", op.num())
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        actual.push_str(&format!(
            "{qstart}\t{sstart}\t{position}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{edits}\n",
            i32::from(fence),
            result.score,
            result.query_start,
            result.query_stop,
            result.subject_start,
            result.subject_stop
        ));
    }
    if let Some(path) = std::env::var_os("LOSAT_FENCE_BOUNDARY_DUMP") {
        std::fs::write(path, &actual).unwrap();
    }
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:802-821
// A huge X-drop (TN-4: `-xdrop_gap 5e8`) keeps NCBI's cell count, but only the
// cells up to `N` exist, and the alignments equal those of an X-drop that the
// extension never reaches.
#[test]
fn huge_xdrop_keeps_ncbi_cell_count_and_touches_only_cells_up_to_n() {
    let encode =
        |text: &[u8]| -> Vec<u8> { text.iter().map(|&c| aa_char_to_ncbistdaa(c)).collect() };
    // Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
    let query = encode(b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ");
    let subject = encode(b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQGGAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQWW");
    let huge = 1_000_000_000;
    let run_traceback = |x_drop: i32| {
        let mut scratch = GapAlignScratch::new();
        let result = blast_gapped_alignment_with_traceback_with_scratch(
            &query,
            &subject,
            20,
            20,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            x_drop,
            &mut scratch,
            None,
        )
        .expect("alignment");
        (
            (
                result.query_start,
                result.query_stop,
                result.subject_start,
                result.subject_stop,
                result.score,
                result.edit_script,
            ),
            scratch.dp_mem.len(),
            scratch.dp_mem_alloc,
        )
    };
    // the comparison side: an X-drop the 16-bit lanes hold as it is (no clamp) and that the
    // extension never reaches (it drops at most about 900 here), so a clamp to a wrong `R`
    // on the huge side would show
    let (bounded, bounded_cells, _) = run_traceback(20_000);
    let (unbounded, cells, alloc) = run_traceback(huge);
    assert_eq!(unbounded, bounded);
    // `GapAlignScratch::new` holds NCBI's first 1000 cells, more than `N + 1`.
    assert_eq!(bounded_cells, 1000);
    assert_eq!(cells, 1000);
    // No NCBI counterpart: the vector kernels run a huge X-drop as the clamp's `R`
    // (`xdrop_simd::x_drop_for_lanes`, `LOSAT_X_DP32`) and do not touch `dp_mem_alloc`; the
    // scalar port, which answers without them (other targets, `LOSAT_X_DP32=0`) and for the
    // calls they decline, keeps NCBI's count (below).
    let vector = crate::utils::xdrop_simd::x_kernels_run() && crate::utils::xdrop_simd::x_dp32();
    if vector {
        assert_eq!(alloc, 1000);
    } else {
        assert!(alloc >= huge as usize + 3, "{alloc}");
    }
    // The scalar ports themselves with the huge X-drop (the left extensions of the traceback
    // and of the score-only alignment): NCBI's count, only the cells up to `N`, and the
    // results of the X-drop the extension never reaches.
    let scalar_traceback = |x_drop: i32| {
        let mut scratch = GapAlignScratch::new();
        let mut fence_hit = false;
        let result = align_ex_protein_scalar(
            &query,
            &subject,
            21,
            21,
            BlastpScoreMatrix::Standard(ScoringMatrix::Blosum62),
            11,
            1,
            x_drop,
            true,
            &mut scratch,
            &mut fence_hit,
            false,
        );
        (result, scratch.dp_mem.len(), scratch.dp_mem_alloc)
    };
    let (bounded, _, _) = scalar_traceback(100_000);
    let (unbounded, cells, alloc) = scalar_traceback(huge);
    assert_eq!(unbounded, bounded);
    assert_eq!(cells, 1000);
    assert!(alloc >= huge as usize + 3, "{alloc}");
    let scalar_score_only = |x_drop: i32| {
        let mut scratch = GapAlignScratch::new();
        let result = align_ex_protein_score_only_scalar(
            &query,
            0,
            &subject,
            0,
            21,
            21,
            BlastpScoreMatrix::Standard(ScoringMatrix::Blosum62),
            11,
            1,
            x_drop,
            true,
            &mut scratch,
        );
        (result, scratch.dp_mem.len(), scratch.dp_mem_alloc)
    };
    let (bounded, _, _) = scalar_score_only(100_000);
    let (unbounded, cells, alloc) = scalar_score_only(huge);
    assert_eq!(unbounded, bounded);
    assert_eq!(cells, 1000);
    assert!(alloc >= huge as usize + 3, "{alloc}");
    for mode in [
        BlastpGappedAlignmentMode::Exact,
        BlastpGappedAlignmentMode::Restricted,
    ] {
        let mut scratch = GapAlignScratch::new();
        let unbounded = blastp_score_only_gapped_alignment_with_scratch(
            &query,
            &subject,
            20,
            20,
            ScoringMatrix::Blosum62,
            11,
            1,
            huge,
            mode,
            &mut scratch,
        );
        assert_eq!(scratch.dp_mem.len(), 1000);
        // the vector kernels (the exact mode's dispatch, `xdrop_align`) leave the count as
        // it was, as above
        if vector && matches!(mode, BlastpGappedAlignmentMode::Exact) {
            assert_eq!(scratch.dp_mem_alloc, 1000);
        } else {
            assert!(scratch.dp_mem_alloc >= huge as usize + 3);
        }
        let bounded = blastp_score_only_gapped_alignment(
            &query,
            &subject,
            20,
            20,
            ScoringMatrix::Blosum62,
            11,
            1,
            20_000,
            mode,
        );
        assert_eq!(unbounded, bounded);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:736-741,766-770
// Both sides of `LOSAT_X_DP32` in one process: with X-drops past the 16-bit lanes (50,000, 1e6
// and 2^31 - 300) on a composition-adjusted matrix (BLOSUM62 x 150, the scale of the redo's
// matrices), the dispatched traceback and score-only extensions (`align_ex_protein`,
// `align_ex_protein_score_only`: the 32-bit lanes) give what the scalar ports give
// (`align_ex_protein_scalar`, `align_ex_protein_score_only_scalar`: what `LOSAT_X_DP32=0`
// runs), left and right of the seed as the extension wrappers call them: offsets, score, edit
// script and fence flag.
#[test]
fn wide_lanes_give_the_scalar_ports_results() {
    let mut scores = [[0i32; 28]; 28];
    for (q, row) in scores.iter_mut().enumerate() {
        for (s, v) in row.iter_mut().enumerate() {
            let b = blosum62_score_ncbistdaa_direct(q as u8, s as u8);
            // the gap letter keeps BLAST_SCORE_MIN, as the adjusted matrices do
            *v = if b > -100 { 150 * b } else { b };
        }
    }
    let adjusted = AdjustedProteinMatrix { scores };
    let matrix = BlastpScoreMatrix::Adjusted(&adjusted);
    let (gap_open, gap_extend) = (11 * 150, 150);
    let mut state = 0x2545_F491u32;
    let mut next = move |n: u32| {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        (state >> 8) % n
    };
    let runs = crate::utils::xdrop_simd::X32_RUNS.with(|n| n.get());
    for case in 0..10 {
        let len = 60 + next(240) as usize;
        // letters 1..=24 (no gap letter, whose -32768 entries the 32-bit lanes decline)
        let query: Vec<u8> = (0..len).map(|_| 1 + next(24) as u8).collect();
        let mut subject = Vec::with_capacity(len + 8);
        for (i, &c) in query.iter().enumerate() {
            match next(40) {
                0 => subject.extend([c, 1 + next(24) as u8]),
                1 if i > 0 => {}
                2..=6 => subject.push(1 + next(24) as u8),
                _ => subject.push(c),
            }
        }
        let (q_start, s_start) = (len / 2, (len / 2).min(subject.len() - 2));
        for x_drop in [50_000, 1_000_000, i32::MAX - 300] {
            for reverse in [true, false] {
                let (q, s, len1, len2) = if reverse {
                    (&query[..], &subject[..], q_start + 1, s_start + 1)
                } else {
                    (
                        &query[q_start..],
                        &subject[s_start..],
                        query.len() - q_start - 1,
                        subject.len() - s_start - 1,
                    )
                };
                let (mut a, mut b) = (GapAlignScratch::new(), GapAlignScratch::new());
                let (mut fence_a, mut fence_b) = (false, false);
                let got = align_ex_protein(
                    q,
                    s,
                    len1,
                    len2,
                    matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    reverse,
                    &mut a,
                    &mut fence_a,
                    false,
                );
                let want = align_ex_protein_scalar(
                    q,
                    s,
                    len1,
                    len2,
                    matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    reverse,
                    &mut b,
                    &mut fence_b,
                    false,
                );
                assert_eq!(
                    (got, fence_a),
                    (want, fence_b),
                    "case {case} x {x_drop} reverse {reverse}: traceback"
                );
                let (q_base, s_base) = if reverse { (0, 0) } else { (q_start, s_start) };
                let got = align_ex_protein_score_only(
                    &query, q_base, &subject, s_base, len1, len2, matrix, gap_open, gap_extend,
                    x_drop, reverse, &mut a,
                );
                let want = align_ex_protein_score_only_scalar(
                    &query, q_base, &subject, s_base, len1, len2, matrix, gap_open, gap_extend,
                    x_drop, reverse, &mut b,
                );
                assert_eq!(
                    got, want,
                    "case {case} x {x_drop} reverse {reverse}: score only"
                );
            }
        }
    }
    if crate::utils::xdrop_simd::x_kernels_run() && crate::utils::xdrop_simd::x_dp32() {
        let runs = crate::utils::xdrop_simd::X32_RUNS.with(|n| n.get()) - runs;
        assert!(runs > 100, "the 32-bit lanes ran {runs} times");
    }
}
