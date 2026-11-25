use super::*;
use crate::algorithm::tblastn::search_gapped::GappedHsp;
use crate::algorithm::tblastn::stage_d_linking::LinkedHsp;
use crate::algorithm::tblastn::stage_d_results::KappaHspPayload;
use crate::common::GapEditOp::{Del, Ins, Sub};
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;

fn sample_nullb_hsp(frame: i8, edit_script: Vec<GapEditOp>) -> (LinkedHsp, KappaHspPayload) {
    let align_length: usize = edit_script.iter().map(|op| op.num() as usize).sum();
    let gap_letters: usize = edit_script
        .iter()
        .filter(|op| !op.is_sub())
        .map(|op| op.num() as usize)
        .sum();
    let gap_opens = edit_script.iter().filter(|op| !op.is_sub()).count();
    let linked = LinkedHsp {
        hsp: GappedHsp {
            frame,
            score: 86,
            q_start: 10,
            q_end: 14,
            s_start: 14,
            s_end: 18,
            q_gapped_start: 15,
            s_gapped_start: 19,
        },
        context: 0,
        evalue: 4.1e-12,
        num: 1,
        source_index: 0,
    };
    let payload = KappaHspPayload {
        bit_score: 37.6,
        num_ident: 2,
        num_positives: 2,
        report_num_ident: 2,
        report_mismatches: 1,
        align_length,
        mismatches: 1,
        gap_opens,
        gap_letters,
        edit_script,
        matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
        query_left_out: 0,
        nullb: None,
    };
    (linked, payload)
}

#[test]
fn positive_frame_db_rem0_pads_nullb_with_dash() {
    let (mut linked, mut payload) = sample_nullb_hsp(1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut linked, &mut payload, 57, false));
    let n = payload.nullb.expect("nullb report");
    assert_eq!(n.blast_hsp_q_span, 4);
    assert!(!n.scores_invalidated);
    assert!(n.trailing_nullb_dash);
    assert!(!n.translate_na_failure);
    assert!(!n.std_seg_lengths_inaccurate);
    assert_eq!((linked.hsp.q_start, linked.hsp.q_end), (10, 16));
    assert_eq!((linked.hsp.s_start, linked.hsp.s_end), (14, 20));
    assert_eq!(linked.hsp.score, 86);
}

#[test]
fn positive_frame_db_rem2_fails_translation() {
    let (mut linked, mut payload) = sample_nullb_hsp(1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut linked, &mut payload, 59, false));
    let n = payload.nullb.expect("nullb report");
    assert!(n.translate_na_failure);
    assert!(!n.std_seg_lengths_inaccurate);
    assert!(!n.scores_invalidated);
    assert_eq!((linked.hsp.q_start, linked.hsp.q_end), (10, 16));
    assert_eq!((linked.hsp.s_start, linked.hsp.s_end), (14, 20));
}

#[test]
fn positive_frame_subject_rem0_truncates_and_invalidates_scores() {
    let (mut linked, mut payload) = sample_nullb_hsp(1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut linked, &mut payload, 57, true));
    let n = payload.nullb.expect("nullb report");
    assert!(n.scores_invalidated);
    assert_eq!(linked.hsp.score, -1);
    assert_eq!(payload.bit_score, -1.0);
    assert_eq!(linked.evalue, -1.0);
    assert_eq!(payload.edit_script, vec![Sub(1), Ins(3), Del(3), Sub(1)]);
    assert_eq!(payload.align_length, 8);
    assert_eq!((linked.hsp.q_start, linked.hsp.q_end), (10, 15));
    assert_eq!((linked.hsp.s_start, linked.hsp.s_end), (14, 19));
}

#[test]
fn positive_frame_subject_rem2_k2_sets_fuzz_and_inaccurate_lengths() {
    let (mut linked, mut payload) = sample_nullb_hsp(1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut linked, &mut payload, 59, true));
    let n = payload.nullb.expect("nullb report");
    assert!(n.scores_invalidated);
    assert!(n.std_seg_lengths_inaccurate);
    assert!(!n.translate_na_failure);
    assert_eq!(n.last_seg_fuzz, Some((58, 1)));
    assert_eq!((linked.hsp.q_start, linked.hsp.q_end), (10, 15));
}

#[test]
fn positive_frame_subject_rem1_k1_sets_fuzz_and_translate_failure() {
    let (mut linked, mut payload) = sample_nullb_hsp(1, vec![Sub(1), Ins(4), Del(4), Sub(1)]);
    assert!(prepare_nullb_hsp(&mut linked, &mut payload, 58, true));
    let n = payload.nullb.expect("nullb report");
    assert!(n.scores_invalidated);
    assert!(!n.std_seg_lengths_inaccurate);
    assert!(n.translate_na_failure);
    assert_eq!(n.last_seg_fuzz, Some((57, 0)));
    assert_eq!((linked.hsp.q_start, linked.hsp.q_end), (10, 16));
}

#[test]
fn negative_frame_db_rem0_wraps_tail_and_rem2_turns_sub_into_ins() {
    let (mut l0, mut p0) = sample_nullb_hsp(-1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut l0, &mut p0, 57, false));
    assert_eq!(p0.nullb.expect("nullb").wrapped_minus_tail, 2);
    assert_eq!(p0.edit_script, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert_eq!(p0.gap_letters, 8);
    assert_eq!(p0.gap_opens, 3);
    assert_eq!((l0.hsp.s_start, l0.hsp.s_end), (14, 18));

    let (mut l2, mut p2) = sample_nullb_hsp(-1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut l2, &mut p2, 59, false));
    assert_eq!(p2.nullb.expect("nullb").wrapped_minus_tail, 0);
    assert_eq!(p2.edit_script, vec![Sub(1), Ins(3), Del(3), Ins(2)]);
    assert_eq!(p2.gap_letters, 8);
    assert_eq!(p2.gap_opens, 3);
    assert_eq!((l2.hsp.s_start, l2.hsp.s_end), (14, 18));
}

#[test]
fn negative_frame_subject_drops_last_sub_and_keeps_scores() {
    let (mut linked, mut payload) = sample_nullb_hsp(-1, vec![Sub(1), Ins(3), Del(3), Sub(2)]);
    assert!(prepare_nullb_hsp(&mut linked, &mut payload, 57, true));
    let n = payload.nullb.expect("nullb");
    assert!(!n.scores_invalidated);
    assert_eq!(linked.hsp.score, 86);
    assert_eq!(payload.edit_script, vec![Sub(1), Ins(3), Del(3)]);
    assert_eq!(payload.align_length, 7);
    assert_eq!((linked.hsp.q_start, linked.hsp.q_end), (10, 14));
    assert_eq!((linked.hsp.s_start, linked.hsp.s_end), (14, 18));
}
