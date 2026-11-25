use super::*;
use crate::algorithm::tblastn::search_gapped::GappedHsp;
use crate::algorithm::tblastn::stage_d_linking::LinkedHsp;
use crate::algorithm::tblastn::stage_d_results::KappaHspPayload;
use crate::common::GapEditOp::{Del, Ins, Sub};
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;

fn run_prepare(
    frame: i8,
    edit_script: Vec<GapEditOp>,
    subject_len: usize,
    remaps: bool,
) -> (LinkedHsp, KappaHspPayload, NullbReport) {
    let align_length: usize = edit_script.iter().map(|op| op.num() as usize).sum();
    let gap_letters: usize = edit_script
        .iter()
        .filter(|op| !op.is_sub())
        .map(|op| op.num() as usize)
        .sum();
    let gap_opens = edit_script.iter().filter(|op| !op.is_sub()).count();
    let mut linked = LinkedHsp {
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
    let mut payload = KappaHspPayload {
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
    assert!(prepare_nullb_hsp(
        &mut linked,
        &mut payload,
        subject_len,
        remaps
    ));
    let n = payload.nullb.expect("nullb report");
    (linked, payload, n)
}

fn k2_script() -> Vec<GapEditOp> {
    vec![Sub(1), Ins(3), Del(3), Sub(2)]
}

#[test]
fn positive_frame_db_rem0_and_rem2() {
    let (l0, _, n0) = run_prepare(1, k2_script(), 57, false);
    assert_eq!(n0.blast_hsp_q_span, 4);
    assert!(!n0.scores_invalidated);
    assert!(n0.trailing_nullb_dash);
    assert!(!n0.translate_na_failure);
    assert!(!n0.std_seg_lengths_inaccurate);
    assert_eq!(
        (l0.hsp.q_start, l0.hsp.q_end, l0.hsp.s_start, l0.hsp.s_end),
        (10, 16, 14, 20)
    );
    assert_eq!(l0.hsp.score, 86);

    let (l2, _, n2) = run_prepare(1, k2_script(), 59, false);
    assert!(n2.translate_na_failure);
    assert!(!n2.std_seg_lengths_inaccurate);
    assert!(!n2.scores_invalidated);
    assert_eq!(
        (l2.hsp.q_start, l2.hsp.q_end, l2.hsp.s_start, l2.hsp.s_end),
        (10, 16, 14, 20)
    );
}

#[test]
fn positive_frame_subject_rem0_rem2_k2_and_rem1_k1() {
    let (l0, p0, n0) = run_prepare(1, k2_script(), 57, true);
    assert!(n0.scores_invalidated);
    assert_eq!((l0.hsp.score, p0.bit_score, l0.evalue), (-1, -1.0, -1.0));
    assert_eq!(p0.edit_script, vec![Sub(1), Ins(3), Del(3), Sub(1)]);
    assert_eq!(p0.align_length, 8);
    assert_eq!((l0.hsp.q_end, l0.hsp.s_end), (15, 19));

    let (l2, _, n2) = run_prepare(1, k2_script(), 59, true);
    assert!(n2.scores_invalidated);
    assert!(n2.std_seg_lengths_inaccurate);
    assert!(!n2.translate_na_failure);
    assert_eq!(n2.last_seg_fuzz, Some((58, 1)));
    assert_eq!(l2.hsp.q_end, 15);

    let (l1, _, n1) = run_prepare(1, vec![Sub(1), Ins(4), Del(4), Sub(1)], 58, true);
    assert!(n1.scores_invalidated);
    assert!(!n1.std_seg_lengths_inaccurate);
    assert!(n1.translate_na_failure);
    assert_eq!(n1.last_seg_fuzz, Some((57, 0)));
    assert_eq!(l1.hsp.q_end, 16);
}

#[test]
fn negative_frame_db_and_subject() {
    let (l0, p0, n0) = run_prepare(-1, k2_script(), 57, false);
    assert_eq!(n0.wrapped_minus_tail, 2);
    assert_eq!(p0.edit_script, k2_script());
    assert_eq!((p0.gap_letters, p0.gap_opens), (8, 3));
    assert_eq!((l0.hsp.s_start, l0.hsp.s_end), (14, 18));

    let (l2, p2, n2) = run_prepare(-1, k2_script(), 59, false);
    assert_eq!(n2.wrapped_minus_tail, 0);
    assert_eq!(p2.edit_script, vec![Sub(1), Ins(3), Del(3), Ins(2)]);
    assert_eq!((p2.gap_letters, p2.gap_opens), (8, 3));
    assert_eq!((l2.hsp.s_start, l2.hsp.s_end), (14, 18));

    let (ls, ps, ns) = run_prepare(-1, k2_script(), 57, true);
    assert!(!ns.scores_invalidated);
    assert_eq!(ls.hsp.score, 86);
    assert_eq!(ps.edit_script, vec![Sub(1), Ins(3), Del(3)]);
    assert_eq!(ps.align_length, 7);
    assert_eq!((ls.hsp.q_end, ls.hsp.s_end), (14, 18));
}
