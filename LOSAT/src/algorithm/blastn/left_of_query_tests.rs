use super::*;
use crate::blastinput::app::{FormatChoice, ReportFormat};
use crate::report::archive::options::RequestKind;
use crate::report::tabular::TabularSpec;

fn spec(fields: &str) -> TabularSpec {
    TabularSpec::of_app(
        ReportFormat::Tabular,
        &FormatChoice {
            number: 6,
            spec: fields.to_string(),
            delimiter: String::new(),
        },
        RequestKind::Blastn,
    )
    .expect("valid custom spec")
}

#[test]
fn tabular_exception_matches_ncbi_field_classification() {
    // Fields outside opens_counts do not throw.
    assert!(tabular_exception(&spec("qseqid sseqid evalue bitscore score"), 200).is_none());
    assert!(tabular_exception(&spec("qseqid sseqid qcovs qcovhsp qcovus"), 200).is_none());
    assert!(tabular_exception(&spec("qseqid sseqid frames"), 200).is_none());
    assert!(tabular_exception(&spec("qseqid sseqid nident"), 0).is_none());

    // Subject-only fields inside opens_counts do not throw.
    assert!(tabular_exception(&spec("qseqid sseqid sstart send sstrand sframe"), 200).is_none());
    assert!(tabular_exception(&spec("qseqid sseqid nident"), 200).is_none());
    assert!(tabular_exception(&spec("qseqid sseqid nident sframe"), 200).is_none());

    let check = |field: &str, line: u32, prefix: &[u8]| {
        let exc = tabular_exception(&spec(&format!("qseqid sseqid {field}")), 200).expect(field);
        assert_eq!(exc.site().line, line, "field {field}");
        assert!(exc.get_msg().starts_with(prefix), "field {field}");
    };

    // Length/gap/mismatch/percent fields throw x_GetSeqRightSeg (alnmap.cpp:753).
    for f in ["length", "gaps", "mismatch", "ppos", "pident", "gapopen"] {
        check(f, 753, b"CAlnVec::x_GetSeqRightSeg()");
    }

    // Query/row/frame fields inside opens_counts throw x_GetSeqLeftSeg (alnmap.cpp:732).
    for f in [
        "qstart",
        "qend",
        "qframe",
        "qseq",
        "sseq",
        "positive",
        "btop",
        "sstart frames",
        "nident frames",
    ] {
        check(f, 732, b"CAlnVec::x_GetSeqLeftSeg()");
    }

    // Standard outfmt 6 throws alnmap.cpp:753.
    let std_spec = TabularSpec::std(false);
    let exc = tabular_exception(std_spec, 200).expect("std");
    assert_eq!(exc.site().line, 753);
}

#[test]
fn non_empty_record_query_range_handles_wrapping() {
    let mut hsp = BlastnHsp {
        identity: 100.0,
        length: 200,
        mismatch: 0,
        gapopen: 0,
        q_start: 0,
        q_end: 0,
        s_start: 1,
        s_end: 200,
        e_value: 0.0,
        bit_score: 367.0,
        num_ident: 200,
        query_frame: 1,
        query_length: 300,
        q_idx: 0,
        s_idx: 0,
        raw_score: 400,
        internal_q_offset_0: -702,
        internal_q_end_0: -502,
        internal_s_offset_0: 0,
        internal_s_end_0: 200,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: 0,
    };
    // -702..=-503 does not cross 0 in u32: non-empty with length 200.
    assert_eq!(
        non_empty_record_query_range(&hsp, 0, true),
        Some((4294966594, 4294966793))
    );
    // Shifted by 701 -> start = -1, len = 200 -> wraps across 0 in u32: empty!
    assert_eq!(non_empty_record_query_range(&hsp, 701, true), None);
    // Shifted by 702 -> start = 0, len = 200: normal 0..=199.
    assert_eq!(
        non_empty_record_query_range(&hsp, 702, true),
        Some((0, 199))
    );
    // Exactly touching u32::MAX (start = -200, len = 200 -> to = u32::MAX, to_open wraps to 0): empty!
    hsp.internal_q_offset_0 = -200;
    hsp.internal_q_end_0 = 0;
    assert_eq!(non_empty_record_query_range(&hsp, 0, true), None);
}

#[test]
fn blastn_hsp_scores_handles_left_of_record_and_wrapping() {
    let hsp = BlastnHsp {
        identity: 100.0,
        length: 200,
        mismatch: 0,
        gapopen: 0,
        q_start: 0,
        q_end: 0,
        s_start: 10,
        s_end: 209,
        e_value: 1e-50,
        bit_score: 367.0,
        num_ident: 200,
        query_frame: 1,
        query_length: 300,
        q_idx: 0,
        s_idx: 0,
        raw_score: 400,
        internal_q_offset_0: -1202,
        internal_q_end_0: -1002,
        internal_s_offset_0: 9,
        internal_s_end_0: 209,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: 0,
    };
    let s0 = blastn_hsp_scores(&hsp, 0, true);
    assert_eq!(s0.query_from, u64::from((-1202i32) as u32));
    assert_eq!(s0.query_to, u64::from((-1003i32) as u32));
    assert_eq!(s0.query_sort_start, -1202);
    assert_eq!(s0.subject_from, 9);

    // Shifting by 1201 yields start = -1, len = 200: range wraps across 0 -> empty (1, 0),
    // while signed sort start is min(u32::MAX, 198) as i32 == 198.
    let s1 = blastn_hsp_scores(&hsp, 1201, true);
    assert_eq!((s1.query_from, s1.query_to), (1, 0));
    assert_eq!(s1.query_sort_start, 198);
}

#[test]
fn gapped_insert_subject_scores_throw_matches_qcovs_and_qcovus() {
    use crate::algorithm::blastn::hsp::BlastnHspList;

    let make = |q_off: i32, q_end: i32, s_start: usize, s_end: usize| BlastnHsp {
        identity: 100.0,
        length: (q_end - q_off) as usize,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 100,
        s_start,
        s_end,
        e_value: 1e-20,
        bit_score: 200.0,
        num_ident: (q_end - q_off) as usize,
        query_frame: 1,
        query_length: 300,
        q_idx: 0,
        s_idx: 0,
        raw_score: 200,
        internal_q_offset_0: q_off,
        internal_q_end_0: q_end,
        internal_s_offset_0: s_start - 1,
        internal_s_end_0: s_end,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: 0,
    };

    // Normal HSP first covering subject 1..=200, then class N HSP overlapping subject 50..=150:
    // qcovs checks every HSP -> throws Dense_seg.cpp:235;
    // qcovus skips the subject-covered second HSP -> does not throw.
    let mut list_covered = BlastnHitList::new(10);
    list_covered.update(BlastnHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![make(0, 200, 1, 200), make(-500, -400, 50, 150)],
        best_evalue: 1e-20,
    });
    let exc = gapped_insert_subject_scores_throw(&list_covered, 0, 300, true, false).unwrap();
    assert_eq!(exc.site().line, 235);
    assert!(gapped_insert_subject_scores_throw(&list_covered, 0, 300, false, true).is_none());
    assert!(gapped_insert_subject_scores_throw(&list_covered, 0, 0, true, true).is_none());

    // Disjoint second HSP with class N throws for qcovus too.
    let mut list_disjoint = BlastnHitList::new(10);
    list_disjoint.update(BlastnHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![make(0, 200, 1, 200), make(-500, -400, 300, 400)],
        best_evalue: 1e-20,
    });
    let exc2 = gapped_insert_subject_scores_throw(&list_disjoint, 0, 300, false, true).unwrap();
    assert_eq!(exc2.site().line, 235);
}
