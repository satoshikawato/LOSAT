use super::*;

fn record(seq: &[u8]) -> FastaRecord {
    FastaRecord::new("Query_1", b"q", seq)
}

fn hsp(frame: i32, q: (usize, usize), s: (usize, usize), ops: Option<Vec<GapEditOp>>) -> BlastnHsp {
    BlastnHsp {
        identity: 0.0,
        length: 0,
        mismatch: 0,
        gapopen: 0,
        q_start: 0,
        q_end: 0,
        s_start: 0,
        s_end: 0,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: frame,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
        raw_score: 0,
        internal_q_offset_0: q.0 as i32,
        internal_q_end_0: q.1 as i32,
        internal_s_offset_0: s.0,
        internal_s_end_0: s.1,
        internal_query_context_offset: 0,
        gap_info: ops,
        num_positives: 0,
    }
}

// seqinfosrc_seqdb.cpp:121-133 (cited at `subject_mask_regions`): a stored mask is shown
// when it holds a letter and meets the closed range `[offset, end]` of an HSP (the end one
// past the HSP); overlapping and nested masks (hand-made masks of the M4a-2 edge databases)
// stay as they are, in their order (`x_FillLocList` places each). A subject set's merged
// lowercase is kept as it is.
#[test]
fn database_masks_meeting_an_hsp_are_shown() {
    let masks = [
        MaskedInterval::new(600, 700),
        MaskedInterval::new(100, 200),
        MaskedInterval::new(150, 250),
        MaskedInterval::new(250, 260),
        MaskedInterval::new(300, 300),
        MaskedInterval::new(900, 950),
    ];
    // The mask 600-700 holds 600..=699: it meets an HSP from 699 on, not one from 700 on.
    let hsps = [
        hsp(1, (0, 10), (199, 250), None),
        hsp(1, (0, 10), (700, 760), None),
    ];
    let shown = |masks: &[MaskedInterval], hsps: &[BlastnHsp]| {
        subject_mask_regions(masks, hsps).collect::<Vec<_>>()
    };
    assert!(shown(&masks[..1], &hsps).is_empty());
    let hsps = [
        hsp(1, (0, 10), (199, 250), None),
        hsp(1, (0, 10), (699, 701), None),
    ];
    // 100-200, 150-250 and 250-260 meet the first HSP, 600-700 the second; the empty 300-300
    // and 900-950 none.
    assert_eq!(
        shown(&masks, &hsps),
        [
            MaskedInterval::new(600, 700),
            MaskedInterval::new(100, 200),
            MaskedInterval::new(150, 250),
            MaskedInterval::new(250, 260)
        ]
    );
    // Merged lowercase runs are kept as they are.
    let runs = merged(&masks);
    assert_eq!(
        shown(&runs, &hsps),
        [MaskedInterval::new(100, 260), MaskedInterval::new(600, 700)]
    );
    let hsps = [hsp(1, (0, 10), (940, 960), None)];
    assert_eq!(shown(&runs, &hsps), [MaskedInterval::new(900, 950)]);
}

#[test]
fn plus_strand_rows_follow_the_edit_script() {
    let query = record(b"ACGTTACG");
    let subject = record(b"ACGTACCG");
    let ops = vec![GapEditOp::Sub(4), GapEditOp::Ins(1), GapEditOp::Sub(3)];
    let (q, s) = displayed_rows(
        &hsp(1, (0, 8), (0, 7), Some(ops)),
        query.seq(),
        subject.seq().len(),
        |from, to| subject.seq()[from..to].to_ascii_uppercase(),
    )
    .unwrap();
    assert_eq!(
        (q.as_slice(), s.as_slice()),
        (&b"ACGTTACG"[..], &b"ACGT-ACC"[..])
    );
}

// A minus-strand HSP is shown with the query on its plus strand: both rows of the
// context alignment are reverse-complemented, gaps included.
#[test]
fn minus_strand_rows_are_reverse_complemented() {
    let query = record(b"AACCGGTT");
    let subject = record(b"CCGGT");
    // Context (reverse complement of the query) AACCGGTT[1..6] = ACCGG against CCGG-.
    let ops = vec![GapEditOp::Ins(1), GapEditOp::Sub(4)];
    let (q, s) = displayed_rows(
        &hsp(-1, (1, 6), (0, 4), Some(ops)),
        query.seq(),
        subject.seq().len(),
        |from, to| subject.seq()[from..to].to_ascii_uppercase(),
    )
    .unwrap();
    assert_eq!((q.as_slice(), s.as_slice()), (&b"CCGGT"[..], &b"CCGG-"[..]));
}

// The XML rows are written into buffers that the caller reuses: stale contents are cleared
// and the minus-strand reversal in place gives the rows of `displayed_rows`.
#[test]
fn rows_into_reused_buffers_equal_displayed_rows() {
    let query = record(b"AACCGGTT");
    let subject = record(b"CCGGT");
    let (mut q, mut s) = (b"stale row".to_vec(), b"x".to_vec());
    for (frame, window, ops) in [
        (-1, (1, 6), vec![GapEditOp::Ins(1), GapEditOp::Sub(4)]),
        (1, (0, 5), vec![GapEditOp::Sub(5)]),
    ] {
        let hsp = hsp(frame, window, (0, 4 + usize::from(frame > 0)), Some(ops));
        let letters = |from: usize, to: usize| subject.seq()[from..to].to_ascii_uppercase();
        let once = displayed_rows(&hsp, query.seq(), subject.seq().len(), letters).unwrap();
        displayed_rows_into(
            &hsp,
            query.seq(),
            subject.seq().len(),
            letters,
            &mut q,
            &mut s,
        )
        .unwrap();
        assert_eq!((q.clone(), s.clone()), once);
    }
}

#[test]
fn masked_positions_are_lowercase_in_the_row_direction() {
    let mut row = b"AC-GT".to_vec();
    let mut masks = crate::report::query_anchored::HspMasks::default();
    let shown = subject_mask_seek([MaskedInterval { start: 7, end: 9 }].into_iter());
    set_subject_row_masks(&mut masks, &row, 10, -1, &shown);
    masks.lower_row(1, &mut row);
    assert_eq!(row, b"Ac-gT");
}

// A context in BLASTNA (either letter case, no uppercase copy) gives the block of the
// letters it was encoded from.
#[test]
fn karlin_of_blastna_equals_karlin_of_letters() {
    for seq in [
        &b"ACGTACGTTTGA"[..],
        b"acgtnACGTRYKMSWBDHVnn-",
        b"NNNN",
        b"",
        b"aaaaCCCCggggTTTTUu",
    ] {
        let blastna = encode_iupac_to_blastna(seq);
        let from_letters = query_ungapped_karlin(seq, 2, -3);
        let from_blastna = query_ungapped_karlin_blastna(&blastna, 2, -3);
        match (from_letters, from_blastna) {
            (None, None) => {}
            (Some(a), Some(b)) => {
                assert_eq!(a.lambda.to_bits(), b.lambda.to_bits());
                assert_eq!(a.k.to_bits(), b.k.to_bits());
                assert_eq!(a.h.to_bits(), b.h.to_bits());
            }
            other => panic!("{seq:?}: {other:?}"),
        }
    }
}

#[test]
fn an_all_n_query_has_no_ungapped_block() {
    assert!(query_ungapped_karlin(b"NNNN", 2, -3).is_none());
    let acgt = query_ungapped_karlin(b"ACGTACGTTTGA", 2, -3).unwrap();
    assert!((acgt.lambda - 0.634).abs() < 5e-4);
}
