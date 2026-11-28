use super::*;

fn record(seq: &[u8]) -> FastaRecord {
    FastaRecord::new("Query_1", b"q", seq)
}

/// The letters of a query searched whole (no `-query_loc`).
fn whole(letters: &[u8]) -> QueryRecord<'_> {
    QueryRecord { letters, offset: 0 }
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
        // The context's length, as the engine sets it: the tests' queries have 8 letters.
        query_length: 8,
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
        whole(query.seq()),
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
        whole(query.seq()),
        subject.seq().len(),
        |from, to| subject.seq()[from..to].to_ascii_uppercase(),
    )
    .unwrap();
    assert_eq!((q.as_slice(), s.as_slice()), (&b"CCGGT"[..], &b"CCGG-"[..]));
}

// remaining-gaps X1b-p: the query letters of an HSP that NCBI keeps outside its query context
// with its Seq-align's query start at or past 0 (`context_letters`, alnvec.cpp:279-300): a
// column whose position on the query's plus strand lies past the query's end shows the gap
// character, the others their letters. Plus strand: offsets past the end; minus strand:
// offsets below 0 (positions past the end of the plus strand); a range that crosses the end
// shows the letters before it. A start below 0 is refused (it never reaches the rows).
#[test]
fn query_letters_past_the_query_end_show_the_gap_character() {
    let query = record(b"AACCGGTT");
    let subject = record(b"ACGTACGTAC");
    let letters = |from: usize, to: usize| subject.seq()[from..to].to_ascii_uppercase();
    let outside = |frame: i32, offset: i32, end: i32, s: (usize, usize), ops| {
        let mut hsp = hsp(frame, (0, 0), s, Some(ops));
        (hsp.internal_q_offset_0, hsp.internal_q_end_0) = (offset, end);
        hsp
    };
    let gapped = || vec![GapEditOp::Sub(2), GapEditOp::Del(1), GapEditOp::Sub(1)];
    for (hsp, start, rows) in [
        (outside(1, 10, 13, (0, 4), gapped()), 10, ("----", "ACGT")),
        (outside(-1, -5, -2, (0, 4), gapped()), 10, ("----", "ACGT")),
        (
            outside(1, 6, 10, (0, 4), vec![GapEditOp::Sub(4)]),
            6,
            ("TT--", "ACGT"),
        ),
        // Context "--AAC" (record positions 9, 8, 7, 6, 5), shown reverse-complemented.
        (
            outside(-1, -2, 3, (0, 5), vec![GapEditOp::Sub(5)]),
            5,
            ("GTT--", "TACGT"),
        ),
    ] {
        assert_eq!(hsp.seqalign_query_start(), start);
        let (q, s) =
            displayed_rows(&hsp, whole(query.seq()), subject.seq().len(), letters).unwrap();
        assert_eq!(
            (q.as_slice(), s.as_slice()),
            (rows.0.as_bytes(), rows.1.as_bytes())
        );
    }
    let left = outside(1, -4, -1, (0, 3), vec![GapEditOp::Sub(3)]);
    assert_eq!(left.seqalign_query_start(), -4);
    assert!(displayed_rows(&left, whole(query.seq()), subject.seq().len(), letters).is_err());
    let right_minus = outside(-1, 9, 12, (0, 3), vec![GapEditOp::Sub(3)]);
    assert_eq!(right_minus.seqalign_query_start(), -4);
    assert!(displayed_rows(
        &right_minus,
        whole(query.seq()),
        subject.seq().len(),
        letters
    )
    .is_err());
}

// remaining-gaps X1b-p review 1 (M1): with `-query_loc` NCBI moves the query row of the
// Seq-align by the interval's start (`RemapToQueryLoc`) and reads the letters of the query's
// Bioseq, the whole record (alnvec.cpp:108-124): an HSP outside its context shows the record's
// letters outside the interval, and the gap character only past the record's end. The record
// GGAACCGGTTAC is searched from position 2 (the interval AACCGGTT, the context of 8 letters):
// past the interval's end on the plus strand and left of the minus context (offsets below 0),
// then a start below 0 that the interval's start moves to the record's first letters (an
// ungapped search's Dense-diag, `BlastnHsp::record_query_start`), on both strands. A start
// that stays below 0 in the record never reaches the rows.
#[test]
fn query_letters_outside_the_interval_are_the_records() {
    let query = record(b"GGAACCGGTTAC");
    let in_record = QueryRecord {
        letters: query.seq(),
        offset: 2,
    };
    let subject = record(b"ACGTACGTAC");
    let letters = |from: usize, to: usize| subject.seq()[from..to].to_ascii_uppercase();
    let outside = |frame: i32, offset: i32, end: i32| {
        let n = (end - offset) as usize;
        let mut hsp = hsp(frame, (0, 0), (0, n), Some(vec![GapEditOp::Sub(n as u32)]));
        (hsp.internal_q_offset_0, hsp.internal_q_end_0) = (offset, end);
        hsp
    };
    for (hsp, start, row) in [
        (outside(1, 6, 10), 8, "TTAC"),
        (outside(1, 8, 11), 10, "AC-"),
        (outside(1, 10, 14), 12, "----"),
        // Context offsets -2..3: record positions 11 down to 7.
        (outside(-1, -2, 3), 7, "GTTAC"),
        (outside(1, -2, 1), 0, "GGA"),
        // Context offsets 7..10: record positions 2 down to 0.
        (outside(-1, 7, 10), 0, "GGA"),
    ] {
        assert_eq!(hsp.record_query_start(in_record.offset, true), start);
        let (q, _) = displayed_rows(&hsp, in_record, subject.seq().len(), letters).unwrap();
        assert_eq!(q.as_slice(), row.as_bytes());
    }
    let left = outside(1, -4, -1);
    assert_eq!(left.record_query_start(2, true), -2);
    assert!(displayed_rows(&left, in_record, subject.seq().len(), letters).is_err());
    // A gapped search's Dense-seg keeps a start below 0 (`CDense_seg::OffsetRow`).
    assert_eq!(outside(1, -2, 1).record_query_start(2, false), -2);
}

// The records the rows read (`QueryRecords::record`): the searched record itself without
// `-query_loc`; with it the input record of the query's input position (records that the
// interval leaves out counted), from the call's first input record, and the interval's start.
#[test]
fn query_records_map_searched_queries_to_their_input_records() {
    use crate::blastinput::seq_range::{Placements, QueryInput, RecordPlacement};
    let input = [record(b"AAAA"), record(b"CC"), record(b"GGGGGG")];
    let searched = [record(b"AA"), record(b"GGG")];
    let ranged = QueryInput {
        placements: Placements(vec![
            RecordPlacement {
                offset: 1,
                length: 4,
            },
            RecordPlacement {
                offset: 1,
                length: 6,
            },
        ]),
        ordinals: vec![0, 2],
        input_lengths: vec![4, 2, 6],
        skipped: vec![false, true, false],
        range_length: Some(2),
    };
    let records = QueryRecords::new(&searched, &ranged, &input, 0);
    let second = records.record(1).unwrap();
    assert_eq!((second.letters, second.offset), (&b"GGGGGG"[..], 1));
    assert!(records.record(2).is_none());
    // A batch's records from input record 2 on.
    let batch = QueryRecords::new(&searched, &ranged, &input[2..], 2);
    assert_eq!(batch.record(1).unwrap().letters, b"GGGGGG");
    assert!(batch.record(0).is_none());
    let plain = QueryInput::whole(&searched);
    let unranged = QueryRecords::new(&searched, &plain, &[], 0)
        .record(1)
        .unwrap();
    assert_eq!((unranged.letters, unranged.offset), (&b"GGG"[..], 0));
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
        let once = displayed_rows(&hsp, whole(query.seq()), subject.seq().len(), letters).unwrap();
        displayed_rows_into(
            &hsp,
            whole(query.seq()),
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
    set_subject_row_masks(
        &mut masks,
        &row,
        10,
        -1,
        &shown,
        crate::report::query_anchored::NO_ROW,
    );
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
