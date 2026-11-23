use super::*;
use crate::ncbi_exception::NCBI_BUILD_ROOT;

// NCBI BLAST+ 2.17.0 (the official Linux x64 release), rows db.x2g.bx.f0616.0,
// db.x2g.bx.f0616.subj.6 and db.rpstbn.sw.nest.0: `CATCH_ALL`'s report of the three exceptions.
#[test]
fn the_exceptions_are_ncbis_reports() {
    // NCBI's reports in LOSAT's words (`crate::wording::losat_text`).
    assert_eq!(
        crate::cli::NativeError::exception(&row_range_exception(RowRangeError::Empty)).message,
        crate::wording::losat_text(&format!(
            "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 227: Error: (CSeqalignException::eInvalidAlignment) CSeq_align::GetSeqRange(): Row is empty\n\n"
        ))
        .as_bytes()
    );
    assert_eq!(
        crate::cli::NativeError::exception(&row_range_exception(
            RowRangeError::SeqidsNotConsistent
        ))
        .message,
        crate::wording::losat_text(&format!(
            "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 186: Error: (CSeqalignException::eInvalidRowNumber) CSeq_align::GetSeqRange(): Row seqids not consistent. Cannot determine range.\n\n"
        ))
        .as_bytes()
    );
    assert_eq!(
        crate::cli::NativeError::exception(&std_seg_lengths_not_accurate()).message,
        crate::wording::losat_text(&format!(
            "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 801: Error: (CSeqalignException::eInvalidInputAlignment) CreateDensegFromStdseg(): Std-seg segment lengths not accurate!\n\n"
        ))
        .as_bytes()
    );
}

fn ungapped(q_start: i32, q_end: i32) -> super::super::results::Hsp {
    super::super::results::Hsp {
        hsp: PreliminaryHsp {
            context: 0,
            frame: 1,
            score: 5,
            q_start,
            q_end,
            q_gapped_start: q_start,
            s_start: 1,
            s_end: 1 + q_end - q_start,
            s_gapped_start: 1,
        },
        num: 1,
        evalue: 1.4,
        bit_score: 5.2,
        identity: 0,
        positive: 0,
        edit_script: Vec::new(),
        composition_method: 0,
    }
}

// The ws gate's S1 rows (db.x2g.bx.s1.*): a subject whose every HSP has length 0 has a
// Seq-align of no row of letters (NCBI's prefetch writes its warning, `-subject`'s mapper
// finds no sequence type); one HSP of another length gives it a row of letters.
#[test]
fn a_subject_of_empty_rows_only_needs_every_hsp_of_length_0() {
    let list = |hsps| super::super::results::HspList {
        query_index: 0,
        oid: 0,
        hsps,
        best_evalue: 1.4,
    };
    assert!(subject_rows_all_empty(&list(vec![ungapped(1, 1)]), false));
    assert!(subject_rows_all_empty(
        &list(vec![ungapped(1, 1), ungapped(4, 4)]),
        false
    ));
    assert!(!subject_rows_all_empty(
        &list(vec![ungapped(1, 1), ungapped(2, 5)]),
        false
    ));
    assert!(!subject_rows_all_empty(&list(Vec::new()), false));
    // A gapped search has no ungapped HSP of length 0.
    assert!(!subject_rows_all_empty(&list(vec![ungapped(1, 1)]), true));
}

// X/SURVEY_X2g_formats.md 3 (set C): -subject throws only with a field of the alignment
// block; -db leaves at the query's fields only for a subject's only HSP with a query field.
// A -subject subject's only HSP of length 0 (decisions X2g-5 and X2g-7): the mapper's
// alignment of letters, no exception.
#[test]
fn the_tabular_row_of_empty_rows_follows_the_mode_and_the_fields() {
    use SubjectLoc::{Interval, Whole};
    use TabularEmptyRow::*;
    assert_eq!(tabular_empty_row(Interval, true, 2, true, false), Throws);
    assert_eq!(tabular_empty_row(Interval, true, 3, false, false), Throws);
    assert_eq!(tabular_empty_row(Interval, false, 2, true, false), Row);
    assert_eq!(tabular_empty_row(Interval, true, 1, true, false), Mapped);
    assert_eq!(tabular_empty_row(Interval, false, 1, false, false), Mapped);
    // The mapped query row wraps (fix3_oob.py's minus-strand pairs): as `-db`'s empty row
    // with a query field, else the Dense-seg's 801 with a field of the block.
    assert_eq!(tabular_empty_row(Interval, true, 1, true, true), Leaves);
    assert_eq!(tabular_empty_row(Interval, true, 1, false, true), Throws);
    assert_eq!(tabular_empty_row(Interval, false, 1, false, true), Row);
    assert_eq!(tabular_empty_row(Whole, true, 1, true, false), Leaves);
    assert_eq!(tabular_empty_row(Whole, false, 1, true, false), Leaves);
    assert_eq!(tabular_empty_row(Whole, true, 2, true, false), Row);
    assert_eq!(tabular_empty_row(Whole, true, 1, false, false), Row);
}

fn record(sequence: &str) -> super::super::input::FastaRecord {
    super::super::input::FastaRecord {
        internal_id: "Query_1".to_string(),
        title: Vec::new(),
        sequence: sequence.as_bytes().to_vec(),
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    }
}

fn mapped(query: &super::super::input::FastaRecord, subject: &str, frame: i8) -> MappedRows {
    let mut hsp = ungapped(1, 1).hsp;
    hsp.frame = frame;
    MappedRows::of(&hsp, query, &record(subject)).expect("an HSP inside its subject")
}

// Decision X2g-7, NCBI BLAST+ 2.17.0 on ws (`X/x2g/fix3`, fix3_grid.py): the rows of the
// mapped alignment of a subject's only HSP of length 0 at query offset 1 and subject offset 1.
#[test]
fn the_mapped_rows_are_ncbis() {
    // x2g_s1 (c01241, frame 1): `c01241_q c01241_s 4 7 2 5 4 TCCC TLWI ... 1`.
    let MappedRows::Letters(plus) = mapped(&record("TTCTCCCTGCCCATA"), "GTLWI", 1) else {
        panic!("letters");
    };
    assert_eq!(
        (plus.q_start, plus.q_end, plus.s_start, plus.s_end),
        (4, 7, 2, 5)
    );
    assert_eq!(plus.rows(), (&b"TCCC"[..], &b"TLWI"[..]));
    assert_eq!((plus.q_frame, plus.length()), (1, 4));
    // Its minus-frame twin (an HSP of frame -1, whose empty row starts at 12):
    // `11 8 2 5 4 CCCT TLWI ... -2`, the row two below that start.
    let MappedRows::Letters(minus) = mapped(&record("TATGGGCAGGGAGAA"), "GTLWI", -1) else {
        panic!("letters");
    };
    assert_eq!((minus.q_start, minus.q_end, minus.q_frame), (11, 8, -2));
    assert_eq!(minus.rows(), (&b"CCCT"[..], &b"TLWI"[..]));
    // x2g_s1o (o00016): one nucleotide past the query's end, the gap character.
    let MappedRows::Letters(past) = mapped(&record("TGCTCGGCATGGCTG"), "ETANLXXXXXXXXX", 1) else {
        panic!("letters");
    };
    assert_eq!((past.q_start, past.q_end, past.length()), (4, 16, 13));
    assert_eq!(past.rows().0, b"TCGGCATGGCTG-");
    // x2g_s1w: the minus twin would start before the query's first nucleotide.
    assert!(matches!(
        mapped(&record("CAGCCATGCCGAGCA"), "ETANLXXXXXXXXX", -1),
        MappedRows::QueryRowWraps
    ));
    // x2g_s1q: -query_loc 1-15 of 43 nucleotides; the row's 16th is the record's, which the
    // cut keeps around the interval (`query_letters_around`: the subject's 14 letters and one).
    let mut cut = record("TGCTCGGCATGGCTGACGTACGTACGTACGTACGTACGTACGT");
    cut.keep_around(0, 15, 15);
    cut.cut_in_place(0, 15, false);
    let MappedRows::Letters(read) = mapped(&cut, "ETANLXXXXXXXXX", 1) else {
        panic!("letters");
    };
    assert_eq!((read.q_start, read.q_end), (4, 16));
    assert_eq!(read.rows().0, b"TCGGCATGGCTGA");
    // Without the kept letters the record does not hold the 16th: no rows (refused before).
    let mut bare = record("TGCTCGGCATGGCTGACGTACGTACGTACGTACGTACGTACGT");
    bare.cut_in_place(0, 15, false);
    assert!(MappedRows::of(&ungapped(1, 1).hsp, &bare, &record("ETANLXXXXXXXXX")).is_none());
    // An interval starting past the record's first letter (review 4 S2): `-query_loc 5-11` of
    // 12 nucleotides, a minus-frame empty row at 8 of the record; the row ends two below it
    // and starts one before the interval (the record's 4th letter, kept around it).
    let mut minus_cut = record("TTTAGGGTACGT");
    minus_cut.keep_around(4, 11, 6);
    minus_cut.cut_in_place(4, 11, false);
    let MappedRows::Letters(shifted) = mapped(&minus_cut, "GTLWI", -1) else {
        panic!("letters");
    };
    assert_eq!((shifted.q_start, shifted.q_end), (7, 4));
    assert_eq!(shifted.rows().0, b"CCCT");
    assert_eq!(
        (0..12)
            .map(|at| minus_cut.record_letter(at))
            .collect::<Option<Vec<u8>>>(),
        Some(b"TTTAGGGTACGT".to_vec())
    );
    assert_eq!(minus_cut.record_letter(12), None);
    // The farthest a minus row reads before the interval (review 5 N3): an HSP at the end of
    // frame -1 of `-query_loc 7-9` (its empty row starts at the interval's first nucleotide)
    // and at the first letter of a subject of 5 letters reads from 6 (5 + 1) letters before
    // the interval: `query_letters_around`'s one more than the subject's letters.
    let farthest = |around: usize| {
        let mut cut = record("ACGTACGGGTT");
        cut.keep_around(6, 9, around);
        cut.cut_in_place(6, 9, false);
        let mut hsp = ungapped(1, 1).hsp;
        (hsp.frame, hsp.s_start, hsp.s_end) = (-1, 0, 0);
        MappedRows::of(&hsp, &cut, &record("GTLWI"))
    };
    let Some(MappedRows::Letters(far)) = farthest(6) else {
        panic!("letters");
    };
    assert_eq!((far.q_start, far.q_end, far.length()), (5, 1, 5));
    assert_eq!(far.rows().0, b"TACGT");
    assert!(farthest(5).is_none());
}

// Review 4 M1 (ws, NCBI BLAST+ 2.17.0, `X/x2g/rev4`): `x_ConvertSegment`'s stop wraps, and the
// mapper rewrites an empty row, only at the first or second letter of a subject without a row
// of letters, and only at the first of a protein subject (X2g-4).
#[test]
fn the_mapper_rewrites_an_empty_row_only_where_its_stop_wraps() {
    assert!(mapper_wraps_empty_row(0, false));
    assert!(mapper_wraps_empty_row(1, false));
    assert!(!mapper_wraps_empty_row(2, false));
    assert!(!mapper_wraps_empty_row(7, false));
    assert!(mapper_wraps_empty_row(0, true));
    assert!(!mapper_wraps_empty_row(1, true));
    let at = |s_start: i32| {
        let mut hsp = ungapped(1, 1);
        hsp.hsp.s_start = s_start;
        hsp.hsp.s_end = s_start;
        hsp
    };
    let list = |hsps| super::super::results::HspList {
        query_index: 0,
        oid: 0,
        hsps,
        best_evalue: 1.4,
    };
    use SubjectLoc::{Interval, Whole};
    assert_eq!(Interval.of_list(&list(vec![at(1)]), false), Interval);
    assert_eq!(Interval.of_list(&list(vec![at(3)]), false), Whole);
    assert_eq!(Interval.of_list(&list(vec![at(3), at(5)]), false), Whole);
    assert_eq!(Interval.of_list(&list(vec![at(3), at(0)]), false), Interval);
    // A subject with letters keeps the run's Seq-loc.
    assert_eq!(
        Interval.of_list(&list(vec![at(3), ungapped(2, 5)]), false),
        Interval
    );
    assert_eq!(Whole.of_list(&list(vec![at(1)]), false), Whole);
}

// X/SURVEY_X2g_formats.md 1: p0095 (frame -1, raw query row (119, 118)) maps to 118..118;
// a plus frame keeps its start.
#[test]
fn the_mapped_nucleotide_follows_the_strand() {
    use crate::report::serial::seqalign::mapped_empty_query_nucleotide;
    assert_eq!(mapped_empty_query_nucleotide(119, -1), 118);
    assert_eq!(mapped_empty_query_nucleotide(81, 1), 81);
    assert_eq!(mapped_empty_query_nucleotide(0, 1), 0);
}

// Line PS item PS12d (`PR/ps12d/scripts/gcmp.sh`, NCBI BLAST+ 2.17.0): RPSTBLASTN's
// Smith-Waterman alignment of nest_p.fna that starts with 95 gaps in the profile's row (`Ins`
// first: the query's letters against the gap). Its Std-seg has no location in that row's first
// segment: `GetSeqRange` of the subject row throws `Row seqids not consistent`, the query row's is
// taken; a gap at the end or in the middle leaves both. A `Del` first is the query row's.
#[test]
fn a_row_that_starts_with_a_gap_has_no_range() {
    use crate::common::GapEditOp::{Del, Ins, Sub};
    let gapped = |ops| super::super::results::Hsp {
        edit_script: ops,
        ..ungapped(0, 285)
    };
    let starts = gapped(vec![Ins(95), Sub(190)]);
    assert_eq!(row_range(&starts, AlignRow::Query), None);
    assert_eq!(
        row_range(&starts, AlignRow::Subject),
        Some(RowRangeError::SeqidsNotConsistent)
    );
    assert_eq!(
        rows_range(&starts),
        Some(RowRangeError::SeqidsNotConsistent)
    );
    let query_first = gapped(vec![Del(3), Sub(190)]);
    assert_eq!(
        row_range(&query_first, AlignRow::Query),
        Some(RowRangeError::SeqidsNotConsistent)
    );
    assert_eq!(row_range(&query_first, AlignRow::Subject), None);
    for ops in [
        vec![Sub(190), Ins(95)],
        vec![Sub(90), Del(5), Sub(100)],
        Vec::new(),
    ] {
        assert_eq!(rows_range(&gapped(ops.clone())), None, "{ops:?}");
    }
    // An HSP of length 0: both rows empty, the query's first.
    assert_eq!(rows_range(&ungapped(4, 4)), Some(RowRangeError::Empty));
    assert_eq!(
        row_range(&ungapped(4, 4), AlignRow::Subject),
        Some(RowRangeError::Empty)
    );
}
