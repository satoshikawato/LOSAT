//! The ASN.1 text (outfmt 8), binary (outfmt 9) and JSON (outfmt 12) documents of `write_query` against NCBI
//! BLAST+ 2.17.0's own output for the same HSPs, byte for byte.
//!
//! The expected files in `tests/fixtures/seqalign/` are NCBI's output, unchanged
//! (outfmt-20261019 `ncbi-samples/out/` and the OF3 design probes `S/tmp/probes/`):
//! - `blastn-megablast_subject_fmt{8,9,12}.out`: `blastn -task megablast -query fixture/qnucl.fa
//!   -subject fixture/snucl.fa -num_threads 1 -outfmt 8|9|12` (one plus-strand and one
//!   minus-strand HSP with gaps on both rows; e-value below 1e-180 written as it is here,
//!   the second query has no hits).
//! - `blastp_db_fmt{8,9,12}.out`: `blastp -query fixture/qprot.fa -db db/sprot -num_threads 1
//!   -outfmt 8|9|12` (two documents, BL_ORD_ID subjects, composition adjustment, positives).
//! - `blastp_subject_fmt{8,9,12}.out`: the same with `-subject fixture/sprot.fa` (local subject
//!   IDs, the strands of the subject remap).
//! - `blastp_dblong_fmt{8,9,12}.out`, `blastp_dbnosp_fmt{8,9,12}.out`: the same search against the
//!   same sequences in databases made with a long `-title` (with spaces and quotes: the
//!   ASN.1 text wraps at a space; without spaces: a hard break at column 78).
//!
//! The records are those of NCBI's HSPs (the starts, lengths and scores printed in the
//! files; an e-value or bit score as its 15 printed digits, which give back the same double).

use super::encoder::{Encoder, SerialError};
use super::{write_query, SeqAlignError, SerialFormat};
use crate::blastinput::app::Program;
use crate::common::GapEditOp::{self, Del, Ins, Sub};
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;

macro_rules! fixture {
    ($name:literal) => {
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/seqalign/",
            $name
        ))
    };
}

/// An HSP record of the tests (frames: query, subject).
#[allow(clippy::too_many_arguments)]
fn hsp<'a>(
    subject: u32,
    (q_offset, s_offset): (i32, i32),
    frames: (i8, i8),
    score: i32,
    evalue: f64,
    bit_score: f64,
    (num_ident, num_positives, comp_adjust): (i32, i32, u8),
    ops: &'a [GapEditOp],
) -> HspRecord<'a> {
    let q_len: u32 = ops
        .iter()
        .map(|op| match op {
            Sub(n) | Ins(n) => *n,
            Del(_) => 0,
        })
        .sum();
    let s_len: u32 = ops
        .iter()
        .map(|op| match op {
            Sub(n) | Del(n) => *n,
            Ins(_) => 0,
        })
        .sum();
    HspRecord {
        subject,
        q_offset,
        q_end: q_offset + q_len as i32,
        s_offset,
        s_end: s_offset + s_len as i32,
        q_frame: frames.0,
        s_frame: frames.1,
        score,
        num: 0,
        evalue,
        bit_score,
        num_ident,
        num_positives,
        comp_adjust,
        ops,
    }
}

fn query(id: &str, searched_length: i32) -> QueryRecord<'_> {
    QueryRecord {
        id: SeqId::Local(id),
        searched_length,
        offset: 0,
    }
}

fn subject(id: SeqId<'_>, searched_length: i32) -> SubjectRecord<'_> {
    SubjectRecord {
        id,
        searched_length,
        offset: 0,
        use_this_seqid: &[],
    }
}

const MB_HSP1: [GapEditOp; 5] = [Sub(200), Ins(4), Sub(100), Del(2), Sub(119)];
const MB_HSP2: [GapEditOp; 9] = [
    Sub(164),
    Ins(4),
    Sub(1),
    Ins(2),
    Sub(94),
    Del(1),
    Sub(4),
    Del(2),
    Sub(146),
];

/// `blastn -task megablast -subject`: Query_1 (480 nt) against Subject_1 (plus strand) and
/// Subject_2 (minus strand); the minus-strand HSP starts at offset 30 of the reverse
/// complement (its first query start is 480 - 30 - 164 = 286).
fn megablast_documents(format: SerialFormat) -> Vec<u8> {
    let run = RunRecord {
        program: Program::Blastn,
        db_title: None,
        gapped: true,
    };
    let subjects = [
        subject(SeqId::Local("Subject_1"), 500),
        subject(SeqId::Local("Subject_2"), 500),
    ];
    let hsps = [
        hsp(
            0,
            (30, 45),
            (1, 1),
            329,
            2.60916226603389e-178,
            608.668029852878,
            (394, 0, 0),
            &MB_HSP1,
        ),
        hsp(
            1,
            (30, 40),
            (-1, 1),
            242,
            6.01914206611342e-130,
            448.009510099483,
            (361, 0, 0),
            &MB_HSP2,
        ),
    ];
    let mut buf = Vec::new();
    assert_eq!(
        write_query(
            format,
            &run,
            &query("Query_1", 480),
            &subjects,
            &hsps,
            &mut buf
        ),
        Ok(true)
    );
    // Query_2 has no hits: no document.
    assert_eq!(
        write_query(
            format,
            &run,
            &query("Query_2", 90),
            &subjects,
            &[],
            &mut buf
        ),
        Ok(false)
    );
    buf
}

const BP_Q1_HSP1: [GapEditOp; 3] = [Sub(60), Ins(3), Sub(77)];
const BP_Q1_HSP2: [GapEditOp; 5] = [Sub(40), Del(4), Sub(55), Ins(2), Sub(42)];
const BP_Q2_HSP1: [GapEditOp; 1] = [Sub(6)];
const BP_Q2_HSP2: [GapEditOp; 1] = [Sub(9)];

/// `blastp -db`: Query_1 (140 aa) against OIDs 0 and 1, Query_2 (45 aa) against OID 0 twice.
fn blastp_db_documents(format: SerialFormat, db_title: &[u8]) -> Vec<u8> {
    blastp_documents(format, Some(db_title))
}

/// `blastp -db` (`db_title`) or `-subject` (`None`, Subject_1 and Subject_2: the same HSPs,
/// with the strands of the subject remap).
fn blastp_documents(format: SerialFormat, db_title: Option<&[u8]>) -> Vec<u8> {
    let run = RunRecord {
        program: Program::Blastp,
        db_title,
        gapped: true,
    };
    let subjects = match db_title {
        Some(_) => [
            subject(SeqId::BlOrdId(0), 137),
            subject(SeqId::BlOrdId(1), 150),
        ],
        None => [
            subject(SeqId::Local("Subject_1"), 137),
            subject(SeqId::Local("Subject_2"), 150),
        ],
    };
    let query1 = [
        hsp(
            0,
            (0, 0),
            (0, 0),
            652,
            5.30325651018409e-94,
            255.758355778158,
            (125, 130, 2),
            &BP_Q1_HSP1,
        ),
        hsp(
            1,
            (1, 1),
            (0, 0),
            477,
            2.8220561436493e-67,
            188.348429992622,
            (100, 110, 2),
            &BP_Q1_HSP2,
        ),
    ];
    let query2 = [
        hsp(
            0,
            (32, 34),
            (0, 0),
            21,
            1.51524398110052,
            12.6974233743084,
            (4, 4, 2),
            &BP_Q2_HSP1,
        ),
        hsp(
            0,
            (23, 0),
            (0, 0),
            20,
            1.88521902451112,
            12.3122237983911,
            (3, 7, 2),
            &BP_Q2_HSP2,
        ),
    ];
    let mut buf = Vec::new();
    for (query, hsps) in [
        (query("Query_1", 140), &query1[..]),
        (query("Query_2", 45), &query2[..]),
    ] {
        assert_eq!(
            write_query(format, &run, &query, &subjects, hsps, &mut buf),
            Ok(true)
        );
    }
    buf
}

#[track_caller]
fn assert_same_bytes(actual: &[u8], expected: &[u8]) {
    if actual != expected {
        let first = actual
            .iter()
            .zip(expected)
            .position(|(a, b)| a != b)
            .unwrap_or(actual.len().min(expected.len()));
        panic!(
            "documents differ at byte {first} (lengths {} and {}):\n--- written\n{}\n--- NCBI\n{}",
            actual.len(),
            expected.len(),
            String::from_utf8_lossy(
                &actual[first.saturating_sub(200)..actual.len().min(first + 200)]
            ),
            String::from_utf8_lossy(
                &expected[first.saturating_sub(200)..expected.len().min(first + 200)]
            ),
        );
    }
}

#[test]
fn blastn_megablast_subject_documents_equal_ncbi() {
    assert_same_bytes(
        &megablast_documents(SerialFormat::AsnText),
        fixture!("blastn-megablast_subject_fmt8.out"),
    );
    assert_same_bytes(
        &megablast_documents(SerialFormat::AsnBinary),
        fixture!("blastn-megablast_subject_fmt9.out"),
    );
    assert_same_bytes(
        &megablast_documents(SerialFormat::Json),
        fixture!("blastn-megablast_subject_fmt12.out"),
    );
}

#[test]
fn blastp_db_documents_equal_ncbi() {
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::AsnText, b"fixture/sprot.fa"),
        fixture!("blastp_db_fmt8.out"),
    );
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::AsnBinary, b"fixture/sprot.fa"),
        fixture!("blastp_db_fmt9.out"),
    );
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::Json, b"fixture/sprot.fa"),
        fixture!("blastp_db_fmt12.out"),
    );
}

/// A `-subject` search's Seq-aligns are remapped to the subjects (`s_RemapToSubjectLoc`):
/// the protein rows' `unknown` strands become `plus`.
#[test]
fn blastp_subject_documents_equal_ncbi() {
    assert_same_bytes(
        &blastp_documents(SerialFormat::AsnText, None),
        fixture!("blastp_subject_fmt8.out"),
    );
    assert_same_bytes(
        &blastp_documents(SerialFormat::AsnBinary, None),
        fixture!("blastp_subject_fmt9.out"),
    );
    assert_same_bytes(
        &blastp_documents(SerialFormat::Json, None),
        fixture!("blastp_subject_fmt12.out"),
    );
}

/// A title past column 78 wraps before its last space at a column up to 78; the doubled
/// quotes of `"quoted"` stay together; JSON escapes the quotes and does not wrap.
#[test]
fn long_database_title_wraps_at_a_space() {
    let title: &[u8] = b"a long database title with spaces so that the label string wraps past column seventy eight \"quoted\" end";
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::AsnText, title),
        fixture!("blastp_dblong_fmt8.out"),
    );
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::AsnBinary, title),
        fixture!("blastp_dblong_fmt9.out"),
    );
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::Json, title),
        fixture!("blastp_dblong_fmt12.out"),
    );
}

/// A title without spaces breaks at column 78.
#[test]
fn long_database_title_without_spaces_breaks_at_column_78() {
    let title: &[u8] = b"fixture/a_very_long_path_without_any_space_characters_to_force_a_hard_break_at_column_78_of_the_line.fa";
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::AsnText, title),
        fixture!("blastp_dbnosp_fmt8.out"),
    );
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::AsnBinary, title),
        fixture!("blastp_dbnosp_fmt9.out"),
    );
    assert_same_bytes(
        &blastp_db_documents(SerialFormat::Json, title),
        fixture!("blastp_dbnosp_fmt12.out"),
    );
}

/// The score list of a linked HSP (`num > 1`: `sum_n`, then the e-value as `sum_e`), an
/// e-value below `SMALLEST_EVALUE` written as 0, and the conditions of the optional scores
/// (blast_seqalign.cpp:1162-1231).
#[test]
fn score_list_follows_ncbi_conditions() {
    let run = RunRecord {
        program: Program::Blastp,
        db_title: None,
        gapped: true,
    };
    let subjects = [subject(SeqId::Local("Subject_1"), 100)];
    let ops = [Sub(10)];
    let mut linked = hsp(0, (0, 0), (0, 0), 30, 1e-200, 20.5, (0, 0, 0), &ops);
    linked.num = 2;
    let mut buf = Vec::new();
    write_query(
        SerialFormat::AsnText,
        &run,
        &query("Query_1", 10),
        &subjects,
        &[linked],
        &mut buf,
    )
    .unwrap();
    let text = String::from_utf8(buf).unwrap();
    let ids: Vec<&str> = text
        .lines()
        .filter_map(|line| line.trim().strip_prefix("id str \""))
        .map(|rest| rest.trim_end_matches("\","))
        .collect();
    assert_eq!(
        ids,
        [
            "score",
            "blast_score",
            "sum_n",
            "sum_e",
            "bit_score",
            "num_ident",
            "hsp_percent_coverage"
        ]
    );
    assert!(text.contains("value real { 0, 10, 0 }"), "{text}");
    // 100 % coverage is not raised by 0.5; num_ident 0 is written, num_positives 0 is not.
    assert!(text.contains("value real { 1, 10, 2 }"), "{text}");
    assert!(text.contains("value int 0"), "{text}");
}

/// HSPs of score 0 have no Seq-align (blast_seqalign.cpp:672-674); a query with only such
/// HSPs writes nothing.
#[test]
fn score_zero_hsps_write_nothing() {
    let run = RunRecord {
        program: Program::Blastn,
        db_title: None,
        gapped: true,
    };
    let subjects = [subject(SeqId::Local("Subject_1"), 100)];
    let ops = [Sub(10)];
    let zero = hsp(0, (0, 0), (1, 1), 0, 1.0, 1.0, (10, 0, 0), &ops);
    let mut buf = b"kept".to_vec();
    for format in [
        SerialFormat::AsnText,
        SerialFormat::AsnBinary,
        SerialFormat::Json,
    ] {
        assert_eq!(
            write_query(
                format,
                &run,
                &query("Query_1", 10),
                &subjects,
                &[zero],
                &mut buf
            ),
            Ok(false)
        );
    }
    assert_eq!(buf, b"kept");
}

/// The ungapped searches whose Seq-aligns would be Dense-diags (BLASTN, BLASTP; neither runs
/// an ungapped search) are refused, not written.
#[test]
fn untranslated_ungapped_runs_are_refused() {
    let subjects = [subject(SeqId::Local("Subject_1"), 100)];
    let ops = [Sub(10)];
    let one = hsp(0, (0, 0), (1, 0), 30, 1.0, 1.0, (10, 0, 0), &ops);
    for (program, gapped) in [(Program::Blastn, false), (Program::Blastp, false)] {
        let run = RunRecord {
            program,
            db_title: None,
            gapped,
        };
        let mut buf = Vec::new();
        assert_eq!(
            write_query(
                SerialFormat::AsnText,
                &run,
                &query("Query_1", 10),
                &subjects,
                &[one],
                &mut buf
            ),
            Err(SeqAlignError::Unsupported(program))
        );
        assert!(buf.is_empty());
    }
}

/// `WriteEscapedChar`: quotes and backslashes escaped, bytes below 0x20 as `\u00XX`
/// (upper-case hex), other bytes copied; `x_WriteString` stops at a NUL byte.
#[test]
fn json_strings_are_escaped_as_ncbi() {
    let mut buf = Vec::new();
    let mut enc = super::json::JsonEncoder::new(&mut buf);
    enc.visible_str(b"a\"b\\c\x01\x1f\x7f\xc3\xa9d\0tail");
    assert_eq!(enc.error(), None);
    assert_eq!(buf, b"\"a\\\"b\\\\c\\u0001\\u001F\x7f\xc3\xa9d\"");
}

/// A byte outside the visible range is an error of the ASN.1 text encoder (the callers
/// refuse such titles before the search).
#[test]
fn asn_text_records_a_non_visible_byte() {
    let mut buf = Vec::new();
    let mut enc = super::asn_text::AsnTextEncoder::new(&mut buf);
    enc.visible_str(b"a\tb");
    assert_eq!(enc.error(), Some(SerialError::NonVisibleChar(b'\t')));
}

/// `COStreamBuffer::PutInt4`.
#[test]
fn integers_are_decimal() {
    for (value, text) in [
        (0, "0"),
        (-1, "-1"),
        (42, "42"),
        (i32::MAX, "2147483647"),
        (-2147483647, "-2147483647"),
    ] {
        let mut digits = super::asn_text::itoa_i32(value);
        assert_eq!(digits.text(), text.as_bytes());
    }
}

/// LOSAT's refusals of a database for the Seq-annot formats (`check_database`): a database
/// named `n/a` or `SRA` (any case), and for the ASN.1 text and binary a title byte outside 0x20..0x7E.
#[test]
fn databases_that_ncbi_describes_otherwise_are_refused() {
    use super::seqalign::check_database;
    for name in ["n/a", "N/A", "sra", "SRA"] {
        for format in [
            SerialFormat::AsnText,
            SerialFormat::AsnBinary,
            SerialFormat::Json,
        ] {
            let error = check_database(format, name, b"title", Program::Blastn).unwrap_err();
            assert!(
                error
                    .to_string()
                    .ends_with("is not supported by LOSAT's BLASTN"),
                "{error}"
            );
        }
    }
    assert!(check_database(SerialFormat::AsnText, "db/x", b"a\tb", Program::Blastp).is_err());
    assert!(check_database(SerialFormat::AsnText, "db/x", b"\xc3\xa9", Program::Blastp).is_err());
    assert!(check_database(SerialFormat::AsnBinary, "db/x", b"a\tb", Program::Blastp).is_err());
    assert!(check_database(
        SerialFormat::AsnBinary,
        "db/x",
        b"\xc3\xa9",
        Program::Blastp
    )
    .is_err());
    assert!(check_database(
        SerialFormat::AsnBinary,
        "db/x",
        b"plain title",
        Program::Blastp
    )
    .is_ok());
    for title in [&b"a\tb"[..], b"\xc3\xa9", b"plain title"] {
        assert!(check_database(SerialFormat::Json, "db/x", title, Program::Blastp).is_ok());
    }
    assert!(check_database(
        SerialFormat::AsnText,
        "db/nash",
        b"fixture/sprot.fa",
        Program::Blastp
    )
    .is_ok());
}

/// The mapper of a `-subject` search's remap leaves out a first and a last Std-seg with a gap
/// row (seq_align_mapper_base.cpp:1568-1580, `std_seg_kept`); a database search keeps every
/// segment. (BLAST's edit scripts start and end with a substitution; this script does not.)
#[test]
fn subject_remap_drops_std_segs_that_start_or_end_with_a_gap() {
    let ops = [Del(2), Sub(10), Ins(1), Sub(5), Ins(3)];
    let one = hsp(0, (0, 0), (1, 0), 30, 1.0, 1.0, (10, 12, 0), &ops);
    let segments = |db_title: Option<&'static [u8]>| {
        let run = RunRecord {
            program: Program::Blastx,
            db_title,
            gapped: true,
        };
        let subjects = [subject(SeqId::Local("Subject_1"), 100)];
        let mut buf = Vec::new();
        assert_eq!(
            write_query(
                SerialFormat::AsnText,
                &run,
                &query("Query_1", 300),
                &subjects,
                &[one],
                &mut buf
            ),
            Ok(true)
        );
        let text = String::from_utf8(buf).unwrap();
        (
            text.matches("ids {").count(),
            text.matches("empty local").count(),
        )
    };
    assert_eq!(segments(None), (3, 1));
    assert_eq!(segments(Some(b"db")), (5, 3));
}

/// An ungapped translated search writes one `diags` Seq-align per subject (no `dim`, no
/// `score`), one Std-seg with its scores per HSP, the minus frame's interval counted from the
/// end of the searched letters and moved by the ranges' starts
/// (blast_seqalign.cpp:1348-1361, 1406-1447).
#[test]
fn ungapped_translated_hsps_are_diags_seq_aligns() {
    let run = RunRecord {
        program: Program::Tblastx,
        db_title: None,
        gapped: false,
    };
    let subjects = [
        SubjectRecord {
            offset: 5,
            ..subject(SeqId::Local("Subject_1"), 60)
        },
        subject(SeqId::Local("Subject_2"), 90),
    ];
    let query = QueryRecord {
        offset: 20,
        ..query("Query_1", 100)
    };
    let mut first = hsp(0, (0, 2), (-1, 2), 30, 1.0, 20.0, (8, 9, 0), &[]);
    first.q_end = 10;
    first.s_end = 12;
    first.num = 2;
    let mut second = first;
    second.q_frame = 3;
    let mut third = first;
    third.subject = 1;
    third.num = 0;
    let mut buf = Vec::new();
    assert_eq!(
        write_query(
            SerialFormat::AsnText,
            &run,
            &query,
            &subjects,
            &[first, second, third],
            &mut buf
        ),
        Ok(true)
    );
    let text = String::from_utf8(buf).unwrap();
    assert_eq!(text.matches("type diags,").count(), 2);
    assert_eq!(text.matches("scores {").count(), 3);
    assert_eq!(text.matches("sum_n").count(), 2);
    // The Seq-align level has neither `dim` nor `score` (6 spaces in).
    assert!(!text.contains("\n      dim ") && !text.contains("\n      score {"));
    // Query frame -1, offsets 0..10 of 100 letters: 70..99, plus the -query_loc start 20.
    assert!(
        text.contains("from 90,\n              to 119,\n              strand minus,"),
        "{text}"
    );
    // Query frame +3: 3*0+3-1 = 2 .. 3*10+3-2 = 31, plus 20.
    assert!(text.contains("from 22,\n              to 51,\n              strand plus,"));
    // Subject frame +2: 3*2+2-1 = 7 .. 3*12+2-2 = 36, plus the -subject_loc start 5.
    assert!(text.contains("from 12,\n              to 41,\n              strand plus,"));
    // An empty edit script and no gap: no `empty` row.
    assert!(!text.contains("empty"));
}

/// Under a positive list every Seq-align of a subject with strings ends with the `ext`
/// `use_this_seqid`: a gapped one per HSP, a `diags` one per subject; a subject without strings
/// has none (blast_seqalign.cpp:1234-1245; probes E01, E12 of line DI-8).
#[test]
fn a_subject_with_strings_writes_the_use_this_seqid_ext() {
    const EXT: &str = "ext{{typestr\"use_this_seqid\",data{{labelstr\"SEQIDS\",num2,\
         datastrs{\"seqid:P69905.2\",\"seqid:P01942.2\"}}}}}";
    let strings = [b"seqid:P69905.2".to_vec(), b"seqid:P01942.2".to_vec()];
    let ops = [Sub(30)];
    for gapped in [true, false] {
        let run = RunRecord {
            program: if gapped {
                Program::Blastp
            } else {
                Program::Tblastx
            },
            db_title: Some(b"db"),
            gapped,
        };
        let subjects = [
            SubjectRecord {
                use_this_seqid: &strings,
                ..subject(SeqId::Local("Subject_1"), 100)
            },
            subject(SeqId::Local("Subject_2"), 100),
        ];
        let document = |subject_index| {
            let frames = if gapped { (0, 0) } else { (1, 1) };
            let hsps = [hsp(
                subject_index,
                (0, 0),
                frames,
                100,
                1e-5,
                50.0,
                (30, 30, 0),
                &ops,
            )];
            let mut buf = Vec::new();
            assert!(write_query(
                SerialFormat::AsnText,
                &run,
                &query("Query_1", 100),
                &subjects,
                &hsps,
                &mut buf,
            )
            .unwrap());
            String::from_utf8(buf).unwrap().replace([' ', '\n'], "")
        };
        let named = document(0);
        assert!(named.contains(EXT), "{named}");
        assert!(!document(1).contains("use_this_seqid"));
    }
}
