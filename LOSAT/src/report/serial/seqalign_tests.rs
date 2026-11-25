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

use super::encoder::{Encoder, ReplacedChar, SerialError, StreamPosition, Traced};
use super::{write_query, SerialFormat};
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
        q_left_out: 0,
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
        program_type: None,
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
        program_type: None,
    };
    let subjects = match db_title {
        Some(_) => [
            subject(SeqId::BlOrdId(crate::seqdb::Oid::of_unrestricted(0)), 137),
            subject(SeqId::BlOrdId(crate::seqdb::Oid::of_unrestricted(1)), 150),
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
        program_type: None,
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
        program_type: None,
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

/// The HSPs of `blastp -ungapped -comp_based_stats 0 -query fixture/q4.faa` (Query_1, 109 aa)
/// against six subjects of `PajaWSV.faa`: (subject, query start, subject start, length, score,
/// e-value, bit score, identities, positives). Subject 2 has two HSPs. The e-values and bit
/// scores as their 15 printed digits.
const BP_UNGAPPED_HSPS: [(u32, i32, i32, i32, i32, f64, f64, i32, i32); 7] = [
    (
        0,
        16,
        796,
        26,
        43,
        651688149968565e-15,
        229109485291396e-13,
        9,
        12,
    ),
    (
        1,
        25,
        146,
        69,
        40,
        171302022762502e-14,
        215166600089639e-13,
        15,
        26,
    ),
    (
        2,
        4,
        32,
        23,
        39,
        236412288450052e-14,
        210518971689053e-13,
        8,
        11,
    ),
    (
        2,
        74,
        252,
        12,
        35,
        857630603424144e-14,
        19192845808671e-12,
        6,
        8,
    ),
    (
        3,
        41,
        971,
        22,
        39,
        236412288450052e-14,
        210518971689053e-13,
        8,
        13,
    ),
    (
        4,
        9,
        1640,
        19,
        38,
        326270345375192e-14,
        205871343288467e-13,
        7,
        10,
    ),
    (
        5,
        15,
        961,
        19,
        38,
        326270345375192e-14,
        205871343288467e-13,
        6,
        14,
    ),
];

/// The ordinal ids (`-db`) and the record numbers (`-subject`) of the six subjects.
const BP_UNGAPPED_SUBJECTS: [u32; 6] = [4, 83, 24, 22, 82, 11];

/// `blastp -ungapped` against a database (`db_title`) or with `-subject` (`None`): one `diags`
/// Seq-align per subject, a Dense-diag per HSP, `unknown` strands on both rows, and no subject
/// remap of the strands (BLASTP `-subject` writes `unknown`, unlike its gapped Dense-segs).
fn blastp_ungapped_documents(format: SerialFormat, db_title: Option<&[u8]>) -> Vec<u8> {
    let run = RunRecord {
        program: Program::Blastp,
        db_title,
        gapped: false,
        program_type: None,
    };
    let subjects: Vec<SubjectRecord<'_>> = BP_UNGAPPED_SUBJECTS
        .iter()
        .map(|&number| match db_title {
            Some(_) => subject(
                SeqId::BlOrdId(crate::seqdb::Oid::of_unrestricted(number - 1)),
                1000,
            ),
            None => subject(SeqId::Local(SUBJECT_NAMES[number as usize]), 1000),
        })
        .collect();
    let hsps: Vec<HspRecord<'_>> = BP_UNGAPPED_HSPS
        .iter()
        .map(
            |&(subject, q_start, s_start, len, score, evalue, bit_score, ident, positives)| {
                HspRecord {
                    subject,
                    q_offset: q_start,
                    q_end: q_start + len,
                    q_left_out: 0,
                    s_offset: s_start,
                    s_end: s_start + len,
                    q_frame: 0,
                    s_frame: 0,
                    score,
                    num: 0,
                    evalue,
                    bit_score,
                    num_ident: ident,
                    num_positives: positives,
                    comp_adjust: 0,
                    ops: &[],
                }
            },
        )
        .collect();
    let mut buf = Vec::new();
    assert_eq!(
        write_query(
            format,
            &run,
            &query("Query_1", 109),
            &subjects,
            &hsps,
            &mut buf
        ),
        Ok(true)
    );
    buf
}

/// `Subject_N` of the record numbers of `BP_UNGAPPED_SUBJECTS`.
const SUBJECT_NAMES: [&str; 84] = {
    let mut names = [""; 84];
    names[4] = "Subject_4";
    names[83] = "Subject_83";
    names[24] = "Subject_24";
    names[22] = "Subject_22";
    names[82] = "Subject_82";
    names[11] = "Subject_11";
    names
};

/// The Dense-diag Seq-aligns of an ungapped BLASTP search, in the three encodings, against
/// NCBI BLAST+ 2.17.0 (`blastp -query q4.faa -subject s.faa|-db sdb -ungapped -comp_based_stats 0
/// -outfmt 8|9|12`; the score list without `comp_adjustment_method`, `unknown` strands).
#[test]
fn blastp_ungapped_dendiag_documents_equal_ncbi() {
    assert_same_bytes(
        &blastp_ungapped_documents(SerialFormat::AsnText, None),
        fixture!("blastp_ungapped_subject_fmt8.out"),
    );
    assert_same_bytes(
        &blastp_ungapped_documents(SerialFormat::AsnBinary, None),
        fixture!("blastp_ungapped_subject_fmt9.out"),
    );
    assert_same_bytes(
        &blastp_ungapped_documents(SerialFormat::Json, None),
        fixture!("blastp_ungapped_subject_fmt12.out"),
    );
    assert_same_bytes(
        &blastp_ungapped_documents(SerialFormat::AsnText, Some(b"in/s.faa")),
        fixture!("blastp_ungapped_db_fmt8.out"),
    );
    assert_same_bytes(
        &blastp_ungapped_documents(SerialFormat::AsnBinary, Some(b"in/s.faa")),
        fixture!("blastp_ungapped_db_fmt9.out"),
    );
    assert_same_bytes(
        &blastp_ungapped_documents(SerialFormat::Json, Some(b"in/s.faa")),
        fixture!("blastp_ungapped_db_fmt12.out"),
    );
}

/// The HSPs of `blastn -task megablast -ungapped -query mix_q.fa -subject mix_s.fa` (Query_1,
/// 1000 nt; Subject_1 with four HSPs, Subject_2 with two): (subject, query start in the query's
/// context, subject start, length, minus strand, score, e-value, bit score, identities). The
/// minus-strand HSPs start at `1000 - end` in the output; the context's offset is `end - length`.
const BN_UNGAPPED_HSPS: [(u32, i32, i32, i32, bool, i32, f64, f64, i32); 6] = [
    (
        0,
        99,
        99,
        302,
        false,
        302,
        177234191655786e-183,
        581339129949787e-12,
        302,
    ),
    (
        0,
        300,
        500,
        200,
        true,
        200,
        192662373350162e-124,
        3852249544608e-10,
        200,
    ),
    (
        0,
        749,
        799,
        152,
        false,
        152,
        116563164172467e-96,
        292935930701277e-12,
        152,
    ),
    (
        0,
        120,
        1000,
        140,
        false,
        140,
        102802205877556e-89,
        269863674761396e-12,
        140,
    ),
    (
        1,
        400,
        0,
        300,
        true,
        300,
        254754008593196e-182,
        577493753959806e-12,
        300,
    ),
    (
        1,
        600,
        380,
        200,
        false,
        200,
        192662373350162e-124,
        3852249544608e-10,
        200,
    ),
];

/// `blastn -ungapped` with `-subject`: the Dense-diags of plus-strand and minus-strand HSPs
/// (`minus, plus` strands; the query start of a minus HSP is `1000 - end`), in the three
/// encodings, against NCBI BLAST+ 2.17.0 (`-outfmt 8|9|12`).
#[test]
fn blastn_ungapped_dendiag_documents_equal_ncbi() {
    let run = RunRecord {
        program: Program::Blastn,
        db_title: None,
        gapped: false,
        program_type: None,
    };
    let subjects = [
        subject(SeqId::Local("Subject_1"), 1848),
        subject(SeqId::Local("Subject_2"), 1000),
    ];
    let hsps: Vec<HspRecord<'_>> = BN_UNGAPPED_HSPS
        .iter()
        .map(
            |&(subject, q_start, s_start, len, minus, score, evalue, bit_score, ident)| HspRecord {
                subject,
                q_offset: q_start,
                q_end: q_start + len,
                q_left_out: 0,
                s_offset: s_start,
                s_end: s_start + len,
                q_frame: if minus { -1 } else { 1 },
                s_frame: 1,
                score,
                num: 0,
                evalue,
                bit_score,
                num_ident: ident,
                num_positives: 0,
                comp_adjust: 0,
                ops: &[],
            },
        )
        .collect();
    for (format, name) in [
        (
            SerialFormat::AsnText,
            &fixture!("blastn_ungapped_mix_fmt8.out")[..],
        ),
        (
            SerialFormat::AsnBinary,
            &fixture!("blastn_ungapped_mix_fmt9.out")[..],
        ),
        (
            SerialFormat::Json,
            &fixture!("blastn_ungapped_mix_fmt12.out")[..],
        ),
    ] {
        let mut buf = Vec::new();
        assert_eq!(
            write_query(
                format,
                &run,
                &query("Query_1", 1000),
                &subjects,
                &hsps,
                &mut buf
            ),
            Ok(true)
        );
        assert_same_bytes(&buf, name);
    }
}

/// A query without an HSP writes nothing in an ungapped search (`results.HasAlignments()`).
#[test]
fn an_ungapped_query_without_hsps_writes_nothing() {
    let run = RunRecord {
        program: Program::Blastn,
        db_title: None,
        gapped: false,
        program_type: None,
    };
    let mut buf = Vec::new();
    assert_eq!(
        write_query(
            SerialFormat::AsnText,
            &run,
            &query("Query_1", 10),
            &[],
            &[],
            &mut buf
        ),
        Ok(false)
    );
    assert!(buf.is_empty());
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

/// A byte outside the visible range is an error of the ASN.1 text encoder (its callers refuse
/// such strings, or write the document again through `Traced`).
#[test]
fn asn_text_records_a_non_visible_byte() {
    let mut buf = Vec::new();
    let mut enc = super::asn_text::AsnTextEncoder::new(&mut buf);
    enc.visible_str(b"a\tb");
    assert_eq!(enc.error(), Some(SerialError::NonVisibleChar(b'\t')));
    assert_eq!(enc.take_replaced(), Vec::new());
}

/// Through `Traced`, each byte outside the visible range is `#` and a `ReplacedChar` at the
/// line where NCBI checks it: before the byte's wrap (NCBI 2.17.0 convert2blastmask, probe RG
/// K/k4/probe case c_txt_w62: the two bytes of the first `é` at the 63rd and 64th place of a
/// line that wraps at the second are both on line 4, the next `é` on line 5). A flush before
/// the end keeps the document.
#[test]
fn asn_text_with_frames_replaces_each_byte_at_its_line() {
    // `  m "` is 5 columns: the wrap comes before the 75th byte of the string (column 79).
    let string = [
        b"a".repeat(73),
        "\u{e9}\u{e9}".as_bytes().to_vec(),
        b"b".to_vec(),
    ]
    .concat();
    let mut sink = Vec::new();
    let mut buf = Vec::new();
    let mut enc = Traced::new(super::asn_text::AsnTextEncoder::new(&mut buf), 0);
    enc.begin_document("T");
    enc.begin_sequence("T");
    enc.member("n", 0);
    enc.int(1);
    enc.end_member();
    enc.flush_to(&mut sink, false).unwrap();
    enc.member("m", 1);
    enc.visible_str(&string);
    enc.end_member();
    enc.member("l", 2);
    enc.begin_list(false);
    enc.element();
    enc.visible_str(b"x");
    enc.element();
    enc.visible_str(b"\x7fy");
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
    assert_eq!(enc.error(), None);
    let replaced = enc.take_replaced_chars();
    enc.flush_to(&mut sink, true).unwrap();
    let text = String::from_utf8(sink).unwrap();
    assert_eq!(
        text,
        format!(
            "T ::= {{\n  n 1,\n  m \"{}\n####b\",\n  l {{\n    \"x\",\n    \"#y\"\n  }}\n}}\n",
            "a".repeat(73)
        )
    );
    let at = |byte: u8, line: u64, trace: &str, string: &[u8]| ReplacedChar {
        trace: trace.to_string(),
        byte,
        position: StreamPosition::Line(line),
        string: std::sync::Arc::from(string),
    };
    assert_eq!(
        replaced,
        vec![
            at(0xc3, 3, "T.m", &string),
            at(0xa9, 3, "T.m", &string),
            at(0xc3, 4, "T.m", &string),
            at(0xa9, 4, "T.m", &string),
            at(0x7f, 7, "T.l.E", b"\x7fy"),
        ]
    );
    assert_eq!(
        replaced[4].message(),
        b"T.l.E\nBad char [0x7F] in string at line 7\n\x7fy".to_vec()
    );
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

/// LOSAT's refusal of a database for the Seq-annot formats (`check_database`): a database named
/// `n/a` or `SRA` (any case), whatever its title.
#[test]
fn databases_that_ncbi_describes_otherwise_are_refused() {
    use super::seqalign::check_database;
    for name in ["n/a", "N/A", "sra", "SRA"] {
        for format in [
            SerialFormat::AsnText,
            SerialFormat::AsnBinary,
            SerialFormat::Json,
        ] {
            let error = check_database(format, name, Program::Blastn).unwrap_err();
            assert!(
                error
                    .to_string()
                    .ends_with("is not supported by LOSAT's BLASTN"),
                "{error}"
            );
        }
    }
    for format in [
        SerialFormat::AsnText,
        SerialFormat::AsnBinary,
        SerialFormat::Json,
    ] {
        assert!(check_database(format, "db/x", Program::Blastp).is_ok());
        assert!(check_database(format, "db/nash", Program::Blastp).is_ok());
    }
}

/// One stream that a test's report and diagnostics both write to (`2>&1`).
#[derive(Clone, Default)]
struct Combined(std::rc::Rc<std::cell::RefCell<Vec<u8>>>);

impl std::io::Write for Combined {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.borrow_mut().extend_from_slice(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl crate::report::ostream::ReportWrite for Combined {}

/// The documents of `blastp_db_documents` written by `write_query_to` with the database title
/// `title`: the report alone, the posts alone, and both in one stream.
fn blastp_db_documents_posted(format: SerialFormat, title: &[u8]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let run = RunRecord {
        program: Program::Blastp,
        db_title: Some(title),
        gapped: true,
        program_type: None,
    };
    let subjects = [
        subject(SeqId::BlOrdId(crate::seqdb::Oid::of_unrestricted(0)), 137),
        subject(SeqId::BlOrdId(crate::seqdb::Oid::of_unrestricted(1)), 150),
    ];
    let hsps = [
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
    let write = |writer: &mut dyn crate::report::ostream::ReportWrite,
                 diagnostics: &mut dyn std::io::Write| {
        let mut buf = Vec::new();
        for id in ["Query_1", "Query_2"] {
            super::write_query_to(
                format,
                &run,
                "blastp",
                &query(id, 140),
                &subjects,
                &hsps,
                &mut buf,
                writer,
                diagnostics,
            )
            .unwrap();
        }
    };
    let (mut report, mut posts) = (Vec::new(), Vec::new());
    write(&mut report, &mut posts);
    let combined = Combined::default();
    write(&mut combined.clone(), &mut combined.clone());
    let both = combined.0.borrow().clone();
    (report, posts, both)
}

/// A database title with bytes outside 0x20..0x7E in the ASN.1 text and binary Seq-annots: `#`
/// for each, and NCBI's post of each byte of each query's document, whose position counts from
/// that document's start (one object stream per document), with the program's prefix; a short
/// document is still in NCBI's buffer when it posts, so the posts come before its bytes
/// (NCBI 2.17.0, `blastn -db` of a title `t\xe9tle db`: `Bad char [0xE9] in string at line 25`,
/// `byte 197`, before each document). The JSON stream posts nothing (it escapes the bytes below
/// 0x20 and writes these as they are).
#[test]
fn a_database_title_outside_the_visible_range_is_written_with_ncbis_posts() {
    let title = b"caf\xc3\xa9 \x7f";
    let shown = b"caf## #";
    for format in [SerialFormat::AsnText, SerialFormat::AsnBinary] {
        let (report, posts, both) = blastp_db_documents_posted(format, title);
        // The report is the plain writer's with `#` in place of each such byte.
        let (plain, no_posts, _) = blastp_db_documents_posted(format, shown);
        assert!(no_posts.is_empty());
        assert_same_bytes(&report, &plain);
        let document = &plain[..plain.len() / 2];
        let at = document
            .windows(shown.len())
            .position(|window| window == shown)
            .unwrap();
        let position = |offset: usize| match format {
            SerialFormat::AsnText => format!(
                "line {}",
                1 + document[..at].iter().filter(|&&byte| byte == b'\n').count()
            ),
            _ => format!("byte {}", at + offset),
        };
        let mut post = Vec::new();
        for (byte, offset) in [(0xC3u8, 3), (0xA9, 4), (0x7F, 6)] {
            post.extend_from_slice(
                format!(
                    "Error: (803.7) [blastp] Seq-annot.desc.E.user.data.E.label.str\nBad char [0x{byte:X}] in string at {}\n",
                    position(offset)
                )
                .as_bytes(),
            );
            post.extend_from_slice(title);
            post.push(b'\n');
        }
        assert_eq!(posts, [&post[..], &post[..]].concat());
        // The two documents differ only in the query's id (`Query_1`, `Query_2`).
        assert_eq!(
            both,
            [&post[..], document, &post[..], &plain[document.len()..]].concat(),
            "{format:?}"
        );
    }
    let (report, posts, _) = blastp_db_documents_posted(SerialFormat::Json, title);
    assert!(posts.is_empty());
    assert!(report.windows(title.len()).any(|window| window == title));
}

/// The Std-segs of a translated HSP are one per segment of its records with either subject: a
/// `-subject` report's records come mapped (`mapped_script`), so the writer drops none.
#[test]
fn std_segs_are_written_for_every_segment_of_the_records() {
    let ops = [Del(2), Sub(10), Ins(1), Sub(5), Ins(3)];
    let one = hsp(0, (0, 0), (1, 0), 30, 1.0, 1.0, (10, 12, 0), &ops);
    let segments = |db_title: Option<&'static [u8]>| {
        let run = RunRecord {
            program: Program::Blastx,
            db_title,
            gapped: true,
            program_type: None,
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
    assert_eq!(segments(None), (5, 3));
    assert_eq!(segments(Some(b"db")), (5, 3));
}

/// The mapper of a `-subject` search leaves out the first and the last Std-seg when it has a gap
/// row, once each (seq_align_mapper_base.cpp:1568-1580): `mapped_script` gives the kept operations
/// and the letters left out at each end (the query's, which the coverage still counts, in
/// `query_left_out`); an X-drop script (a substitution at each end) and an ungapped HSP (no
/// script) keep everything.
#[test]
fn mapped_script_leaves_out_an_end_segment_with_a_gap_row_once() {
    use super::seqalign::{mapped_script, MappedScript};
    let ops = [Del(2), Sub(10), Ins(1), Sub(5), Ins(3)];
    assert_eq!(
        mapped_script(&ops),
        MappedScript {
            kept: 1..4,
            start: [0, 2],
            end: [3, 0],
        }
    );
    assert!(mapped_script(&ops).drops(&ops));
    // The coverage still counts the query letters left out: the end `Ins(3)`.
    assert_eq!(mapped_script(&ops).query_left_out(), 3);
    // The Smith-Waterman alignment of the probe: 95 subject letters before the query's first.
    let sw = [Del(95), Sub(190)];
    assert_eq!(
        mapped_script(&sw),
        MappedScript {
            kept: 1..2,
            start: [0, 95],
            end: [0, 0],
        }
    );
    assert_eq!(mapped_script(&sw).query_left_out(), 0);
    // Only the first segment goes, though the next one has a gap row too.
    let two = [Ins(4), Del(2), Sub(7)];
    assert_eq!(
        mapped_script(&two),
        MappedScript {
            kept: 1..3,
            start: [4, 0],
            end: [0, 0],
        }
    );
    assert_eq!(mapped_script(&two).query_left_out(), 4);
    let xdrop = [Sub(10), Ins(1), Sub(5), Del(2), Sub(3)];
    let kept = mapped_script(&xdrop);
    assert_eq!(kept.kept, 0..5);
    assert!(!kept.drops(&xdrop));
    assert_eq!((kept.start, kept.end), ([0, 0], [0, 0]));
    assert_eq!(kept.query_left_out(), 0);
    assert_eq!(mapped_script(&[]).kept, 0..0);
    // A script of one gap segment (never made) keeps nothing.
    assert_eq!(mapped_script(&[Del(3)]).kept, 1..1);
}

/// The end gaps of each row (`row_end_gaps`): the row's gap operations before its first letter
/// and after its last (a `Del` in the query's row, an `Ins` in the subject's), which NCBI's
/// `GetAlignLengths` and `x_GetNumGaps` leave out of the query row's counts; the inner gaps and
/// the other row's are not.
#[test]
fn row_end_gaps_are_the_rows_gap_operations_outside_its_letters() {
    use super::seqalign::{gapped_row_range, row_end_gaps, AlignRow, RowEndGaps, RowRangeError};
    // A gap start in the query's row (the probe's Smith-Waterman alignment).
    let query_start = [Del(95), Sub(190)];
    assert_eq!(
        row_end_gaps(&query_start, AlignRow::Query),
        RowEndGaps {
            leading: 95,
            trailing: 0,
            opens: 1,
        }
    );
    assert_eq!(
        row_end_gaps(&query_start, AlignRow::Subject),
        RowEndGaps::default()
    );
    assert_eq!(
        gapped_row_range(&query_start, AlignRow::Query),
        Some(RowRangeError::SeqidsNotConsistent)
    );
    assert_eq!(gapped_row_range(&query_start, AlignRow::Subject), None);
    // A gap start and a gap end in the subject's row: the query row has letters in the first and
    // the last operation, so its deletions are inner gaps.
    let ends = [Ins(3), Del(2), Sub(10), Del(1), Sub(4), Del(5), Ins(6)];
    assert_eq!(
        row_end_gaps(&ends, AlignRow::Subject),
        RowEndGaps {
            leading: 3,
            trailing: 6,
            opens: 2,
        }
    );
    assert_eq!(row_end_gaps(&ends, AlignRow::Query), RowEndGaps::default());
    assert_eq!(row_end_gaps(&ends, AlignRow::Subject).span(31), (3, 25));
    assert_eq!(
        gapped_row_range(&ends, AlignRow::Subject),
        Some(RowRangeError::SeqidsNotConsistent)
    );
    assert_eq!(gapped_row_range(&ends, AlignRow::Query), None);
    // A gap end in the query's row (two operations), after an inner gap that stays.
    let query_end = [Sub(10), Del(1), Sub(4), Del(5), Del(2)];
    assert_eq!(
        row_end_gaps(&query_end, AlignRow::Query),
        RowEndGaps {
            leading: 0,
            trailing: 7,
            opens: 2,
        }
    );
    assert_eq!(row_end_gaps(&query_end, AlignRow::Query).columns(), 7);
    assert_eq!(row_end_gaps(&query_end, AlignRow::Query).span(22), (0, 15));
    assert_eq!(
        row_end_gaps(&query_end, AlignRow::Subject),
        RowEndGaps::default()
    );
    assert_eq!(gapped_row_range(&query_end, AlignRow::Query), None);
    // An X-drop script has none: its inner gaps stay in the counts.
    let xdrop = [Sub(10), Ins(1), Sub(5), Del(2), Sub(3)];
    for row in [AlignRow::Query, AlignRow::Subject] {
        assert_eq!(row_end_gaps(&xdrop, row), RowEndGaps::default());
        assert_eq!(row_end_gaps(&xdrop, row).span(21), (0, 21));
    }
    assert_eq!(row_end_gaps(&[], AlignRow::Query), RowEndGaps::default());
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191,226-230
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1568-1580
/// `reported_std_seg_row_range` / `reported_std_seg_rows_range`: `empty_query_rewritten` clears
/// `RowRangeError::Empty` on the query row alone; any leading gap left in the reported edit
/// script (`SeqidsNotConsistent`, e.g. after `-subject` dropped one end segment of `[Del, Ins]` or
/// `[Ins, Del]`) is reported for `-db` and `-subject` alike.
#[test]
fn reported_std_seg_row_range_handles_empty_rewrite_and_mapped_leading_gaps() {
    use super::seqalign::{
        mapped_script, reported_std_seg_row_range, reported_std_seg_rows_range, AlignRow,
        RowRangeError,
    };
    // Empty ungapped rows: `-subject` rewrites the query row only.
    assert_eq!(
        reported_std_seg_row_range(4, 4, &[], AlignRow::Query, false),
        Some(RowRangeError::Empty)
    );
    assert_eq!(
        reported_std_seg_row_range(4, 4, &[], AlignRow::Query, true),
        None
    );
    assert_eq!(
        reported_std_seg_row_range(4, 4, &[], AlignRow::Subject, true),
        Some(RowRangeError::Empty)
    );
    assert_eq!(
        reported_std_seg_rows_range(4, 4, &[], true),
        Some(RowRangeError::Empty)
    );

    // `[Del(2), Ins(1), Sub(10)]` mapped once -> `[Ins(1), Sub(10)]` (subject starts with gap).
    let del_ins = [Del(2), Ins(1), Sub(10)];
    let mapped = &del_ins[mapped_script(&del_ins).kept];
    assert_eq!(mapped, &[Ins(1), Sub(10)]);
    for rewritten in [false, true] {
        assert_eq!(
            reported_std_seg_row_range(0, 11, mapped, AlignRow::Query, rewritten),
            None
        );
        assert_eq!(
            reported_std_seg_row_range(0, 11, mapped, AlignRow::Subject, rewritten),
            Some(RowRangeError::SeqidsNotConsistent)
        );
        assert_eq!(
            reported_std_seg_rows_range(0, 11, mapped, rewritten),
            Some(RowRangeError::SeqidsNotConsistent)
        );
    }

    // `[Ins(1), Del(2), Sub(10)]` mapped once -> `[Del(2), Sub(10)]` (query starts with gap).
    let ins_del = [Ins(1), Del(2), Sub(10)];
    let mapped = &ins_del[mapped_script(&ins_del).kept];
    assert_eq!(mapped, &[Del(2), Sub(10)]);
    for rewritten in [false, true] {
        assert_eq!(
            reported_std_seg_row_range(1, 11, mapped, AlignRow::Query, rewritten),
            Some(RowRangeError::SeqidsNotConsistent)
        );
        assert_eq!(
            reported_std_seg_row_range(1, 11, mapped, AlignRow::Subject, rewritten),
            None
        );
        assert_eq!(
            reported_std_seg_rows_range(1, 11, mapped, rewritten),
            Some(RowRangeError::SeqidsNotConsistent)
        );
    }

    // `[Del(2), Sub(10)]` mapped once -> `[Sub(10)]` (no leading gap remains).
    let single = [Del(2), Sub(10)];
    let mapped = &single[mapped_script(&single).kept];
    assert_eq!(mapped, &[Sub(10)]);
    for rewritten in [false, true] {
        assert_eq!(reported_std_seg_rows_range(0, 10, mapped, rewritten), None);
    }
}

/// The query frame of a translated query's section (`HspStrands::master_frame`): the HSP's own,
/// or, where the query row starts with a gap, the minus frame of the record's length (NCBI
/// BLAST+ 2.17.0 shows queries of 570, 571 and 572 nucleotides under -1, -2 and -3).
#[test]
fn a_query_row_that_starts_with_a_gap_takes_the_frame_of_the_records_length() {
    use crate::report::query_anchored::HspStrands;
    for frame in [1, 2, 3, -1, -2, -3] {
        assert_eq!(HspStrands::master_frame(frame, false, 570), frame);
    }
    assert_eq!(HspStrands::master_frame(1, true, 570), -1);
    assert_eq!(HspStrands::master_frame(2, true, 571), -2);
    assert_eq!(HspStrands::master_frame(-1, true, 572), -3);
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
        program_type: None,
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
            program_type: None,
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

/// The bytes replaced in one string share one copy of it, and the lines go on from the last
/// replaced byte (a string of 20,000 bytes, every other one DEL: 10,000 posts of NCBI, each with
/// the whole string).
#[test]
fn replaced_bytes_of_a_long_string_share_it() {
    let string: Vec<u8> = (0..20_000)
        .map(|at| if at % 2 == 0 { b'a' } else { 0x7f })
        .collect();
    let mut buf = Vec::new();
    let mut enc = Traced::new(super::asn_text::AsnTextEncoder::new(&mut buf), 0);
    enc.begin_document("T");
    enc.begin_sequence("T");
    enc.member("m", 0);
    enc.visible_str(&string);
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
    let replaced = enc.take_replaced_chars();
    assert_eq!(replaced.len(), 10_000);
    assert!(replaced
        .iter()
        .all(|byte| std::sync::Arc::ptr_eq(&byte.string, &replaced[0].string)));
    assert_eq!(&replaced[0].string[..], &string[..]);
    // `  m "` and 73 more bytes fill line 2 to column 78; each later line holds 78 bytes.
    let lines: Vec<u64> = replaced
        .iter()
        .map(|byte| match byte.position {
            StreamPosition::Line(line) => line,
            StreamPosition::Byte(_) => unreachable!(),
        })
        .collect();
    assert_eq!(lines[0], 2);
    assert!(lines.windows(2).all(|pair| pair[0] <= pair[1]));
    let text_lines = buf.iter().filter(|&&byte| byte == b'\n').count() as u64;
    assert!(lines[9_999] <= text_lines);
}

/// NCBI's insertions of a Seq-annot into the report's stream (`document_pieces`, the writes of its
/// object stream's buffer) cover the document, and for JSON the model of the stream's calls
/// (`Traced`) and the one taken back from the bytes (`out_buffer::json_pieces`, which the side
/// files and outfmt 15 use) give the same pieces, over a document of several buffers.
#[test]
fn the_pieces_of_a_seq_annot_cover_it_and_both_json_models_agree() {
    let run = RunRecord {
        program: Program::Blastp,
        db_title: Some(b"a database"),
        gapped: true,
        program_type: None,
    };
    let subjects: Vec<_> = (0..40)
        .map(|i| {
            subject(
                SeqId::BlOrdId(crate::seqdb::Oid::of_unrestricted(i)),
                137 + i as i32,
            )
        })
        .collect();
    let hsps: Vec<_> = (0..40)
        .map(|i| {
            hsp(
                i,
                (i as i32, i as i32),
                (0, 0),
                652,
                5.30325651018409e-94,
                255.758355778158,
                (125, 130, 2),
                &BP_Q1_HSP2,
            )
        })
        .collect();
    let query = query("Query_1", 140);
    for format in [
        SerialFormat::AsnText,
        SerialFormat::AsnBinary,
        SerialFormat::Json,
    ] {
        let mut bytes = Vec::new();
        assert!(super::write_query(format, &run, &query, &subjects, &hsps, &mut bytes).unwrap());
        let pieces = super::seqalign::query_pieces(format, &run, &query, &subjects, &hsps);
        assert!(pieces.len() > 2, "{format:?} {}", bytes.len());
        assert_eq!(pieces.iter().sum::<usize>(), bytes.len(), "{format:?}");
        if format == SerialFormat::AsnText {
            // The separator after the document is a write of its own (`WriteSeparator`).
            assert_eq!(pieces.last(), Some(&1));
        }
        if format == SerialFormat::Json {
            assert_eq!(
                crate::report::serial::out_buffer::json_pieces(&bytes),
                pieces
            );
        }
    }
}

// Line DBID pfoff-4: BLAST XML's `catch` of an HSP whose subject's lookup throws writes the
// Seq-align's identities as `GetAlnScores` reads them: the count, or -1 where the Seq-align has
// no `num_ident` (a count below 0; psiblast's 0, which NCBI sets to -1 first), the one rule of
// the Seq-annot formats too (`seq_align_num_ident`).
#[test]
fn the_identities_of_a_seq_align_as_ncbi_reads_them() {
    use super::seqalign::{read_num_ident, seq_align_num_ident};
    use crate::report::record::ProgramType;
    let psiblast = Some(ProgramType::PsiBlast);
    assert_eq!(read_num_ident(66, None), 66);
    assert_eq!(read_num_ident(0, None), 0);
    assert_eq!(read_num_ident(-1, None), -1);
    assert_eq!(read_num_ident(0, psiblast), -1);
    assert_eq!(read_num_ident(5, psiblast), 5);
    assert_eq!(read_num_ident(0, Some(ProgramType::PsiTblastn)), 0);
    assert_eq!(read_num_ident(i64::from(i32::MAX) + 1, None), i32::MAX);
    // The score the Seq-annot formats write: none where XML's `catch` reads -1.
    assert_eq!(seq_align_num_ident(0, psiblast), None);
    assert_eq!(seq_align_num_ident(0, Some(ProgramType::RpsBlast)), Some(0));
    assert_eq!(seq_align_num_ident(-1, None), None);
}

/// A PSSM search's Seq-annot names its program type (`Blast Type`: psiblast 10, psitblastn 11,
/// blast_format.cpp:666-673), and psiblast's Seq-aligns leave out a `num_ident` of 0
/// (blast_seqalign.cpp:1490-1496); blastp's write it (`score_list_follows_ncbi_conditions`).
#[test]
fn pssm_searches_name_their_type_and_psiblast_omits_zero_identities() {
    use crate::report::record::ProgramType;
    let subjects = [subject(SeqId::Local("Subject_1"), 100)];
    let ops = [Sub(10)];
    let write = |program_type: Option<ProgramType>, program: Program, num_ident: i32| {
        let run = RunRecord {
            program,
            db_title: None,
            gapped: true,
            program_type,
        };
        let frames = if program == Program::Tblastn {
            (0, 1)
        } else {
            (0, 0)
        };
        let one = hsp(0, (0, 0), frames, 30, 1e-5, 20.5, (num_ident, 0, 0), &ops);
        let mut buf = Vec::new();
        write_query(
            SerialFormat::AsnText,
            &run,
            &query("Query_1", 10),
            &subjects,
            &[one],
            &mut buf,
        )
        .unwrap();
        String::from_utf8(buf).unwrap()
    };
    let psiblast = write(Some(ProgramType::PsiBlast), Program::Blastp, 0);
    assert!(psiblast.contains("str \"psiblast\""), "{psiblast}");
    assert!(psiblast.contains("int 10"), "{psiblast}");
    assert!(!psiblast.contains("num_ident"), "{psiblast}");
    let psiblast_ident = write(Some(ProgramType::PsiBlast), Program::Blastp, 7);
    assert!(psiblast_ident.contains("num_ident"), "{psiblast_ident}");
    let blastp = write(None, Program::Blastp, 0);
    assert!(blastp.contains("str \"blastp\""), "{blastp}");
    assert!(blastp.contains("num_ident"), "{blastp}");
    let tblastn = write(Some(ProgramType::PsiTblastn), Program::Tblastn, 0);
    assert!(tblastn.contains("str \"psitblastn\""), "{tblastn}");
    assert!(tblastn.contains("int 11"), "{tblastn}");
    assert!(tblastn.contains("num_ident"), "{tblastn}");
}

/// The Dense-seg of a translated HSP's rows (`translated_row_segments`, NCBI's
/// `CreateDensegFromStdseg`) gives each row the strand of its first segment with one in every
/// segment: a minus-frame query row that starts with a gap (whose first Std-seg has an empty
/// location) is on the minus strand, so outfmt 0 counts its positions down; the protein row has
/// no strand, and an X-drop script's rows keep theirs.
#[test]
fn a_translated_rows_leading_gap_takes_the_rows_strand() {
    use super::seqalign::{translated_row_segments, Segments, Strand};
    let strands = |frame: i32, ops: Vec<GapEditOp>| {
        let hit = crate::common::Hit {
            identity: 100.0,
            length: 285,
            mismatch: 0,
            gapopen: 1,
            q_start: 1,
            q_end: 570,
            s_start: 1,
            s_end: 285,
            e_value: 1e-100,
            bit_score: 353.0,
            num_ident: 190,
            query_frame: frame,
            query_length: 570,
            q_idx: 0,
            s_idx: 0,
            raw_score: 906,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: true,
            gap_info: Some(ops),
            num_positives: 190,
        };
        let mut segments = Segments::default();
        translated_row_segments(&hit, Some((frame as i8, 570)), None, &mut segments);
        segments.strands().to_vec()
    };
    let (minus, plus, none) = (Strand::Minus, Strand::Plus, Strand::Unknown);
    assert_eq!(
        strands(-1, vec![Del(95), Sub(190)]),
        [minus, none, minus, none]
    );
    assert_eq!(
        strands(1, vec![Del(95), Sub(190)]),
        [plus, none, plus, none]
    );
    assert_eq!(
        strands(-1, vec![Sub(100), Ins(5), Sub(85)]),
        [minus, none, minus, none, minus, none]
    );
}
