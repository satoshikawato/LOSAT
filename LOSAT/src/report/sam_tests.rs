use super::*;
use crate::common::GapEditOp;
use crate::report::seqid::SeqId;

fn g(value: f64) -> String {
    let mut out = Vec::new();
    push_g(&mut out, value);
    String::from_utf8(out).unwrap()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2316-2340 (`%g`).
#[test]
fn g_prints_six_significant_digits_as_c_does() {
    assert_eq!(g(0.0), "0");
    assert_eq!(g(1e-180), "1e-180");
    assert_eq!(g(2.60916e-178), "2.60916e-178");
    assert_eq!(g(0.8216932), "0.821693");
    assert_eq!(g(608.6682), "608.668");
    assert_eq!(g(55395.12), "55395.1");
    assert_eq!(g(100000.0), "100000");
    assert_eq!(g(999999.5), "1e+06");
    assert_eq!(g(1e6), "1e+06");
    assert_eq!(g(123456789.0), "1.23457e+08");
    assert_eq!(g(0.0001), "0.0001");
    assert_eq!(g(0.00001), "1e-05");
    assert_eq!(g(0.5), "0.5");
    assert_eq!(g(2.0), "2");
    assert_eq!(g(1.5e-5), "1.5e-05");
    assert_eq!(g(-0.25), "-0.25");
    assert_eq!(g(f64::INFINITY), "inf");
}

#[test]
fn the_specification_words_are_exact() {
    let parse = SamSpec::parse;
    assert_eq!(parse(""), SamSpec::default());
    assert_eq!(
        parse("SQ"),
        SamSpec {
            subject_as_reference: false,
            sequence: true
        }
    );
    assert_eq!(
        parse("SR SQ"),
        SamSpec {
            subject_as_reference: true,
            sequence: true
        }
    );
    for plain in ["sq", "foo SQ,SR", "std", "SQSR", "S R"] {
        assert_eq!(parse(plain), SamSpec::default(), "{plain}");
    }
    assert!(parse("foo SQ").sequence);
    assert!(!parse("SQ  SQ").subject_as_reference);
}

#[test]
fn a_reverse_complement_keeps_the_ambiguity_codes() {
    let mut letters = b"ACGTRYKMBVDHNSW".to_vec();
    reverse_complement(&mut letters);
    assert_eq!(letters, b"WSNDHBVKMRYACGT");
}

/// A document of one query with the given HSP over a subject of length 100; the query
/// has 60 letters.
fn document(
    spec: SamSpec,
    cl: Option<&[u8]>,
    hsp: HspRecord<'_>,
    letters: &mut dyn SamLetters,
) -> String {
    let run = RunRecord {
        program: crate::blastinput::app::Program::Blastn,
        rps: false,
        db_title: None,
        gapped: true,
    };
    let query = SamSequence {
        record: SamId::Query(QueryRecord {
            id: SeqId::Local("Query_1"),
            searched_length: 60,
            offset: 0,
        }),
        key: 0,
        index: 0,
        length: 60,
    };
    let subject = SamSequence {
        record: SamId::Subject(SubjectRecord {
            id: SeqId::Local("Subject_1"),
            searched_length: 100,
            offset: 0,
            use_this_seqid: &[],
        }),
        key: 0,
        index: 0,
        length: 100,
    };
    let mut document = SamDocument::new(spec, "2.17.0+", cl.map(<[u8]>::to_vec));
    document
        .query(&run, &query, &[subject], &[hsp], letters)
        .unwrap();
    let mut out = Vec::new();
    document.finish(&mut out).unwrap();
    String::from_utf8(out).unwrap()
}

fn hsp(q: (i32, i32), s: (i32, i32), frame: i8, ident: i32, ops: &[GapEditOp]) -> HspRecord<'_> {
    HspRecord {
        subject: 0,
        q_offset: q.0,
        q_end: q.1,
        s_offset: s.0,
        s_end: s.1,
        q_frame: frame,
        s_frame: 1,
        score: 40,
        num: 0,
        evalue: 2.5e-7,
        bit_score: 75.8,
        num_ident: ident,
        num_positives: 0,
        comp_adjust: 0,
        ops,
    }
}

struct Fixed;

impl SamLetters for Fixed {
    fn letters(
        &mut self,
        row: SamRow,
        _: usize,
        from: usize,
        to: usize,
        out: &mut Vec<u8>,
    ) -> Result<()> {
        let alphabet: &[u8] = if row == SamRow::Query {
            b"AACCGGTTAC"
        } else {
            b"ACGTRYACGT"
        };
        out.extend((from..to).map(|at| alphabet[at % alphabet.len()]));
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:179-298 (AddRow) and
// cigar_formatter.cpp:157-352 (x_FormatLine): the target is the subject, the reference the
// query; a plus-strand HSP with a gap of 3 in the subject (D) and one of 2 in the query (I).
#[test]
fn a_plus_hsp_has_the_subject_as_target_and_clips_it() {
    let ops = [
        GapEditOp::Sub(10),
        GapEditOp::Ins(3),
        GapEditOp::Sub(20),
        GapEditOp::Del(2),
        GapEditOp::Sub(5),
    ];
    // The query consumes 10 + 3 + 20 + 5 = 38 letters, the subject 10 + 20 + 2 + 5 = 37.
    let hsp = hsp((4, 42), (20, 57), 1, 30, &ops);
    let text = document(SamSpec::default(), Some(b"blastn -x "), hsp, &mut NoLetters);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "@HD\tVN:1.2\tSO:coordinate\tGO:reference");
    assert_eq!(lines[1], "@SQ\tSN:Query_1\tLN:60");
    assert_eq!(lines[2], "@PG\tID:0\tVN:2.17.0+\tCL:blastn -x \tPN:blastn");
    let fields: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(fields[0], "Subject_1");
    assert_eq!(fields[1], "0");
    assert_eq!(fields[2], "Query_1");
    assert_eq!(fields[3], "5");
    assert_eq!(fields[4], "255");
    // The clips come from the subject range [20, 57) of a length of 100.
    assert_eq!(fields[5], "20H10M3D20M2I5M43H");
    assert_eq!(fields[9], "*");
    assert_eq!(fields[11], "AS:i:40");
    assert_eq!(fields[12], "EV:f:2.5e-07");
    assert_eq!(fields[13], "NM:i:5");
    assert_eq!(fields[14], "PI:f:85.71");
    assert_eq!(fields[15], "BS:f:75.8");
    assert_eq!(lines.len(), 4);
}

#[test]
fn sq_writes_the_target_letters_and_sr_swaps_the_roles() {
    let ops = [GapEditOp::Sub(10), GapEditOp::Ins(3), GapEditOp::Sub(10)];
    let spec = SamSpec::parse("SQ SR");
    let hsp = hsp((4, 27), (20, 40), 1, 20, &ops);
    let text = document(spec, None, hsp, &mut Fixed);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "@HD\tVN:1.2\tGO:query");
    assert_eq!(lines[1], "@SQ\tSN:Subject_1\tLN:100");
    // No CL tag without a command line.
    assert_eq!(lines[2], "@PG\tID:0\tVN:2.17.0+\tPN:blastn");
    let fields: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(fields[0], "Query_1");
    assert_eq!(fields[2], "Subject_1");
    assert_eq!(fields[3], "21");
    // The target is the query: its three letters have no subject letters (I).
    assert_eq!(fields[5], "4H10M3I10M33H");
    let expected: String = (4..27).map(|at| b"AACCGGTTAC"[at % 10] as char).collect();
    assert_eq!(fields[9], expected);
}

#[test]
fn a_minus_hsp_is_reversed_and_complemented() {
    let ops = [GapEditOp::Sub(12)];
    // Query row on the minus strand: the HSP covers the query context offsets 8..20 of 60,
    // so the plus-strand range is 40..52.
    let hsp = hsp((8, 20), (30, 42), -1, 12, &ops);
    let text = document(SamSpec::parse("SQ"), None, hsp, &mut Fixed);
    let lines: Vec<&str> = text.lines().collect();
    let fields: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(fields[1], "16");
    assert_eq!(fields[3], "41");
    // Subject range [30, 42): reverse clip at the back (30H) and the front (58H).
    assert_eq!(fields[5], "58H12M30H");
    assert_eq!(fields[13], "NM:i:0");
    assert_eq!(fields[14], "PI:f:100.00");
    assert_eq!(fields[9].len(), 12);
}

fn query_seq(index: usize, id: &'static str) -> SamSequence<'static> {
    SamSequence {
        record: SamId::Query(QueryRecord {
            id: SeqId::Local(id),
            searched_length: 60,
            offset: 0,
        }),
        key: index as u64,
        index,
        length: 60,
    }
}

fn subject_seq(index: usize, id: &'static str) -> SamSequence<'static> {
    SamSequence {
        record: SamId::Subject(SubjectRecord {
            id: SeqId::Local(id),
            searched_length: 100,
            offset: 0,
            use_this_seqid: &[],
        }),
        key: index as u64,
        index,
        length: 100,
    }
}

fn run_record() -> RunRecord<'static> {
    RunRecord {
        program: crate::blastinput::app::Program::Blastn,
        rps: false,
        db_title: None,
        gapped: true,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1401-1430
// (`SortHspByMasterStartAscending`): the lowest query coordinate, then the e-value; `list::sort`
// is stable, so HSPs equal in both keep the order of the Seq-align.
#[test]
fn the_default_mode_sorts_a_querys_records_stably_by_start_then_evalue() {
    let ops = [GapEditOp::Sub(10)];
    let mut hsps = Vec::new();
    // (query start, e-value, score): input order 0..4
    for (start, evalue, score) in [
        (30, 1e-5, 11),
        (10, 1e-3, 12),
        (10, 1e-9, 13),
        (30, 1e-5, 14),
        (5, 1.0, 15),
    ] {
        hsps.push(HspRecord {
            evalue,
            score,
            ..hsp((start, start + 10), (0, 10), 1, 10, &ops)
        });
    }
    let mut document = SamDocument::new(SamSpec::default(), "2.17.0+", None);
    document
        .query(
            &run_record(),
            &query_seq(0, "Query_1"),
            &[subject_seq(0, "Subject_1")],
            &hsps,
            &mut NoLetters,
        )
        .unwrap();
    let mut out = Vec::new();
    document.finish(&mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    let scores: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(11).unwrap())
        .collect();
    assert_eq!(
        scores,
        ["AS:i:15", "AS:i:13", "AS:i:12", "AS:i:11", "AS:i:14"]
    );
    // `SR` keeps the Seq-align order.
    let mut document = SamDocument::new(SamSpec::parse("SR"), "2.17.0+", None);
    document
        .query(
            &run_record(),
            &query_seq(0, "Query_1"),
            &[subject_seq(0, "Subject_1")],
            &hsps,
            &mut NoLetters,
        )
        .unwrap();
    let mut out = Vec::new();
    document.finish(&mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    let scores: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').nth(11).unwrap())
        .collect();
    assert_eq!(
        scores,
        ["AS:i:11", "AS:i:12", "AS:i:13", "AS:i:14", "AS:i:15"]
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:445-452
// (`AddSequence` drops a reference it has): one `@SQ` per reference of the whole run, in order of
// first appearance, whichever call adds it.
#[test]
fn the_sq_lines_list_each_reference_once_across_calls() {
    let ops = [GapEditOp::Sub(10)];
    let hsps = [
        hsp((0, 10), (0, 10), 1, 10, &ops),
        HspRecord {
            subject: 1,
            ..hsp((20, 30), (5, 15), 1, 10, &ops)
        },
    ];
    let subjects = [subject_seq(0, "Subject_1"), subject_seq(1, "Subject_2")];
    for (spec, expected) in [
        (
            SamSpec::default(),
            vec!["@SQ\tSN:Query_1\tLN:60", "@SQ\tSN:Query_2\tLN:60"],
        ),
        (
            SamSpec::parse("SR"),
            vec!["@SQ\tSN:Subject_1\tLN:100", "@SQ\tSN:Subject_2\tLN:100"],
        ),
    ] {
        let mut document = SamDocument::new(spec, "2.17.0+", None);
        for (index, id) in [(0, "Query_1"), (1, "Query_2"), (2, "Query_1")] {
            // the third call names the first query's reference key again
            let mut query = query_seq(index, id);
            query.key = if index == 2 { 0 } else { index as u64 };
            document
                .query(&run_record(), &query, &subjects, &hsps, &mut NoLetters)
                .unwrap();
        }
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let sq: Vec<&str> = text.lines().filter(|l| l.starts_with("@SQ")).collect();
        assert_eq!(sq, expected);
        assert_eq!(text.lines().filter(|l| !l.starts_with('@')).count(), 6);
        // a second finish writes nothing
        let mut again = Vec::new();
        document.finish(&mut again).unwrap();
        assert!(again.is_empty());
    }
}

#[test]
fn an_empty_document_writes_nothing_and_finish_empties_it() {
    let mut document = SamDocument::new(SamSpec::default(), "2.17.0+", None);
    let mut out = Vec::new();
    document.finish(&mut out).unwrap();
    assert!(out.is_empty());
    assert!(!document.has_records());
}
