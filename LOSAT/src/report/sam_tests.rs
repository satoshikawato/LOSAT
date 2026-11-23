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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:67-82
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1149-1208
// (`CreateDensegFromDendiag`): an ungapped search's Seq-align of a subject holds a Dense-diag per
// HSP, which the formatter joins into one Dense-seg: one record per subject. Its tags are the
// first HSP's (the first of each score, sam_formatter.cpp:247-279), PI over the summed lengths;
// the sort key is the lesser of `GetSeqStart` and `GetSeqStop` of row 0 by the strand of the first
// segment (Dense_seg.cpp:165-238), ties broken by the last HSP's e-value (`s_GetBlastScore` keeps
// the last one, align_format_util.cpp:229-270). NCBI BLAST+ 2.17.0 sorts the same way (a subject
// whose lowest query start is lower sorts after one with a lower key; the tie above).
#[test]
fn an_ungapped_subject_is_one_record_of_its_joined_hsps() {
    let run = RunRecord {
        gapped: false,
        ..run_record()
    };
    let at = |subject: u32, q: (i32, i32), s: (i32, i32), frame: i8, evalue: f64, score: i32| {
        HspRecord {
            subject,
            evalue,
            score,
            ..hsp(q, s, frame, q.1 - q.0, &[])
        }
    };
    let hsps = [
        // Subject_1: key min(30, 20 + 5 - 1) = 24, though an HSP starts at 2.
        HspRecord {
            num_ident: 18,
            ..at(0, (30, 50), (0, 20), 1, 1e-9, 21)
        },
        at(0, (2, 12), (40, 50), 1, 1e-3, 22),
        at(0, (20, 25), (70, 75), 1, 1e-2, 23),
        // Subject_2: key min(10, 40 + 10 - 1) = 10, last e-value 1e-4.
        at(1, (10, 20), (5, 15), 1, 1e-10, 31),
        at(1, (40, 50), (30, 40), 1, 1e-4, 32),
        // Subject_3: the first HSP on the minus strand (query start 60 - 30 = 30), so the key
        // is min(40, 30 + 20 - 1) = 40 and the flag 16.
        at(2, (10, 30), (50, 70), -1, 1e-7, 41),
        at(2, (40, 48), (10, 18), 1, 1e-4, 42),
        // Subject_4: key 10, e-value 1e-8: before Subject_2, whose first e-value is lower.
        at(3, (10, 20), (60, 70), 1, 1e-8, 51),
    ];
    let subjects = [
        subject_seq(0, "Subject_1"),
        subject_seq(1, "Subject_2"),
        subject_seq(2, "Subject_3"),
        subject_seq(3, "Subject_4"),
    ];
    let records = |spec: SamSpec| {
        let mut document = SamDocument::new(spec, "2.17.0+", None);
        document
            .query(
                &run,
                &query_seq(0, "Query_1"),
                &subjects,
                &hsps,
                &mut NoLetters,
            )
            .unwrap();
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(str::to_string)
            .collect::<Vec<_>>()
    };
    let lines = records(SamSpec::default());
    let names: Vec<&str> = lines
        .iter()
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(names, ["Subject_4", "Subject_2", "Subject_1", "Subject_3"]);
    let fields: Vec<&str> = lines[2].split('\t').collect();
    // POS: the lowest query start + 1; the clips from the union of the subject rows [0, 75).
    assert_eq!(fields[1..6], ["0", "Query_1", "3", "255", "35M25H"]);
    assert_eq!(
        fields[11..],
        ["AS:i:21", "EV:f:1e-09", "NM:i:0", "PI:f:51.43", "BS:f:75.8"]
    );
    // Reversed: the subject rows [10, 70) clip 100 - 70 in front and 10 behind.
    let fields: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(fields[1..6], ["16", "Query_1", "31", "255", "30H28M10H"]);
    // `SR` keeps the Seq-align order: one record per subject.
    let lines = records(SamSpec::parse("SR"));
    let names: Vec<&str> = lines
        .iter()
        .map(|l| l.split('\t').nth(2).unwrap())
        .collect();
    assert_eq!(names, ["Subject_1", "Subject_2", "Subject_3", "Subject_4"]);
}

// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:267-276 and
// c++/src/objtools/format/cigar_formatter.cpp:164, :324-346: the aligned length that PI divides
// by is `len`, an `int`, from `GetAlignLength` with false (a `TSeqPos` sum) and the CIGAR count a
// `TSeqPos` sum. NCBI BLAST+ 2.17.0 on a 16-letter unit repeated 11,600 and 16,400 times against
// itself (`-ungapped -dust no`): `2152959968M` with `PI:f:-0.01`, and `8392672M`
// (4303359968 - 2^32) with `PI:f:3.13`.
#[test]
fn an_ungapped_record_sums_its_lengths_in_32_bits() {
    let run = RunRecord {
        gapped: false,
        ..run_record()
    };
    let unit = 1 << 30;
    let record = |copies: i32, ident: i32| {
        let hsps: Vec<HspRecord<'_>> = (0..copies)
            .map(|_| HspRecord {
                num_ident: ident,
                ..hsp((0, unit), (0, unit), 1, ident, &[])
            })
            .collect();
        let mut document = SamDocument::new(SamSpec::default(), "2.17.0+", None);
        document
            .query(
                &run,
                &query_seq(0, "Query_1"),
                &[subject_seq(0, "Subject_1")],
                &hsps,
                &mut NoLetters,
            )
            .unwrap();
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let line = text
            .lines()
            .find(|l| !l.starts_with('@'))
            .unwrap()
            .to_string();
        let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        (fields[5].clone(), fields[14].clone())
    };
    // 3 * 2^30 = 3221225472 is above 2^31 - 1: as an int it is -1073741824.
    assert_eq!(
        record(3, unit),
        ("3221225472M".to_string(), "PI:f:-100.00".to_string())
    );
    // 5 * 2^30 wraps to 2^30, which the identities then equal.
    assert_eq!(
        record(5, unit),
        ("1073741824M".to_string(), "PI:f:100.00".to_string())
    );
    // 4 * 2^30 wraps to 0: 100 * 1 / 0 is infinity and `min(99.99, inf)` 99.99.
    assert_eq!(record(4, 1), ("0M".to_string(), "PI:f:99.99".to_string()));
}

/// A TBLASTX HSP of subject `subject`: the residues `q` and `s` of the frames `frames`.
fn tbx_hsp(
    subject: u32,
    q: (i32, i32),
    s: (i32, i32),
    frames: (i8, i8),
    evalue: f64,
) -> HspRecord<'static> {
    HspRecord {
        subject,
        q_offset: q.0,
        q_end: q.1,
        s_offset: s.0,
        s_end: s.1,
        q_frame: frames.0,
        s_frame: frames.1,
        score: 40,
        num: 0,
        evalue,
        bit_score: 75.8,
        num_ident: 12,
        num_positives: 20,
        comp_adjust: 0,
        ops: &[],
    }
}

fn tbx_run() -> RunRecord<'static> {
    RunRecord {
        program: crate::blastinput::app::Program::Tblastx,
        rps: false,
        db_title: None,
        gapped: false,
    }
}

/// The document of one query (600 letters) with `hsps` over the subjects `Subject_1`, `Subject_2`
/// (700 letters each): the document and the number of Seq-aligns the formatter could not convert.
fn tbx_document(hsps: &[HspRecord<'_>]) -> (String, usize) {
    let query = SamSequence {
        record: SamId::Query(QueryRecord {
            id: SeqId::Local("Query_1"),
            searched_length: 600,
            offset: 0,
        }),
        key: 0,
        index: 0,
        length: 600,
    };
    let subject = |index: usize, id: &'static str| SamSequence {
        record: SamId::Subject(SubjectRecord {
            id: SeqId::Local(id),
            searched_length: 700,
            offset: 0,
            use_this_seqid: &[],
        }),
        key: index as u64,
        index,
        length: 700,
    };
    let subjects = [subject(0, "Subject_1"), subject(1, "Subject_2")];
    let mut document = SamDocument::new(SamSpec::default(), "2.17.0+", None)
        .for_program(crate::blastinput::app::Program::Tblastx);
    let failed = document
        .query(&tbx_run(), &query, &subjects, hsps, &mut NoLetters)
        .unwrap();
    let mut out = Vec::new();
    document.finish(&mut out).unwrap();
    (String::from_utf8(out).unwrap(), failed)
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:740-900 (`CreateDensegFromStdseg`)
// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:178-299 (`AddRow`)
// A TBLASTX subject with one HSP: nucleotide ranges (3 per residue) of the Std-seg, one record with
// the hard clips of the subject around its range, no score tags (the Seq-align holds none) and
// `@PG` names the options' program.
#[test]
fn a_translated_subject_with_one_hsp_has_one_record_without_score_tags() {
    let (text, failed) = tbx_document(&[tbx_hsp(0, (10, 40), (20, 50), (1, 1), 1e-9)]);
    assert_eq!(failed, 0);
    assert_eq!(
        text,
        "@HD\tVN:1.2\tSO:coordinate\tGO:reference\n\
         @SQ\tSN:Query_1\tLN:600\n\
         @PG\tID:0\tVN:2.17.0+\tPN:tblastx\n\
         Subject_1\t0\tQuery_1\t31\t255\t60H90M550H\t*\t0\t0\t*\t*\tNM:i:0\n"
    );
}

// The CIGAR's count is the sum of the segments' lengths, not the span: two HSPs of one strand pair
// are two segments of one M run, and the ranges are their unions.
#[test]
fn two_hsps_on_one_strand_add_up_in_one_m_count_over_the_union_of_their_ranges() {
    let (text, failed) = tbx_document(&[
        tbx_hsp(0, (10, 40), (20, 50), (1, 1), 1e-9),
        tbx_hsp(0, (100, 120), (200, 220), (1, 1), 1e-5),
    ]);
    assert_eq!(failed, 0);
    assert!(
        text.ends_with("Subject_1\t0\tQuery_1\t31\t255\t60H150M40H\t*\t0\t0\t*\t*\tNM:i:0\n"),
        "{text}"
    );
}

// A subject on the other strand of the subject: flag 16 and the clips swapped. The row's starts are
// those of the frame's interval (`ungapped_interval`).
#[test]
fn a_subject_on_the_minus_strand_is_a_reverse_record() {
    let (text, failed) = tbx_document(&[tbx_hsp(0, (10, 40), (20, 50), (1, -2), 1e-9)]);
    assert_eq!(failed, 0);
    assert!(
        text.ends_with("Subject_1\t16\tQuery_1\t31\t255\t61H90M549H\t*\t0\t0\t*\t*\tNM:i:0\n"),
        "{text}"
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:815-830
// Std-segs that put a row on different strands cannot be a Dense-seg: no record and no `@SQ` for
// that Seq-align, which `query` counts; the other subject is written; a document of failures only
// is empty (`Flush` writes nothing without a record).
#[test]
fn a_subject_whose_hsps_change_a_strand_has_no_record_and_is_counted() {
    let mixed = [
        tbx_hsp(0, (10, 40), (20, 50), (1, 1), 1e-9),
        tbx_hsp(0, (100, 120), (200, 220), (1, -1), 1e-5),
        tbx_hsp(1, (10, 40), (20, 50), (1, 1), 1e-9),
    ];
    let (text, failed) = tbx_document(&mixed);
    assert_eq!(failed, 1);
    assert!(text.contains("Subject_2\t0\tQuery_1\t31\t255\t"), "{text}");
    assert!(!text.contains("Subject_1"), "{text}");
    let (text, failed) = tbx_document(&mixed[..2]);
    assert_eq!(failed, 1);
    assert_eq!(text, "");
    // The query's strand changes too.
    let (text, failed) = tbx_document(&[
        tbx_hsp(0, (10, 40), (20, 50), (1, 1), 1e-9),
        tbx_hsp(0, (100, 120), (200, 220), (-1, 1), 1e-5),
    ]);
    assert_eq!((text.as_str(), failed), ("", 1));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:829-845
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1402-1430
// The sort reads the lowest query start over the Std-segs, and breaks a tie by the e-value of the
// FIRST Std-seg (the Seq-align has no scores of its own): Subject_1's first HSP has the greater
// e-value (its later one is smaller than Subject_2's), so Subject_2 comes first.
#[test]
fn a_translated_sort_ties_by_the_first_hsps_evalue() {
    let (text, failed) = tbx_document(&[
        tbx_hsp(0, (10, 40), (20, 50), (1, 1), 1e-5),
        tbx_hsp(0, (100, 120), (200, 220), (1, 1), 1e-30),
        tbx_hsp(1, (10, 40), (20, 50), (1, 1), 1e-10),
    ]);
    assert_eq!(failed, 0);
    let subjects: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap())
        .collect();
    assert_eq!(subjects, ["Subject_2", "Subject_1"]);
    // A lower start comes before a lower e-value.
    let (text, _) = tbx_document(&[
        tbx_hsp(0, (100, 120), (200, 220), (1, 1), 1e-30),
        tbx_hsp(1, (10, 40), (20, 50), (1, 1), 1.0),
    ]);
    let subjects: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').next().unwrap())
        .collect();
    assert_eq!(subjects, ["Subject_2", "Subject_1"]);
}
