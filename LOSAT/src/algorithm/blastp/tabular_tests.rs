use super::*;
use crate::blastinput::app::{FormatChoice, ReportFormat};
use crate::common::Hit;

fn hit(q: (usize, usize), s: (usize, usize), s_idx: u32) -> Hit {
    Hit {
        identity: 50.0,
        length: q.1 - q.0 + 1,
        mismatch: 2,
        gapopen: 0,
        q_start: q.0,
        q_end: q.1,
        s_start: s.0,
        s_end: s.1,
        e_value: 1.0e-10,
        bit_score: 40.5,
        num_ident: q.1 - q.0 - 1,
        query_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx,
        raw_score: 90,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 0,
    }
}

fn hit_list(hsp_lists: Vec<Vec<Hit>>, count: usize) -> BlastpHitList {
    BlastpHitList {
        hsplist_count: count,
        hsplist_max: hsp_lists.len(),
        worst_evalue: 0.0,
        low_score: 0,
        heapified: false,
        hsplist_array: hsp_lists
            .into_iter()
            .map(|hits| super::super::hsp::BlastpHspList {
                oid: hits[0].s_idx,
                query_index: 0,
                hsps: hits.into_iter().map(BlastpHsp::from_hit).collect(),
                best_evalue: 0.0,
            })
            .collect(),
        hsplist_current: 0,
        num_hits: 0,
    }
}

fn report(
    number: i32,
    fields: &str,
    app: BlastpApp,
    hit_lists: &[Option<BlastpHitList>],
) -> (String, String) {
    let format = ReportFormat::from_number(number).unwrap();
    let spec = TabularSpec::of_app(
        format,
        &FormatChoice {
            number,
            spec: fields.to_string(),
            delimiter: String::new(),
        },
        app.request_kind(),
    )
    .unwrap();
    let ids: Vec<crate::algorithm::common::tabular::SubjectId> =
        vec![Some(Arc::from(&b"s1"[..])), Some(Arc::from(&b"s2"[..]))];
    let length = |s_idx: usize| 200 + s_idx;
    let tab = BlastpTabular {
        spec: &spec,
        app,
        program_version: "PSIBLAST 2.17.0+",
        iteration: (app == BlastpApp::Psiblast).then_some(1),
        matrix: ScoringMatrix::Blosum62,
        query_ids: &[Arc::from(&b"q1"[..]), Arc::from(&b"q2"[..])],
        query_records: &[],
        query_titles: &[Arc::from(&b"q1 first"[..]), Arc::from(&b"q2"[..])],
        searched_lengths: &[100, 50],
        range_length: None,
        ungapped: false,
        database: "User specified sequence set (Input: s.faa)",
        subjects: TabularSubjects::Records {
            records: &[],
            ids: &ids,
            accessions: &[],
            parsed: &[],
            length: &length,
        },
        unsearched: &[false, true],
        hits_found: None,
        prolog: true,
        epilog: Some(2),
        taxdb_found: false,
        names: None,
        engine_counts: &crate::algorithm::blastp::seqalign::EngineCounts::default(),
        query_other_coding: None,
    };
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let before = [crate::report::diag_log::Posts::of_lines(
        b"Warning: q1\n".to_vec(),
    )];
    let mut warnings = QueryWarnings::of_posts(&before, None, &mut err);
    write_blastp_tabular(&tab, hit_lists, &mut out, Some(&mut warnings), None, None).unwrap();
    drop(warnings);
    (
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1061-1091 (frames 1, `sstrand` N/A of a protein subject), c++/src/algo/blast/format/blastfmtutil.cpp
// :629-640 (`qcovs`, the union of the subject's HSPs; `qcovus` is blastn's only), tabular.cpp
// :1264-1298 (outfmt 7's lines: psiblast's round, no count for an unsearched batch) and the
// taxonomy database's warning of each formatter (tabular.cpp:1093-1104).
#[test]
fn rows_and_comment_lines_are_ncbis() {
    let lists = [
        Some(hit_list(
            vec![
                vec![hit((1, 40), (1, 40), 0), hit((31, 60), (101, 130), 0)],
                vec![hit((5, 14), (5, 14), 1)],
            ],
            1,
        )),
        None,
    ];
    let (out, err) = report(
        6,
        "qseqid sseqid slen qcovs qcovhsp qcovus sstrand frames stitle",
        BlastpApp::Blastp,
        &lists,
    );
    assert_eq!(
        out,
        "q1\ts1\t200\t60\t40\tN/A\tN/A\t1/1\tN/A\nq1\ts1\t200\t60\t30\tN/A\tN/A\t1/1\tN/A\n"
    );
    assert_eq!(err, "Warning: q1\n");
    let (out, err) = report(7, "qseqid sseqid sscinames", BlastpApp::Psiblast, &lists);
    assert_eq!(
        out,
        "# PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: q1 first\n# Database: User specified sequence set (Input: s.faa)\n# Fields: query id, subject id, subject sci names\n# 2 hits found\nq1\ts1\tN/A\nq1\ts1\tN/A\n\
             # PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: q2\n# Database: User specified sequence set (Input: s.faa)\n# BLAST processed 2 queries\n"
    );
    let warning = taxdb_warning_line("psiblast");
    assert_eq!(err, format!("Warning: q1\n{warning}{warning}{warning}"));
    let (out, _) = report(20, "qseqid evalue", BlastpApp::Blastp, &lists);
    assert_eq!(out, "qseqid,evalue\nq1,1.00e-10\nq1,1.00e-10\n");
}

// NCBI 2.17.0 (line RG F/d36 probes bp.f20at, bp.f6at, ps.f20at, ps.f6at, ps.f7comma,
// rps.f10at, rps.f6at): blastp's reports take `delim=`, psiblast's and rpsblast's keep the
// format's delimiter, and their `6 delim=@` and `7 delim=,` are the default reports; the value
// is the command line's or the library's as given.
#[test]
fn the_spec_is_the_apps_formatters() {
    let spec = custom_tabular_spec("20 delim=@ qseqid evalue", BlastpApp::Blastp).unwrap();
    assert_eq!(
        (spec.kind, spec.delimiter.as_str()),
        (TabularKind::CsvWithHeader, "@")
    );
    assert_eq!(
        custom_tabular_spec("6 delim=@", BlastpApp::Blastp).map(|spec| spec.delimiter.as_str()),
        Some("@")
    );
    assert_eq!(
        custom_tabular_spec("10 delim=@ qseqid", BlastpApp::Blastp)
            .map(|spec| spec.delimiter.as_str()),
        Some("@")
    );
    for app in [BlastpApp::Psiblast, BlastpApp::Rpsblast] {
        let spec = custom_tabular_spec("20 delim=@ qseqid evalue", app).unwrap();
        assert_eq!(
            (spec.kind, spec.delimiter.as_str(), spec.fields.len()),
            (TabularKind::CsvWithHeader, ",", 2),
            "{app:?}"
        );
        let spec = custom_tabular_spec("10 delim=@ qseqid", app).unwrap();
        assert_eq!(spec.delimiter, ",", "{app:?}");
        assert!(custom_tabular_spec("6 delim=@", app).is_none(), "{app:?}");
        assert!(custom_tabular_spec("7 delim=,", app).is_none(), "{app:?}");
        assert!(custom_tabular_spec("6", app).is_none(), "{app:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:773-780 (an ungapped HSP's
// `align_length` is the length of its query range): an HSP of length 0 is one row whose
// `qend` is `qstart - 1`, with 0.000 identity and no coverage; `qcovs` leaves it out.
#[test]
fn an_hsp_of_length_zero_is_a_row_that_covers_nothing() {
    let mut empty = hit((3, 4), (24, 25), 0);
    empty.length = 0;
    empty.q_end = 2;
    empty.s_end = 23;
    empty.num_ident = 0;
    empty.mismatch = 0;
    empty.identity = 0.0;
    let lists = [
        Some(hit_list(vec![vec![hit((31, 60), (101, 130), 0), empty]], 1)),
        None,
    ];
    let (out, _) = report(
        6,
        "qseqid sseqid pident length qstart qend sstart send qcovs qcovhsp",
        BlastpApp::Blastp,
        &lists,
    );
    assert_eq!(
        out,
        "q1\ts1\t93.333\t30\t31\t60\t101\t130\t30\t30\nq1\ts1\t0.000\t0\t3\t2\t24\t23\t30\t0\n"
    );
}
