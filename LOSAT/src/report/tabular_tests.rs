use super::*;

/// The spec of an application whose formatter gets the custom delimiter (blastp's).
fn blastp_spec(format: ReportFormat, choice: &FormatChoice) -> Option<TabularSpec> {
    TabularSpec::of_app(
        format,
        choice,
        crate::report::archive::options::RequestKind::Blastp,
    )
}

#[test]
fn fields_follow_ncbis_rules() {
    assert_eq!(parse_fields(""), STD_FIELDS.to_vec());
    assert_eq!(parse_fields("bogus"), STD_FIELDS.to_vec());
    assert_eq!(
        parse_fields("qseqid  sseqid qseqid -qseqid evalue"),
        vec![SubjectSeqId, Evalue]
    );
    let mut std_less = STD_FIELDS.to_vec();
    std_less.retain(|field| *field != Evalue);
    std_less.push(QueryCovHsp);
    assert_eq!(parse_fields("std -evalue qcovhsp"), std_less);
    for (index, (name, field, _)) in FIELDS.iter().enumerate() {
        assert_eq!(*field as usize, index);
        assert_eq!(TabularField::from_name(name), Some(*field));
    }
}

// Line DBID P1: with NCBI's prefetch off, the rows read the subject's Bioseq (and look its
// taxids up) only for these fields (review 1 of P1, field by field against NCBI 2.17.0); P2
// (review 1): in the translated programs also for the subject's frame with a coordinate; L1
// (review 2 of P2): `nident` without another field of the block only on rows with identities.
#[test]
fn rows_read_the_subject_for_its_ids_taxonomy_titles_strand_and_letters() {
    use crate::blastinput::app::Program;
    let spec = |fields: &str| TabularSpec {
        kind: TabularKind::Tabular,
        fields: parse_fields(fields),
        delimiter: "\t".into(),
    };
    let fetch = |fields: &str, no_fetch: bool| {
        spec(fields).subject_fetch(if no_fetch {
            Program::Blastn
        } else {
            Program::Blastp
        })
    };
    for (fields, translated) in [
        ("qseqid sframe", SubjectFetch::Never),
        ("qseqid frames qcovs", SubjectFetch::Never),
        ("qseqid sstart sframe", SubjectFetch::Every),
        ("qseqid qframe frames", SubjectFetch::Every),
        ("qseqid length sframe", SubjectFetch::Every),
        ("qseqid qstart qframe slen", SubjectFetch::Never),
        ("qseqid nident sframe", SubjectFetch::WithIdentities),
        ("qseqid nident frames sstart", SubjectFetch::Every),
    ] {
        for program in [Program::Blastx, Program::Tblastn, Program::Tblastx] {
            assert_eq!(spec(fields).subject_fetch(program), translated, "{fields}");
        }
        let protein = match translated {
            SubjectFetch::Every if !fields.contains("nident") => SubjectFetch::Never,
            fetch => fetch,
        };
        assert_eq!(
            spec(fields).subject_fetch(Program::Blastp),
            protein,
            "{fields}"
        );
        assert_eq!(
            spec(fields).subject_fetch(Program::Blastn),
            SubjectFetch::Never,
            "{fields}"
        );
    }
    for fields in [
        "qseqid evalue",
        "qseqid qacc qlen slen qstart qend sstart send evalue bitscore score length gaps gapopen",
        "qseqid qframe sframe frames qcovs qcovhsp qcovus",
    ] {
        assert_eq!(fetch(fields, false), SubjectFetch::Never, "{fields}");
    }
    for fields in [
        "std",
        "qseqid sseqid",
        "qseqid sallacc",
        "qseqid staxid",
        "qseqid sskingdoms",
        "qseqid salltitles",
        "qseqid sstrand",
        "qseqid sseq",
        "qseqid btop",
        "qseqid ppos",
        "qseqid nident sseqid",
    ] {
        assert_eq!(fetch(fields, false), SubjectFetch::Every, "{fields}");
        assert_eq!(fetch(fields, true), SubjectFetch::Every, "{fields}");
    }
    for fields in ["qseqid mismatch", "qseqid pident", "qseqid nident length"] {
        assert_eq!(fetch(fields, false), SubjectFetch::Every, "{fields}");
        assert_eq!(fetch(fields, true), SubjectFetch::Never, "{fields}");
    }
    for fields in ["qseqid nident", "nident evalue qcovs slen"] {
        assert_eq!(
            fetch(fields, false),
            SubjectFetch::WithIdentities,
            "{fields}"
        );
        assert_eq!(fetch(fields, true), SubjectFetch::Never, "{fields}");
    }
    assert!(!SubjectFetch::WithIdentities.reads(0) && SubjectFetch::WithIdentities.reads(1));
    assert!(SubjectFetch::Every.reads(0) && !SubjectFetch::Never.reads(1));
}

// Line DBID tabfields: a table reads the sequence that NCBI finds by the subject's Seq-align
// id where its rows read the Bioseq or letters, or its length (`slen`); NCBI 2.17.0 writes
// every other list the same for a subject it finds as another sequence or not at all (probe
// `blastdb-ids-20261021/R/tabfields`).
#[test]
fn tables_read_the_found_sequence_for_the_bioseq_its_letters_and_slen() {
    use crate::blastinput::app::Program;
    let spec = |fields: &str| TabularSpec {
        kind: TabularKind::Tabular,
        fields: parse_fields(fields),
        delimiter: "\t".into(),
    };
    let all = [
        Program::Blastp,
        Program::Blastn,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
    ];
    assert!(all.iter().all(|&program| reads_subject(None, program)));
    for fields in [
        "qseqid evalue",
        "qseqid qlen qstart qend sstart send evalue bitscore score length gaps gapopen qcovhsp",
        "qseqid qcovs qcovus",
    ] {
        for program in all {
            assert!(!spec(fields).reads_subject(program), "{fields} {program:?}");
            assert!(!reads_subject(Some(&spec(fields)), program), "{fields}");
        }
    }
    for fields in [
        "std",
        "qseqid sseqid",
        "qseqid slen",
        "evalue stitle",
        "qseqid sseq",
    ] {
        for program in all {
            assert!(spec(fields).reads_subject(program), "{fields} {program:?}");
        }
    }
    // BLASTN's formatter counts the identities without the letters; the translated
    // programs read the subject for its frame.
    for fields in ["qseqid pident", "qseqid nident", "qseqid mismatch length"] {
        assert!(!spec(fields).reads_subject(Program::Blastn), "{fields}");
        assert!(spec(fields).reads_subject(Program::Blastp), "{fields}");
    }
    assert!(spec("qseqid sstart sframe").reads_subject(Program::Tblastn));
    assert!(!spec("qseqid sstart sframe").reads_subject(Program::Blastp));
}

#[test]
fn delimiter_conflicts_are_ncbis() {
    assert!(custom_delimiter_conflict("qseqid staxids", ";")
        .unwrap()
        .contains("(staxids)"));
    assert!(custom_delimiter_conflict("std -frames", "/").is_some());
    assert!(custom_delimiter_conflict("salltitles", "<>").is_some());
    assert!(custom_delimiter_conflict("qseqid", ";").is_none());
    assert!(custom_delimiter_conflict("  ", ";").is_none());
}

#[test]
fn default_reports_have_no_spec() {
    let choice = |spec: &str, delimiter: &str| FormatChoice {
        number: 6,
        spec: spec.to_string(),
        delimiter: delimiter.to_string(),
    };
    assert!(blastp_spec(ReportFormat::Tabular, &choice("", "")).is_none());
    assert!(blastp_spec(ReportFormat::Tabular, &choice("std", "\t")).is_none());
    assert!(blastp_spec(ReportFormat::TabularWithComments, &choice("", "")).is_none());
    let csv = blastp_spec(ReportFormat::CommaSeparatedValues, &choice("", "")).unwrap();
    assert_eq!((csv.delimiter.as_str(), csv.fields.len()), (",", 12));
    let piped = blastp_spec(ReportFormat::Tabular, &choice("std", "|")).unwrap();
    assert_eq!(piped.delimiter, "|");
}

// NCBI 2.17.0 (line RG F/d36 probes ps.f10at, ps.f20at, ps.f6at, ps.pssm.f10at, bp.f20at,
// bp.f6at, rps.f10at, rps.f6at, rpst.f10at): `-outfmt "10 delim=@ qseqid sseqid evalue"`
// writes `,` in psiblast, rpsblast and rpstblastn, `@` in blastp; `6 delim=@` is the default
// outfmt 6 of the first three.
#[test]
fn only_the_apps_that_set_the_delimiter_take_it() {
    use crate::report::archive::options::RequestKind;
    let choice = |number: i32, spec: &str| FormatChoice {
        number,
        spec: spec.to_string(),
        delimiter: "@".to_string(),
    };
    let csv = choice(10, "qseqid sseqid evalue");
    let std6 = choice(6, "");
    let header = choice(20, "qseqid");
    for app in [
        RequestKind::Blastn,
        RequestKind::Blastp,
        RequestKind::Blastx,
        RequestKind::Tblastn,
        RequestKind::PsiTblastn,
        RequestKind::Tblastx,
    ] {
        let spec = TabularSpec::of_app(ReportFormat::CommaSeparatedValues, &csv, app).unwrap();
        assert_eq!(
            (spec.kind, spec.delimiter.as_str(), spec.fields.len()),
            (TabularKind::Csv, "@", 3),
            "{app:?}"
        );
        let spec = TabularSpec::of_app(ReportFormat::Tabular, &std6, app).unwrap();
        assert_eq!(spec.delimiter, "@", "{app:?}");
    }
    for app in [
        RequestKind::Psiblast,
        RequestKind::Rpsblast,
        RequestKind::Rpstblastn,
    ] {
        let spec = TabularSpec::of_app(ReportFormat::CommaSeparatedValues, &csv, app).unwrap();
        assert_eq!(
            (spec.kind, spec.delimiter.as_str(), spec.fields.len()),
            (TabularKind::Csv, ",", 3),
            "{app:?}"
        );
        let spec = TabularSpec::of_app(ReportFormat::CommaSeparatedValuesWithHeader, &header, app);
        assert_eq!(spec.unwrap().delimiter, ",", "{app:?}");
        assert!(TabularSpec::of_app(ReportFormat::Tabular, &std6, app).is_none());
        let commented = choice(7, "qseqid");
        let spec = TabularSpec::of_app(ReportFormat::TabularWithComments, &commented, app);
        assert_eq!(spec.unwrap().delimiter, "\t", "{app:?}");
    }
}

#[test]
fn coverage_merges_touching_ranges() {
    assert_eq!(union_length(&mut [(10, 19), (0, 9), (30, 39)]), 30);
    assert_eq!(union_length(&mut [(0, 9), (5, 14)]), 15);
    assert_eq!(coverage_percent(0, 100), -1);
    assert_eq!(coverage_percent(50, 100), 50);
    assert_eq!(coverage_percent(995, 1000), 99);
    assert_eq!(coverage_percent(1000, 1000), 100);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1397-1447
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3162-3205
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:737-747
// The values that NCBI BLAST+ 2.17.0 writes. BLASTN: x1bp fix 2's d3mK probe, a 200-letter HSP
// with -query_loc 1203-1000000 of a record of 2202 letters (998798 as given, 1000 searched).
// BLASTP and TBLASTN -ungapped -comp_based_stats 0: rows db.x1bp.cov.bp.6 and db.x1bp.cov.tbn.6
// (q_plain, 142 residues, with -query_loc 1-1000000), the review's probe D without the range and
// gapped. BLASTX: tab.bx.cov.6 (qn2, 135 letters; a gapped search counts a frame's residues
// against a third of them), with -query_loc a range of 300 as given.
/// The `qcovhsp` of one HSP: its subject's rule (`QueryCoveragePerHsp::of_subject`), then the
/// row's division.
fn query_coverage_per_hsp(
    ungapped: bool,
    hsps: usize,
    span: i32,
    context_length: i32,
    seqalign_span: i32,
    query_length: i32,
) -> i32 {
    QueryCoveragePerHsp::of_subject(ungapped, hsps, context_length, query_length)
        .of_hsp(span, seqalign_span)
}

#[test]
fn query_coverage_per_hsp_reads_the_seq_align_that_the_formatter_gets() {
    // A subject's single ungapped HSP counts against the range as given; several, or a gapped
    // HSP, against the searched letters.
    assert_eq!(query_coverage_per_hsp(true, 1, 200, 1000, 200, 998_798), 0);
    assert_eq!(query_coverage_per_hsp(true, 2, 200, 1000, 200, 998_798), 20);
    assert_eq!(
        query_coverage_per_hsp(false, 1, 200, 1000, 200, 998_798),
        20
    );
    // BLASTP: 120 residues; TBLASTN: 142.
    assert_eq!(query_coverage_per_hsp(true, 1, 120, 142, 120, 1_000_000), 0);
    assert_eq!(query_coverage_per_hsp(true, 1, 120, 142, 120, 142), 85);
    assert_eq!(
        query_coverage_per_hsp(false, 1, 120, 142, 120, 1_000_000),
        85
    );
    assert_eq!(query_coverage_per_hsp(true, 1, 142, 142, 142, 1_000_000), 0);
    assert_eq!(query_coverage_per_hsp(true, 1, 142, 142, 142, 142), 100);
    // BLASTX: residues of the frame in the HSP, the Seq-align's nucleotides.
    assert_eq!(query_coverage_per_hsp(false, 2, 6, 45, 18, 135), 13);
    assert_eq!(query_coverage_per_hsp(false, 1, 15, 45, 45, 135), 33);
    assert_eq!(query_coverage_per_hsp(true, 2, 6, 135, 18, 135), 4);
    assert_eq!(query_coverage_per_hsp(true, 1, 6, 135, 18, 135), 13);
    assert_eq!(query_coverage_per_hsp(false, 1, 140, 140, 420, 420), 100);
    assert_eq!(query_coverage_per_hsp(false, 1, 15, 45, 45, 300), 33);
    assert_eq!(query_coverage_per_hsp(true, 2, 6, 135, 18, 300), 4);
    assert_eq!(query_coverage_per_hsp(true, 1, 6, 135, 18, 300), 6);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:626-633 and 693-700: without the taxonomy database's names, `staxids`
// brings `N/A` fillers, quoted in a CSV (the singular names only for a non-zero taxid, the
// plural ones once per taxid).
#[test]
fn taxonomy_names_get_ncbis_fillers() {
    let row = |format: ReportFormat,
               fields: &str,
               taxonomy: &TabularTaxonomy,
               names: Option<&TabularNames>| {
        let choice = FormatChoice {
            number: 10,
            spec: fields.to_string(),
            delimiter: String::new(),
        };
        let spec = blastp_spec(format, &choice).unwrap();
        let writer = TabularWriter::new(&spec, "BLASTN 2.17.0+", None);
        let subject = TabularSubject {
            id: b"s",
            accession: b"s",
            accession_version: b"s",
            gi: 0,
            all_ids: &[],
            length: 10,
            titles: None,
            taxonomy,
            names,
            nucleotide: true,
            qcovs: -1,
            qcovus: -1,
            lookup: SubjectLookup::Found,
            index: 0,
        };
        let hsp = TabularHsp {
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            evalue: 1e-5,
            bit_score: 20.0,
            score: 10,
            length: 10,
            gap_opens: 0,
            gaps: 0,
            num_ident: 10,
            q_frame: 1,
            s_frame: 1,
            rows: None,
            ops: &[],
            qcovhsp: 100,
        };
        let mut out = Vec::new();
        writer.row(
            &mut out,
            &mut TabularRows::default(),
            &TabularQuery::local(b"q", 10),
            &subject,
            None,
            &hsp,
        );
        String::from_utf8(out).unwrap()
    };
    let two = TabularTaxonomy {
        taxid: 9606,
        taxids: vec![9606, 10090],
    };
    let fields =
        "ssciname scomname sblastname sskingdom sscinames scomnames sblastnames sskingdoms";
    assert_eq!(
        row(
            ReportFormat::CommaSeparatedValues,
            &format!("{fields} staxids"),
            &two,
            None
        ),
        "\"N/A\",\"N/A\",N/A,N/A,\"N/A;N/A\",\"N/A;N/A\",\"N/A\",N/A,9606;10090\n"
    );
    assert_eq!(
        row(
            ReportFormat::Tabular,
            &format!("{fields} staxids"),
            &two,
            None
        ),
        "N/A\tN/A\tN/A\tN/A\tN/A;N/A\tN/A;N/A\tN/A\tN/A\t9606;10090\n"
    );
    assert_eq!(
        row(ReportFormat::CommaSeparatedValues, fields, &two, None),
        "N/A,N/A,N/A,N/A,\"N/A\",\"N/A\",\"N/A\",N/A\n"
    );
    let none = TabularTaxonomy {
        taxid: 0,
        taxids: vec![0],
    };
    assert_eq!(
        row(
            ReportFormat::CommaSeparatedValues,
            "ssciname sscinames staxid staxids",
            &none,
            None
        ),
        "N/A,\"N/A\",N/A,0\n"
    );
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:612-624 and 669-680 with the names found (RESULTS.md §3 of line OF
    // taxdb): the plural blast names and super kingdoms are sorted sets, `-` and
    // `unclassified` are dropped, an empty common name stays empty in the plural list and is
    // `N/A` alone; a taxid not found gets the fillers.
    let entry = |taxid, sci: &str, common: &str, blast: &str, kingdom: &str| {
        let names = crate::seqdb::taxdb::TaxNames {
            taxid,
            scientific_name: sci.into(),
            common_name: common.into(),
            blast_name: blast.into(),
            super_kingdom: Ok(kingdom.into()),
        };
        (taxid, (Some(names), Vec::new()))
    };
    let found = TabularNames {
        program: "blastp".into(),
        open: Vec::new(),
        opened: std::cell::Cell::new(false),
        answers: [
            entry(9606, "Homo sapiens", "human", "primates", "Eukaryota"),
            entry(7777, "Quote \"x\", y; z", "", "-", "Eukaryota"),
            entry(32630, "synthetic construct", "", "unclassified", "-"),
        ]
        .into_iter()
        .collect(),
        report: None,
    };
    let three = TabularTaxonomy {
        taxid: 9606,
        taxids: vec![7777, 9606, 10090, 32630],
    };
    assert_eq!(
        row(
            ReportFormat::CommaSeparatedValues,
            &format!("{fields} staxids"),
            &three,
            Some(&found)
        ),
        "\"Homo sapiens\",\"human\",\"primates\",Eukaryota,\"Quote \\\"x\\\", y; z;Homo sapiens;N/A;synthetic construct\",\";human;N/A;\",\"primates\",Eukaryota,7777;9606;10090;32630\n"
    );
    let quote = TabularTaxonomy {
        taxid: 7777,
        taxids: vec![7777],
    };
    assert_eq!(
        row(ReportFormat::Tabular, fields, &quote, Some(&found)),
        "Quote \"x\", y; z\tN/A\tN/A\tEukaryota\tQuote \"x\", y; z\t\tN/A\tEukaryota\n"
    );
}

// NCBI writes the texts of opening at the first lookup of the report, and a lookup's own
// texts at each lookup: the plural names' lookups before the singular one's.
#[test]
fn row_messages_follow_the_lookups() {
    let spec = blastp_spec(
        ReportFormat::Tabular,
        &FormatChoice {
            number: 6,
            spec: "ssciname sscinames".to_string(),
            delimiter: String::new(),
        },
    )
    .unwrap();
    let names = TabularNames {
        program: "blastn".into(),
        open: vec!["open".into()],
        opened: std::cell::Cell::new(false),
        answers: [
            (99990, (None, vec!["offset".into()])),
            (9606, (None, Vec::new())),
        ]
        .into_iter()
        .collect(),
        report: None,
    };
    let none = TabularTaxonomy::default();
    assert!(names.row_messages(&spec, &none).is_empty());
    let last = TabularTaxonomy {
        taxid: 99990,
        taxids: vec![9606, 99990],
    };
    assert_eq!(
        String::from_utf8(names.row_messages(&spec, &last)).unwrap(),
        "Error: [blastn] open\nError: [blastn] offset\nError: [blastn] offset\n"
    );
    assert_eq!(
        String::from_utf8(names.row_messages(&spec, &last)).unwrap(),
        "Error: [blastn] offset\nError: [blastn] offset\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908 (the frames only with a field of `SetCounts`' block, `num_ident >
// 0` for `nident`), tabular.cpp:1003-1021 (the positives of the rows by the formatter's
// matrix, `GetAsciiProteinMatrix`: U and C count as C-C, a J pair not, `*` with `*` is 1)
// and tabular.cpp:1293-1298 (psiblast's round).
#[test]
fn protein_rows_frames_and_rounds_are_ncbis() {
    let spec = |format: ReportFormat, fields: &str| {
        blastp_spec(
            format,
            &FormatChoice {
                number: 6,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap()
    };
    let taxonomy = TabularTaxonomy::default();
    let subject = TabularSubject {
        id: b"s",
        accession: b"s",
        accession_version: b"s",
        gi: 0,
        all_ids: &[],
        length: 10,
        titles: None,
        taxonomy: &taxonomy,
        names: None,
        nucleotide: false,
        qcovs: -1,
        qcovus: -1,
        lookup: SubjectLookup::Found,
        index: 0,
    };
    let row = |spec: &TabularSpec, num_ident: usize, rows: Option<(&[u8], &[u8])>| {
        let writer = TabularWriter::new(spec, "BLASTP 2.17.0+", Some(ScoringMatrix::Blosum62));
        let hsp = TabularHsp {
            q_start: 1,
            q_end: 6,
            s_start: 1,
            s_end: 6,
            evalue: 1e-5,
            bit_score: 20.0,
            score: 10,
            length: 6,
            gap_opens: 0,
            gaps: 0,
            num_ident,
            q_frame: 1,
            s_frame: 1,
            rows,
            ops: &[],
            qcovhsp: 100,
        };
        let mut out = Vec::new();
        writer.row(
            &mut out,
            &mut TabularRows::default(),
            &TabularQuery::local(b"q", 6),
            &subject,
            None,
            &hsp,
        );
        String::from_utf8(out).unwrap()
    };
    let frames = spec(ReportFormat::Tabular, "qseqid frames sframe");
    assert_eq!(row(&frames, 3, None), "q\t0/0\t0\n");
    let counted = spec(ReportFormat::Tabular, "nident frames sframe");
    assert_eq!(row(&counted, 0, None), "0\t0/0\t0\n");
    assert_eq!(row(&counted, 3, None), "3\t1/1\t1\n");
    let positives = spec(ReportFormat::Tabular, "nident positive mismatch sstrand");
    assert_eq!(
        row(&positives, 1, Some((b"UCJ*KA", b"CUJ*IA"))),
        "3\t5\t3\tN/A\n"
    );
    assert_eq!(row(&positives, 1, Some((b"UJ", b"UI"))), "1\t1\t5\tN/A\n");
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:444-449 and
    // :482-485 (`ppos` and `positive` print the count of the positives, cited at
    // `TabularWriter::new`): either field alone keeps the matrix (U-C, C-U, J-J, *-*, A-A
    // of 6; K-I is -3), so a writer that kept it for one of them only would print 50.00 / 3.
    let percent = spec(ReportFormat::Tabular, "ppos");
    assert_eq!(row(&percent, 1, Some((b"UCJ*KA", b"CUJ*IA"))), "83.33\n");
    let count = spec(ReportFormat::Tabular, "positive");
    assert_eq!(row(&count, 1, Some((b"UCJ*KA", b"CUJ*IA"))), "5\n");
    let commented = spec(ReportFormat::TabularWithComments, "qseqid");
    let writer = TabularWriter::new(
        &commented,
        "PSIBLAST 2.17.0+",
        Some(ScoringMatrix::Blosum62),
    );
    let mut out = Vec::new();
    writer
        .round_query_header(&mut out, Some(1), b"q1 title", b"db", Some(0))
        .unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "# PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: q1 title\n# Database: db\n# 0 hits found\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:219-255 (`sseqid`, `sacc`, `saccver`, `sgi`: the subject's own ids) and
// :387-439 (`salltitles` joined by `<>`, `stitle` the first title, an empty one `N/A`); a local
// query's id is its `qseqid`, `qacc` and `qaccver` alike.
#[test]
fn id_columns_are_the_subject_s_own() {
    let spec = blastp_spec(
        ReportFormat::Tabular,
        &FormatChoice {
            number: 6,
            spec: "qseqid qacc qaccver sseqid sacc saccver sgi stitle salltitles".to_string(),
            delimiter: String::new(),
        },
    )
    .unwrap();
    let taxonomy = TabularTaxonomy::default();
    let titles: [&[u8]; 2] = [b"first title", b""];
    let subject = TabularSubject {
        id: b"gb|AE000111.1|",
        accession: b"AE000111",
        accession_version: b"AE000111.1",
        gi: 42,
        all_ids: &[],
        length: 10,
        titles: Some(&titles),
        taxonomy: &taxonomy,
        names: None,
        nucleotide: false,
        qcovs: -1,
        qcovus: -1,
        lookup: SubjectLookup::Found,
        index: 0,
    };
    let hsp = TabularHsp {
        q_start: 1,
        q_end: 6,
        s_start: 1,
        s_end: 6,
        evalue: 1e-5,
        bit_score: 20.0,
        score: 10,
        length: 6,
        gap_opens: 0,
        gaps: 0,
        num_ident: 6,
        q_frame: 1,
        s_frame: 1,
        rows: None,
        ops: &[],
        qcovhsp: 100,
    };
    let mut out = Vec::new();
    TabularWriter::new(&spec, "BLASTP 2.17.0+", None).row(
        &mut out,
        &mut TabularRows::default(),
        &TabularQuery::local(b"q1", 6),
        &subject,
        None,
        &hsp,
    );
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "q1\tq1\tq1\tgb|AE000111.1|\tAE000111\tAE000111.1\t42\tfirst title\tfirst title<>N/A\n"
    );
}

/// The rows of `subjects` (failed or found, in order) of one query under `fields`, with a
/// nucleotide HSP of 450 letters (qcovs 99, qcovhsp 50) and query length 900.
fn failed_rows(fields: &str, subjects: &[TabularSubject<'_>]) -> String {
    let spec = blastp_spec(
        ReportFormat::Tabular,
        &FormatChoice {
            number: 6,
            spec: fields.to_string(),
            delimiter: String::new(),
        },
    )
    .unwrap();
    let hsp = TabularHsp {
        q_start: 1,
        q_end: 450,
        s_start: 1,
        s_end: 450,
        evalue: 0.0,
        bit_score: 812.0,
        score: 900,
        length: 450,
        gap_opens: 0,
        gaps: 0,
        num_ident: 450,
        q_frame: 1,
        s_frame: 1,
        rows: None,
        ops: &[],
        qcovhsp: 50,
    };
    let writer = TabularWriter::new(&spec, "BLASTN 2.17.0+", None);
    let mut rows = TabularRows::default();
    let mut out = Vec::new();
    for subject in subjects {
        // The kept subject of a failed row: the last found one, by its index.
        let kept = rows
            .last_found()
            .and_then(|index| subjects.iter().find(|subject| subject.index == index));
        writer.row(
            &mut out,
            &mut rows,
            &TabularQuery::local(b"nq1", 900),
            subject,
            kept.filter(|_| writer.fails(subject)),
            &hsp,
        );
    }
    String::from_utf8(out).unwrap()
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:828-875 and 1374-1378 (the catch of `SetFields` and `SetSubjectId` of `list`):
// NCBI BLAST+ 2.17.0 `blastn -subject` with two titles `Subject_9 bad x\x81y title`,
// `Subject_7 other \x8dz` (probe `remaining-gaps-20261112/T/t7`, `s.bb.blastn.f6wide`), and a
// found subject before a failed one (`s.gb`): `sseqid` Unknown, `sgi` 0, the Seq-align ids added
// to `sall*` row by row, the counts of `x_ResetFields`, `slen` 0, the scores, the query's ids
// and coverages real, the titles and taxonomy of the last found row.
#[test]
fn failed_lookups_write_ncbis_catch_row() {
    let taxonomy = TabularTaxonomy::default();
    let subject = |id: &'static [u8], failed: Option<&'static [u8]>, index| TabularSubject {
        id,
        accession: id,
        accession_version: id,
        gi: 0,
        all_ids: &[],
        length: 1350,
        titles: None,
        taxonomy: &taxonomy,
        names: None,
        nucleotide: true,
        qcovs: 99,
        qcovus: -1,
        lookup: match failed {
            Some(align_id) => SubjectLookup::Failed(FailedLookup {
                full: align_id,
                accession: align_id,
                subject_id: true,
                all_ids: None,
                empty_handle: false,
            }),
            None => SubjectLookup::Found,
        },
        index,
    };
    let wide = "qseqid sseqid evalue bitscore pident length mismatch qstart sstart slen qlen \
                    nident sgi sallseqid sallgi sallacc stitle salltitles staxids sscinames \
                    sstrand qcovs qcovhsp";
    let bad = |n: &'static [u8]| subject(b"", Some(n), 9);
    let fallback = "\t0.0\t812\t0.000\t0\t0\t0\t0\t0\t900\t0\t0\t";
    let tail = "\tN/A\tN/A\tN/A\tN/A\tN/A\t99\t50\n";
    assert_eq!(
        failed_rows(
            wide,
            &[bad(b"Subject_1"), bad(b"Subject_1"), bad(b"Subject_2")]
        ),
        format!(
            "nq1\tUnknown{fallback}Subject_1\t0\tSubject_1{tail}\
                 nq1\tUnknown{fallback}Subject_1;Subject_1\t0;0\tSubject_1;Subject_1{tail}\
                 nq1\tUnknown{fallback}Subject_1;Subject_1;Subject_2\t0;0;0\t\
                 Subject_1;Subject_1;Subject_2{tail}"
        )
    );
    let rows = failed_rows(wide, &[subject(b"s2", None, 1), bad(b"Subject_2")]);
    assert_eq!(
        rows.lines().nth(1).unwrap(),
        format!(
            "nq1\tUnknown{fallback}s2;Subject_2\t0;0\ts2;Subject_2{}",
            tail.trim_end()
        )
    );
    // Without `sseqid`, `sgi`, `sacc` and `saccver` a `-subject` record's lookup does not run
    // (`x_SetSubjectIds` of a record has no `GetSeqIdList`): the row is the found one's.
    assert_eq!(
        failed_rows(
            "qseqid sallseqid evalue",
            &[subject(b"Subject_9", Some(b"Subject_1"), 0)]
        ),
        "nq1\tSubject_9\t0.0\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:828-875 (the catch of the subject block, `SetSubjectId` of `list`, then a return of -1)
// where the load of the subject throws (`DB contains no accession info.`, a batch whose
// prefetch is off; probe of line DBID pfoff-2, `R/pfoff/p2`, `bp.U2.*`, DECISIONS 62 (12)):
// the row keeps the last loaded row's ids, titles and taxonomy, `sall*` add the Seq-align id
// (its GI in `sallgi`), `slen` and the counts are `x_ResetFields`' 0, `qlen` comes only with
// a field of the query block; before any loaded row `sseqid` is `Unknown`, the titles and
// names `N/A` (quoted in CSV).
#[test]
fn a_load_that_throws_writes_the_kept_row() {
    let exception =
        crate::report::db_subject::LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO)
            .exception()
            .unwrap();
    let taxonomy = TabularTaxonomy {
        taxid: 9606,
        taxids: vec![9606],
    };
    let none = TabularTaxonomy::default();
    let titles: [&[u8]; 1] = [b"found title"];
    let all_ids = [DeflineIds {
        id: b"gb|AAA00001.1|".to_vec(),
        gi: 7,
        accession: b"AAA00001".to_vec(),
    }];
    let found = TabularSubject {
        id: b"gb|AAA00001.1|",
        accession: b"AAA00001",
        accession_version: b"AAA00001.1",
        gi: 7,
        all_ids: &all_ids,
        length: 1350,
        titles: Some(&titles),
        taxonomy: &taxonomy,
        names: None,
        nucleotide: true,
        qcovs: 99,
        qcovus: -1,
        lookup: SubjectLookup::Found,
        index: 0,
    };
    let thrown = TabularSubject {
        id: b"",
        accession: b"",
        accession_version: b"",
        gi: 0,
        all_ids: &[],
        length: 0,
        titles: None,
        taxonomy: &none,
        names: None,
        nucleotide: true,
        qcovs: 98,
        qcovus: -1,
        lookup: SubjectLookup::Throws(ThrownLookup {
            full: b"gi|5|gb|BBB00002.1|",
            gi: 5,
            accession: b"BBB00002",
            exception: &exception,
        }),
        index: 1,
    };
    assert!(thrown.throws() && !found.throws());
    assert_eq!(
        failed_rows(
            "qseqid sseqid sacc saccver sgi slen qlen length pident evalue stitle staxids \
                 sstrand qcovs",
            &[found, thrown]
        )
        .lines()
        .nth(1)
        .unwrap(),
        "nq1\tgb|AAA00001.1|\tAAA00001\tAAA00001.1\t7\t0\t900\t0\t0.000\t0.0\tfound title\t9606\tN/A\t98"
    );
    // `qlen` without a field of the query block keeps `x_ResetFields`' 0.
    assert_eq!(
        failed_rows("sseqid qlen slen evalue", &[found, thrown])
            .lines()
            .nth(1)
            .unwrap(),
        "gb|AAA00001.1|\t0\t0\t0.0"
    );
    // `sall*`: the kept lists and the Seq-align id, row by row.
    assert_eq!(
        failed_rows("qseqid sallseqid sallgi sallacc", &[found, thrown, thrown]),
        "nq1\tgb|AAA00001.1|\t7\tAAA00001\n\
             nq1\tgb|AAA00001.1|;gi|5|gb|BBB00002.1|\t7;5\tAAA00001;BBB00002\n\
             nq1\tgb|AAA00001.1|;gi|5|gb|BBB00002.1|;gi|5|gb|BBB00002.1|\t7;5;5\t\
             AAA00001;BBB00002;BBB00002\n"
    );
    // Before any loaded row: `Unknown`, `sgi` 0, the titles `N/A`.
    assert_eq!(
        failed_rows("qseqid sseqid sgi stitle sallseqid evalue", &[thrown]),
        "nq1\tUnknown\t0\tN/A\tgi|5|gb|BBB00002.1|\t0.0\n"
    );
    // The refusals before the report: a first row with `sacc`, `saccver` (the default fields)
    // or `staxid`; a later one has a result.
    let spec = |fields: &str| {
        blastp_spec(
            ReportFormat::Tabular,
            &FormatChoice {
                number: 6,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap()
    };
    let thrown = LookupFailure::thrown();
    assert!(failed_lookup_stop(None, thrown, true).is_some());
    assert!(failed_lookup_stop(None, thrown, false).is_none());
    assert!(failed_lookup_stop(Some(&spec("qseqid saccver")), thrown, true).is_some());
    assert!(failed_lookup_stop(Some(&spec("qseqid saccver")), thrown, false).is_none());
    assert!(failed_lookup_stop(Some(&spec("staxid")), thrown, true).is_some());
    assert!(failed_lookup_stop(Some(&spec("qseqid sseqid")), thrown, true).is_none());
    assert!(may_stop_at_failed_lookup(Some(&spec("staxid ssciname"))));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:828-875 (the catch of the subject block) for the default fields (probe of
// line DBID pfoff-2, `R/pfoff/p2`, `*.U2.dflt.6`; review 2, M3): the one writer of the
// default tables' rows of a subject whose load throws (BLASTP's, BLASTN's, TBLASTN's and
// TBLASTX's, `TabularRows::write_thrown_default`; BLASTX's are the `std` fields of
// `TabularWriter::row`) prints `saccver` of the query's last loaded row
// (`TabularRows::found`), the counts and coordinates 0 and the scores; a new
// query's formatter keeps none. The end of the run (`run_end`) writes the rows before it,
// ends the reads and gives NCBI's exception (exit status 2).
#[test]
fn the_default_row_of_a_throw_keeps_the_last_loaded_subject() {
    let config = crate::report::outfmt6::OutputConfig::ncbi_compat();
    let saccver = |index: usize| ["AAA00001.1", "BBB00002.1"][index];
    let mut rows = TabularRows::default();
    assert!(rows.kept_id().is_none());
    rows.found(1);
    assert_eq!(rows.kept_id(), Some(1));
    let mut out = Vec::new();
    rows.write_thrown_default(&mut out, b"q1", saccver, 1e-30, 55.4, &config)
        .unwrap();
    rows.found(0);
    rows.write_thrown_default(&mut out, b"q1", saccver, 0.5, 120.7, &config)
        .unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "q1\tBBB00002.1\t0.000\t0\t0\t0\t0\t0\t0\t0\t1.00e-30\t55.4\n\
             q1\tAAA00001.1\t0.000\t0\t0\t0\t0\t0\t0\t0\t0.50\t120\n"
    );
    let exception =
        crate::report::db_subject::LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO)
            .exception()
            .unwrap();
    let (mut written, mut ended) = (Vec::new(), false);
    let error = run_end(&mut written, b"rows before\n", &exception, || ended = true);
    assert_eq!(written, b"rows before\n");
    assert!(ended);
    let native = error.downcast::<crate::cli::NativeError>().unwrap();
    assert_eq!(native.exit, 2);
    assert!(
        String::from_utf8_lossy(&native.message).contains("DB contains no accession info."),
        "{}",
        String::from_utf8_lossy(&native.message)
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-898: which rows end the run where the load throws (`slen`, the
// letters after the subject block), which are written (a field of that block), and which
// read nothing; `nident` alone reads only on rows with identities.
#[test]
fn a_rows_fields_say_what_its_throw_does() {
    use crate::blastinput::app::Program;
    use crate::report::db_subject::{AfterThrow, FormatThrow, RunThrow};
    let spec = |fields: &str| {
        blastp_spec(
            ReportFormat::CommaSeparatedValues,
            &FormatChoice {
                number: 10,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap()
    };
    let row = |fields: &str, program, ident| spec(fields).row_throw(program, ident);
    assert_eq!(
        row("qseqid stitle slen", Program::Blastp, 1),
        Some(AfterThrow::WritesRow)
    );
    assert_eq!(
        row("qseqid slen", Program::Blastn, 1),
        Some(AfterThrow::EndsRun)
    );
    assert_eq!(
        row("qseqid sseq", Program::Tblastx, 1),
        Some(AfterThrow::EndsRun)
    );
    assert_eq!(
        row("qseqid nident", Program::Blastx, 1),
        Some(AfterThrow::EndsRun)
    );
    assert_eq!(row("qseqid nident", Program::Blastx, 0), None);
    assert_eq!(row("qseqid nident", Program::Blastn, 1), None);
    assert_eq!(row("qseqid evalue qcovs", Program::Tblastn, 1), None);
    assert_eq!(
        spec("qseqid nident").format_throw(Program::Blastx),
        FormatThrow::Writes(RunThrow {
            rows: SubjectFetch::WithIdentities,
            ..RunThrow::every_row(AfterThrow::EndsRun)
        })
    );
    // `slen` and a field of the alignment block read the subject at every row.
    for every in ["qseqid nident slen", "qseqid nident qstart"] {
        assert_eq!(spec(every).throw_rows(Program::Blastx), SubjectFetch::Every);
        assert_eq!(
            spec(every).format_throw(Program::Blastx),
            FormatThrow::Writes(RunThrow::every_row(AfterThrow::EndsRun)),
            "{every}"
        );
        assert_eq!(
            spec(every).row_throw(Program::Blastx, 0),
            Some(AfterThrow::EndsRun)
        );
    }
    // The organism report (outfmt 18) leaves out the HSPs of a throwing subject and goes on.
    let report = blastp_spec(
        ReportFormat::TaxFormat,
        &FormatChoice {
            number: 18,
            spec: String::new(),
            delimiter: String::new(),
        },
    )
    .unwrap();
    assert_eq!(
        report.format_throw(Program::Blastp),
        FormatThrow::Writes(RunThrow::every_row(AfterThrow::WritesRow))
    );
    assert_eq!(
        spec("qseqid evalue").format_throw(Program::Blastx),
        FormatThrow::Unread
    );
}

// NCBI BLAST+ 2.17.0 `blastn -db` of a sequence with an empty title and a taxid whose
// name has 0x81 after one of Homo sapiens (probe `t.tx_bh.blastn.f6all|f6wide`): without
// `sseqid` the lookup fails in `x_SetSubjectIds`, which clears `m_SubjectIds` first, so the
// Seq-align id is the list alone and `slen` is set; the taxonomy is the last found row's.
#[test]
fn a_database_lookup_failing_in_x_set_subject_ids_keeps_the_length() {
    let homo = TabularTaxonomy {
        taxid: 9606,
        taxids: vec![9606],
    };
    let none = TabularTaxonomy::default();
    let titles: [&[u8]; 1] = [b""];
    let subject = |id: &'static [u8], taxonomy, failed: bool, index| TabularSubject {
        id,
        accession: id,
        accession_version: id,
        gi: 0,
        all_ids: &[],
        length: 660,
        titles: Some(&titles),
        taxonomy,
        names: None,
        nucleotide: true,
        qcovs: 100,
        qcovus: -1,
        lookup: if failed {
            SubjectLookup::Failed(FailedLookup {
                full: b"gnl|BL_ORD_ID|0",
                accession: b"BL_ORD_ID:0",
                subject_id: true,
                all_ids: Some(&[]),
                empty_handle: false,
            })
        } else {
            SubjectLookup::Found
        },
        index,
    };
    let found = subject(b"Homo", &homo, false, 1);
    let failed = subject(b"", &none, true, 0);
    assert_eq!(
        failed_rows("qseqid sallseqid sallacc slen staxids", &[found, failed]),
        "nq1\tHomo\tHomo\t660\t9606\nnq1\tgnl|BL_ORD_ID|0\tBL_ORD_ID:0\t660\t9606\n"
    );
    assert_eq!(
        failed_rows("qseqid sseqid slen staxids", &[found, failed]),
        "nq1\tHomo\t660\t9606\nnq1\tUnknown\t0\t9606\n"
    );
}

// NCBI BLAST+ 2.17.0 `blastn -db` of a `-parse_seqids` sequence with the definition lines
// `gb|AB000001.1|` and `lcl|Subject_9`, an empty title and a taxid whose name has 0x81 (probe
// `T/t7/review` `B.gbs9.blastn.f6wide`): `SetSubjectId` finds the Bioseq's `gb` id, and
// `x_SetSubjectIds` throws at the second line after it pushed the first one's ids, so `sall*`
// are the first line's and the Seq-align id; `slen` is set, the counts are `x_ResetFields`'.
#[test]
fn a_lookup_failing_after_the_first_definition_line_keeps_its_ids() {
    let none = TabularTaxonomy::default();
    let titles: [&[u8]; 2] = [b"", b""];
    let before = [DeflineIds {
        id: b"gb|AB000001.1|".to_vec(),
        gi: 0,
        accession: b"AB000001".to_vec(),
    }];
    let subject = TabularSubject {
        id: b"gb|AB000001.1|",
        accession: b"AB000001",
        accession_version: b"AB000001.1",
        gi: 0,
        all_ids: &before,
        length: 1350,
        titles: Some(&titles),
        taxonomy: &none,
        names: None,
        nucleotide: true,
        qcovs: 99,
        qcovus: -1,
        lookup: SubjectLookup::Failed(FailedLookup {
            full: b"gb|AB000001.1|",
            accession: b"AB000001",
            subject_id: false,
            all_ids: Some(&before),
            empty_handle: false,
        }),
        index: 0,
    };
    assert_eq!(
        failed_rows(
            "qseqid sseqid slen sallseqid sallgi sallacc evalue",
            &[subject; 2]
        ),
        "nq1\tgb|AB000001.1|\t1350\tgb|AB000001.1|;gb|AB000001.1|\t0;0\tAB000001;AB000001\t0.0\n"
            .repeat(2)
    );
    // `saccver` is the Bioseq's: no stop where only `x_SetSubjectIds` throws.
    assert_eq!(
        failed_rows("qseqid saccver sallseqid", &[subject]),
        "nq1\tAB000001.1\tgb|AB000001.1|;gb|AB000001.1|\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:174-185 (`FindBestChoice` of the empty `m_SubjectId`, a null pointer:
// NCBI BLAST+ 2.17.0 exits 255) and :262-268 (`staxid` of an `m_SubjectTaxId` that no row set:
// a different number in each run, probe of line RG T7 `6 qseqid sseqid staxid`).
#[test]
fn failed_lookups_stop_where_ncbi_has_no_result() {
    let spec = |fields: &str| {
        blastp_spec(
            ReportFormat::CommaSeparatedValues,
            &FormatChoice {
                number: 10,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap()
    };
    let stop = |fields: Option<&str>, all_ids: bool, first: bool| {
        let spec = fields.map(spec);
        let failure = LookupFailure {
            subject_id: true,
            all_ids,
            bioseq: false,
            length: false,
        };
        failed_lookup_stop(spec.as_ref(), failure, first).is_some()
    };
    assert!(stop(None, false, false));
    assert!(stop(Some(""), false, false));
    assert!(stop(Some("qseqid sacc"), false, false));
    assert!(stop(Some("qseqid saccver"), true, false));
    assert!(!stop(Some("qseqid sseqid evalue"), false, true));
    assert!(stop(Some("qseqid sseqid staxid"), false, true));
    assert!(!stop(Some("qseqid sseqid staxid"), false, false));
    assert!(!stop(Some("qseqid sallseqid staxid"), false, true));
    assert!(stop(Some("qseqid sallseqid staxid"), true, true));
    // `SetSubjectId` finds the ids where only `x_SetSubjectIds` throws: no null pointer.
    let only_all = LookupFailure {
        subject_id: false,
        all_ids: true,
        bioseq: false,
        length: false,
    };
    assert!(failed_lookup_stop(Some(&spec("qseqid saccver sallseqid")), only_all, false).is_none());
    assert!(failed_lookup_stop(None, only_all, true).is_none());
    assert!(!may_stop_at_failed_lookup(Some(&spec(
        "qseqid sseqid sallseqid evalue"
    ))));
    assert!(may_stop_at_failed_lookup(None));
    let record = |s_idx| (s_idx == 1).then(|| LookupFailure::of_plain(false));
    let result = check_failed_lookups(
        Some(&spec("qseqid sseqid staxid")),
        [vec![0, 1], vec![1, 0]],
        record,
        |s_idx| failed_record_title(s_idx, b"Subject_9"),
        "BLASTN",
    );
    assert!(result
        .unwrap_err()
        .to_string()
        .starts_with("subject record 2 (Subject_9) has a title"));
    assert!(check_failed_lookups(
        Some(&spec("qseqid sseqid staxid")),
        [vec![0, 1], vec![0, 1]],
        record,
        |s_idx| failed_record_title(s_idx, b"x"),
        "BLASTN",
    )
    .is_ok());
}

#[test]
fn btop_counts_runs() {
    use crate::common::GapEditOp::{Del, Sub};
    let mut out = Vec::new();
    write_btop(&mut out, b"ACGT-A", b"ACCTTA", &[Sub(4), Del(1), Sub(1)]);
    assert_eq!(out, b"2GC1-T1");
}

// tabular.cpp:982-1027 (cited at `push_alignment_row` and `write_btop`): a row's end gaps (its
// gap operations of the script before its first letter and after its last, `row_end_gaps`) are
// spaces, its inner gaps `-`; an end-gap column of `btop` is a letter and a space.
#[test]
fn end_gaps_are_spaces_in_the_rows_and_btop() {
    use crate::common::GapEditOp::{Del, Ins, Sub};
    use crate::report::serial::seqalign::AlignRow::{Query, Subject};
    let mut out = Vec::new();
    let script = [Del(3), Sub(2), Del(1), Sub(2), Del(2)];
    push_alignment_row(&mut out, b"---AC-DE--", &script, Query);
    assert_eq!(out, b"   AC-DE  ");
    out.clear();
    push_alignment_row(&mut out, b"MKV-L", &[Sub(3), Del(1), Sub(1)], Query);
    assert_eq!(out, b"MKV-L");
    out.clear();
    let script = [Ins(3), Sub(2), Del(1), Sub(1)];
    push_alignment_row(&mut out, b"---ACKD", &script, Subject);
    assert_eq!(out, b"   ACKD");
    out.clear();
    write_btop(&mut out, b"MQLAC-D", b"---ACKD", &script);
    assert_eq!(out, b"M Q L 2-K1");
    out.clear();
    write_btop(
        &mut out,
        b"AC-DKK",
        b"ACGD--",
        &[Sub(2), Del(1), Sub(1), Ins(2)],
    );
    assert_eq!(out, b"2-G1K K ");
    // Rows whose `-` fill a segment, without a script (BLASTX's mapped row past a query's end,
    // db.x2g.bx.s1o.subj.6w and .10: NCBI's `TCGGCATGGCTG-` and BTOP `...-X`).
    out.clear();
    push_alignment_row(&mut out, b"TCG-", &[], Query);
    assert_eq!(out, b"TCG-");
    out.clear();
    write_btop(&mut out, b"TCG-", b"TANX", &[]);
    assert_eq!(out, b"1CAGN-X");
    // A `-` that is a residue (a PSSM query's gap residue, a database's stored gap character) is
    // a letter of the script, no end gap (M4a-1c LG probe `ts_std0first.f2`: NCBI's qseq
    // `-ASLG...` and BTOP `-M23...` for a query whose first residue is the gap residue).
    out.clear();
    push_alignment_row(&mut out, b"-ASL-", &[Sub(4), Del(1)], Query);
    assert_eq!(out, b"-ASL ");
    out.clear();
    write_btop(&mut out, b"-ASL", b"MASL", &[Sub(4)]);
    assert_eq!(out, b"-M3");
    out.clear();
    write_btop(&mut out, b"MASLK", b"MASL-", &[Sub(5)]);
    assert_eq!(out, b"4K-");
}

// Review S1 of PS15a-2 (`tmp/rev15a2/col`, every field alone against NCBI with a profile whose id
// is `lcl|Query_1`): the RPS programs refuse outfmt 18 and 20 only where NCBI's bytes differ.
#[test]
fn an_rps_report_refuses_a_generated_id_where_ncbi_shows_the_query() {
    use crate::blastinput::app::Program;
    let spec = |number: i32, format: ReportFormat, fields: &str| {
        TabularSpec::of_app(
            format,
            &FormatChoice {
                number,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
            crate::report::archive::options::RequestKind::Rpsblast,
        )
        .unwrap()
    };
    let csv = |fields: &str| spec(20, ReportFormat::CommaSeparatedValuesWithHeader, fields);
    let refuses = |fields: &str, program, ident, taxid| {
        csv(fields).refuses_generated_id(program, ident, taxid)
    };
    for program in [Program::Blastp, Program::Blastx] {
        // The profile's own ids and the scores: NCBI's bytes.
        for fields in [
            "qseqid sseqid sacc saccver sgi sallgi evalue",
            "qseqid bitscore score qcovs qcovus",
        ] {
            assert!(!refuses(fields, program, 5, true), "{fields} {program:?}");
        }
        // What NCBI finds by the id beyond it: refused.
        for fields in [
            "qseqid sallseqid",
            "qseqid sallacc",
            "qseqid stitle",
            "qseqid salltitles",
            "qseqid staxids",
            "qseqid slen",
            "qseqid sseq",
            "qseqid pident",
            "",
        ] {
            assert!(refuses(fields, program, 5, false), "{fields} {program:?}");
        }
        // Taxonomy: only where the profile has a taxid.
        for fields in ["qseqid staxid", "qseqid ssciname", "qseqid sskingdoms"] {
            assert!(!refuses(fields, program, 5, false), "{fields} {program:?}");
            assert!(refuses(fields, program, 5, true), "{fields} {program:?}");
        }
        // `nident` where the block opens for it (identities).
        assert!(!refuses("qseqid nident", program, 0, false));
        assert!(refuses("qseqid nident", program, 1, false));
        // The organism report: only a profile with a taxid; outfmt 6, 7 and 10 are line RG X7's.
        let report = spec(18, ReportFormat::TaxFormat, "");
        assert!(!report.refuses_generated_id(program, 5, false));
        assert!(report.refuses_generated_id(program, 5, true));
        for (number, format) in [
            (6, ReportFormat::Tabular),
            (7, ReportFormat::TabularWithComments),
            (10, ReportFormat::CommaSeparatedValues),
        ] {
            assert!(!spec(number, format, "qseqid stitle").refuses_generated_id(program, 5, true));
        }
    }
    // RPSBLAST: the query's letters and the strand; RPSTBLASTN: the alignment block throws.
    for fields in ["qseqid qseq", "qseqid sstrand", "qseqid qstart length"] {
        assert!(!refuses(fields, Program::Blastp, 5, false), "{fields}");
        assert!(refuses(fields, Program::Blastx, 5, false), "{fields}");
    }
}

/// The exception of a query of the molecule `other` names the coding of its Seq-data (NCBI's
/// object manager converts it to the coding 1).
#[test]
fn the_sequence_exception_names_the_codings_of_the_query() {
    for coding in [2, 8, 10] {
        let message = query_sequence_exception(coding).get_msg().to_vec();
        assert_eq!(
            message,
            format!("Incompatible sequence codings: {coding} -> 1").into_bytes()
        );
    }
}

/// A table row counts only the columns inside the query row's aligned range
/// (`TabularHsp::within_query_range` of its script, NCBI's `GetAlignLengths`): the query row's end
/// gaps of a Smith-Waterman alignment come off `length`, `gaps` and `gapopen`; a subject row's end
/// gaps and the inner gaps of an X-drop script stay.
#[test]
fn table_counts_leave_out_the_query_rows_end_gaps() {
    use crate::common::GapEditOp::{self, Del, Ins, Sub};
    let counts = |length, gap_opens, gaps, ops: &[GapEditOp]| {
        let hsp = TabularHsp {
            q_start: 1,
            q_end: 190,
            s_start: 1,
            s_end: 285,
            evalue: 1e-100,
            bit_score: 380.0,
            score: 990,
            length,
            gap_opens,
            gaps,
            num_ident: 190,
            q_frame: 0,
            s_frame: 0,
            rows: None,
            ops,
            qcovhsp: 100,
        }
        .within_query_range();
        (hsp.length, hsp.gaps, hsp.gap_opens)
    };
    // The probe's alignment: 95 gaps before the query's first letter.
    assert_eq!(counts(285, 1, 95, &[Del(95), Sub(190)]), (190, 0, 0));
    // A gap end in the query's row, with an inner gap that stays.
    assert_eq!(
        counts(20, 2, 6, &[Sub(10), Ins(2), Sub(4), Del(4)]),
        (16, 2, 1)
    );
    // The subject row's end gaps are inside the query row's range.
    assert_eq!(counts(14, 2, 4, &[Ins(3), Sub(10), Ins(1)]), (14, 4, 2));
    // An X-drop script (a substitution at each end) and an ungapped HSP keep their counts.
    let xdrop = [Sub(10), Ins(1), Sub(5), Del(2), Sub(3)];
    assert_eq!(counts(21, 2, 3, &xdrop), (21, 3, 2));
    assert_eq!(counts(190, 0, 0, &[]), (190, 0, 0));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-1028
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1233-1245
// Without a subject-block field, an empty subject handle is not caught in `SetFields`' subject
// block; NCBI throws in `CAlnVec::GetWholeAlnSeqString` on the subject row (or `GetBioseqHandle`
// in `GetFrame` for translated-subject minus-strand HSPs). With a subject-block field (such as
// `sseqid`), NCBI catches the empty handle and prints its fallback row.
#[test]
fn empty_handle_tabular_fields_stop_only_without_subject_block_catch() {
    use crate::blastinput::app::Program;
    let spec = |fields: &str| {
        blastp_spec(
            ReportFormat::Tabular,
            &FormatChoice {
                number: 6,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap()
    };
    let empty = LookupFailure::empty_handle();
    let plain = LookupFailure::of_plain(true);
    // Subject-block field catches the empty handle before the alignment block runs.
    assert!(!spec("sseqid sseq").empty_handle_throws(Program::Blastp));
    assert!(!spec("sseqid qseq").empty_handle_throws(Program::Blastn));
    assert!(!spec("stitle pident").empty_handle_throws(Program::Blastx));
    assert!(!spec("sseqid sframe qstart").empty_handle_throws(Program::Tblastx));
    assert!(
        failed_lookup_stop_for(Some(&spec("sseqid sseq")), empty, true, Program::Blastp).is_none()
    );
    // Without a subject-block field, sequence/scoring fields throw on empty handle.
    for fields in [
        "qseqid sseq",
        "qseqid qseq",
        "qseqid positive",
        "qseqid ppos",
        "qseqid btop",
    ] {
        for program in [
            Program::Blastp,
            Program::Blastn,
            Program::Blastx,
            Program::Tblastn,
            Program::Tblastx,
        ] {
            assert!(
                spec(fields).empty_handle_throws(program),
                "{fields} {program:?}"
            );
            assert!(
                failed_lookup_stop_for(Some(&spec(fields)), empty, false, program).is_some(),
                "{fields} {program:?}"
            );
            assert!(
                failed_lookup_stop_for(Some(&spec(fields)), plain, false, program).is_none(),
                "{fields} {program:?}"
            );
        }
    }
    // `pident`/`nident`/`mismatch` throw outside BLASTN, and do not throw in BLASTN alone.
    for fields in ["qseqid pident", "qseqid nident", "qseqid mismatch"] {
        assert!(!spec(fields).empty_handle_throws(Program::Blastn));
        for program in [
            Program::Blastp,
            Program::Blastx,
            Program::Tblastn,
            Program::Tblastx,
        ] {
            assert!(
                spec(fields).empty_handle_throws(program),
                "{fields} {program:?}"
            );
        }
    }
    // `sframe`/`frames` with an `opens_counts` field throw for translated subjects only.
    for fields in ["qseqid qstart sframe", "qseqid length frames"] {
        assert!(spec(fields).empty_handle_throws(Program::Tblastn));
        assert!(spec(fields).empty_handle_throws(Program::Tblastx));
        assert!(!spec(fields).empty_handle_throws(Program::Blastp));
        assert!(!spec(fields).empty_handle_throws(Program::Blastn));
        assert!(!spec(fields).empty_handle_throws(Program::Blastx));
    }
    // Standalone `sframe`/`frames` never open the alignment block.
    assert!(!spec("qseqid sframe frames").empty_handle_throws(Program::Tblastx));
    assert!(!spec("qseqid sframe frames").empty_handle_throws(Program::Tblastn));
    assert!(may_stop_at_failed_lookup_for(
        Some(&spec("qseqid sseq")),
        Program::Blastp
    ));
    assert!(!may_stop_at_failed_lookup_for(
        Some(&spec("qseqid evalue")),
        Program::Blastp
    ));
}
