use super::*;

#[test]
fn v1_output_format_reads_the_frozen_numbers() {
    use crate::report::OutputFormat;
    let parse = super::v1_blastp::v1_output_format;
    assert_eq!(parse("0").unwrap(), (OutputFormat::Pairwise, None));
    assert_eq!(parse("6").unwrap(), (OutputFormat::Tabular, None));
    assert_eq!(
        parse("7").unwrap(),
        (OutputFormat::TabularWithComments, None)
    );
    // With custom fields
    let (fmt, fields) = parse("6 qaccver saccver").unwrap();
    assert_eq!(fmt, OutputFormat::Tabular);
    assert_eq!(fields, Some("qaccver saccver".to_string()));
    // Invalid format; ABI v1 reads `06` as no format.
    assert!(parse("99").is_err());
    assert!(parse("06").is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:485-651
#[test]
fn fasta_store_allocates_and_releases_handles() {
    let mut store = FastaStore::default();

    let first = store
        .insert(1, "query".to_string(), b">q1\nACGT\n".to_vec())
        .expect("first handle");
    let second = store
        .insert(1, "subject".to_string(), b">s1\nACGT\n".to_vec())
        .expect("second handle");

    assert_eq!(first, 1);
    assert_eq!(second, 2);
    assert!(store.entries.contains_key(&first));
    store.release(first).expect("release existing handle");
    assert!(!store.entries.contains_key(&first));
    assert!(store.release(first).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
#[test]
fn tblastx_web_args_reject_unimplemented_outfmt() {
    let parse = |values: &[&str]| {
        parse_tblastx_args(values, PathBuf::new(), PathBuf::new(), PathBuf::new())
    };

    assert_eq!(parse(&["-outfmt", "6"]).expect("outfmt 6").outfmt, "6");
    assert_eq!(parse(&["-outfmt=6"]).expect("outfmt=6").outfmt, "6");
    assert!(parse(&["-outfmt", "0"]).is_err());
    assert!(parse(&["-outfmt=7"]).is_err());
    assert!(parse(&["-outfmt", "6 qseqid sseqid"]).is_err());
}

// Plan TD-1: ABI v1 is frozen, so its BLASTP checks the options before the output
// format and refuses a tabular field that it does not write, as before S08+ (the CLI
// and ABI v2 follow NCBI: -outfmt first, an unknown token ignored).
#[test]
fn blastp_web_pair_keeps_v1_error_order_and_fields() {
    let args = |outfmt: &str, extra: &[&str]| {
        let mut words = vec!["-outfmt", outfmt];
        words.extend_from_slice(extra);
        parse_blastp_args(&words, PathBuf::new(), PathBuf::new(), PathBuf::new())
            .expect("blastp web args")
    };
    let fasta = ">q\nMKVLAAGIVGLLLAQPAMAAEIPVDPALAV\n";
    let error = |outfmt: &str, extra: &[&str]| {
        engine_error(v1_blastp::run_web_pair(args(outfmt, extra), fasta, fasta).unwrap_err())
    };
    assert_eq!(
        error("6 nosuch", &[]),
        "unsupported blastp tabular field 'nosuch'"
    );
    assert_eq!(
        error("9", &[]),
        "Unsupported output format: 9. Supported: 0, 6, 7"
    );
    assert_eq!(
        error("9", &["-num_threads", "0"]),
        "num_threads must be greater than zero"
    );
    assert!(!error("9", &["-ungapped"]).contains("output format"));
    assert!(error("6 nosuch", &["-matrix", "FOO"]).contains("FOO"));
    assert!(v1_blastp::run_web_pair(args("6 qseqid sseqid std", &[]), fasta, fasta).is_ok());
}

// Plan TD-1: ABI v1 is frozen, so its BLASTN keeps outfmt 6 and 7 and rejects
// outfmt 0, which the engine implements since LOSAT Web session S07, with the error
// that it gave before.
#[test]
fn blastn_web_pair_keeps_rejecting_outfmt_0() {
    let args = |outfmt: &str| {
        parse_blastn_args(
            &["-outfmt", outfmt],
            PathBuf::new(),
            PathBuf::new(),
            PathBuf::new(),
        )
        .expect("blastn web args")
    };
    let fasta = ">q\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let error = v1_blastn::run_web_pair(args("0"), fasta, fasta).unwrap_err();
    assert_eq!(engine_error(error), "unsupported BLASTN output format: 0");
    assert!(v1_blastn::run_web_pair(args("6"), fasta, fasta).is_ok());
    // Plan TD-1: LOSAT's limits come after NCBI's checks and before the records; an
    // empty query gives the empty report of NCBI and of v1 before S07+.
    let with = |extra: &[&str]| {
        let mut words = vec!["-outfmt", "6"];
        words.extend_from_slice(extra);
        parse_blastn_args(&words, PathBuf::new(), PathBuf::new(), PathBuf::new())
            .expect("blastn web args")
    };
    let empty_record = ">s0\n>s1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    assert!(v1_blastn::run_web_pair(with(&[]), "", empty_record).is_ok());
    // An empty query gives the empty report whatever the subject's deflines or LOSAT's
    // limits (as before S07+ and in NCBI); the serial build still checks threads first.
    let tab_defline = ">s0\tx\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    assert!(v1_blastn::run_web_pair(with(&[]), "", tab_defline).is_ok());
    assert!(v1_blastn::run_web_pair(with(&["-evalue", "1e400"]), "", fasta).is_ok());
    #[cfg(not(feature = "parallel"))]
    assert!(v1_blastn::run_web_pair(with(&["-num_threads", "2"]), "", fasta).is_err());
    // A subject without residues gets its warning and is not searched, as in the CLI.
    assert!(v1_blastn::run_web_pair(with(&[]), fasta, empty_record).is_ok());
    // NCBI's CArg_Double reads a signed infinity and NaN, not a bare `inf`, and NCBI
    // searches with them (E2g).
    let error = parse_blastn_args(
        &["-outfmt", "6", "-evalue", "inf"],
        PathBuf::new(),
        PathBuf::new(),
        PathBuf::new(),
    )
    .unwrap_err();
    assert!(error.contains("expected a number"), "{error}");
    assert!(v1_blastn::run_web_pair(with(&["-evalue", "+inf"]), fasta, fasta).is_ok());
    assert!(v1_blastn::run_web_pair(with(&["-evalue=+nan"]), fasta, fasta).is_ok());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:338-340
// Plan TD-1, revised (DW-48): ABI v1 reads its inputs with NCBI's reader, as the CLI: the
// same report for what the reader removes or changes, and every defline that it reads.
#[test]
fn inputs_are_read_by_ncbis_reader() {
    let dna = "GGATCCTTAGCAATCGGCTAAGCTTCCGATACGTTAGGCATCCAGTTCAGGATTACAGGCTTAACCGGT";
    let protein = "MKVLAAGIVGLLLAQPAMAAEIPVDPALAVKLTSGQTVEFRWDNGHYCPIM";
    for (program, letters) in [("blastn", dna), ("tblastx", dna), ("blastp", protein)] {
        let (head, tail) = letters.split_at(30);
        let subject = format!(">s1\n{letters}\n");
        let report = |query: &str| run_pair(program, query, &subject, "6", "");
        let expected = report(&format!(">q1 first\n{letters}\n")).expect("report");
        assert!(!expected.is_empty(), "{program}");
        for query in [
            format!(">q1 first\n{head}-{tail}\n"),
            format!(">q1 first\r\n{head}\r\n\r\n{tail}\r\n"),
            format!(";comment\n>q1 first\n{head} {tail}\n"),
            format!("  \n>q1 first\n{letters}\n>q2\n"),
        ] {
            assert_eq!(
                report(&query).as_ref(),
                Ok(&expected),
                "{program}: {query:?}"
            );
        }
        // A subject without residues among others gets its warning and counts in the
        // statistics of the search (the E-value differs), as in the CLI.
        let subjects = format!(">s0\n{subject}");
        let query = format!(">q1 first\n{letters}\n");
        let with_empty = run_pair(program, &query, &subjects, "6", "").expect("report");
        assert!(with_empty.starts_with(b"q1\ts1\t"), "{program}");
        // White space only is NCBI's `Query is Empty!`; comment lines only are not empty.
        assert_eq!(report(" \n"), Ok(Vec::new()), "{program}");
        let error = report(";comment\n").unwrap_err();
        assert!(error.contains("Empty CBlastQueryVector"), "{error}");
        // A defline that is not ASCII is read as bytes.
        let named = report(&format!(">\u{3042}1\n{letters}\n")).expect("report");
        assert!(named.starts_with("\u{3042}1\t".as_bytes()), "{program}");
    }
}

// BLASTP's handles read the stored bytes as `run_pair` reads the call's, with the labels.
#[test]
fn blastp_handles_read_the_stored_bytes_as_the_call_does() {
    let fasta = b">q\nMKVLAAGIVGLLLAQPAMAAEIPVDPALAVKLTSGQTVEFRWDNGHYCPIM\n".to_vec();
    let text = std::str::from_utf8(&fasta).expect("UTF-8");
    let expected = run_pair("blastp", text, text, "6", "");
    assert!(expected.as_ref().is_ok_and(|report| !report.is_empty()));
    let stored = |label: &str| {
        fasta_store()
            .lock()
            .expect("store")
            .insert(
                count_stored_records(&fasta).expect("records"),
                label.to_string(),
                fasta.clone(),
            )
            .expect("handle")
    };
    let (query, subject) = (stored("query"), stored("subject"));
    assert_eq!(
        run_pair_handles("blastp", query, subject, "6", ""),
        expected
    );
    let mut store = fasta_store().lock().expect("store");
    store.release(query).expect("release");
    store.release(subject).expect("release");
}

// The store counts the records of NCBI's reader: white space and comments have none, and
// the letters before the first defline are a record (`fDLOptional`).
#[test]
fn the_store_counts_the_records_of_ncbis_reader() {
    assert_eq!(count_stored_records(b""), Ok(0));
    assert_eq!(count_stored_records(b" \n\n"), Ok(0));
    assert_eq!(count_stored_records(b";comment\n"), Ok(0));
    assert_eq!(count_stored_records(b">a\nACGT\n>b\nMKV\n"), Ok(2));
    assert_eq!(count_stored_records(b"ACGT\n"), Ok(1));
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:213-219
// TBLASTN (added on the Owner's request, 2026-11-04): both exports return the bytes of the
// shared vector (`tblastn::web::vector`), which `tests/web_v1_tblastn.rs` holds to the CLI.
#[test]
fn tblastn_web_pair_and_handles_return_the_cli_bytes_of_the_shared_vector() {
    use crate::algorithm::tblastn::web::vector;
    let stored = |file: &str, label: &str| {
        let bytes = vector::read(file);
        fasta_store()
            .lock()
            .expect("store")
            .insert(
                // The store takes only inputs with one record or more.
                count_stored_records(&bytes).unwrap_or_default(),
                label.to_string(),
                bytes,
            )
            .ok()
    };
    for case in vector::cases() {
        let extra = case.extra.join("\0");
        let expected = case
            .expected
            .map_err(|message| String::from_utf8(message).expect("UTF-8 message"));
        // `run_pair` names the subject `subject`.
        if case.label == "subject" {
            let (query, subject) = (vector::read(&case.query), vector::read(&case.subject));
            assert_eq!(
                run_pair(
                    "tblastn",
                    std::str::from_utf8(&query).expect("UTF-8 query"),
                    std::str::from_utf8(&subject).expect("UTF-8 subject"),
                    &case.outfmt,
                    &extra
                ),
                expected,
                "{} (run_pair)",
                case.name
            );
        }
        // The handles of the inputs that the store takes (not the blank ones).
        let handles = (
            stored(&case.query, "query"),
            stored(&case.subject, &case.label),
        );
        if let (Some(query_handle), Some(subject_handle)) = handles {
            assert_eq!(
                run_pair_handles(
                    "tblastn",
                    query_handle,
                    subject_handle,
                    &case.outfmt,
                    &extra
                ),
                expected,
                "{} (run_pair_handles)",
                case.name
            );
        }
        let mut store = fasta_store().lock().expect("store");
        for handle in [handles.0, handles.1].into_iter().flatten() {
            store.release(handle).expect("release");
        }
    }
    let query = vector::read("query.faa");
    let subject = vector::read("subject.fna");
    let (query_text, subject_text) = (
        std::str::from_utf8(&query).expect("UTF-8 query"),
        std::str::from_utf8(&subject).expect("UTF-8 subject"),
    );
    // An empty `outfmt` is 6 (`web_outfmt_or_default`), as for the other programs.
    let fmt6 = vector::read("fmt6.out");
    assert_eq!(
        run_pair("tblastn", query_text, subject_text, " ", ""),
        Ok(fmt6)
    );
    // The words that name files or a database are refused before the search.
    assert_eq!(
        run_pair("tblastn", query_text, subject_text, "6", "-db\0nt"),
        Err("-db (a database search) is not supported in the browser yet".to_string())
    );
    assert_eq!(
        run_pair("tblastn", query_text, subject_text, "6", "-out\0x"),
        Err("file/help arguments are unsupported for the TBLASTN web API".to_string())
    );
    // `-help-full` asks for the CLI's help, as `-help` (line RG item E2).
    assert_eq!(
        run_pair("tblastn", query_text, subject_text, "6", "-help-full"),
        Err("file/help arguments are unsupported for the TBLASTN web API".to_string())
    );
    // The serial module has no threads to give (as the CLI built without them).
    #[cfg(not(feature = "parallel"))]
    assert!(run_pair(
        "tblastn",
        query_text,
        subject_text,
        "6",
        "-num_threads\x002"
    )
    .is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
#[test]
fn fasta_store_rejects_empty_and_unknown_handles() {
    let mut store = FastaStore::default();

    assert!(store.insert(0, "empty".to_string(), Vec::new()).is_err());
    assert!(store.release(42).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
#[test]
fn fasta_store_clear_removes_entries_and_resets_handles() {
    let mut store = FastaStore::default();
    let handle = store
        .insert(1, "query".to_string(), b">q1\nACGT\n".to_vec())
        .expect("stored handle");

    assert_eq!(handle, 1);
    store.clear();
    assert!(store.entries.is_empty());

    let next = store
        .insert(1, "query2".to_string(), b">q2\nACGT\n".to_vec())
        .expect("stored handle after clear");
    assert_eq!(next, 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238 and 2570-2572
// `-searchsp` and `-dbsize` (added on gbdraw's request, 2026-11-10, so that a comparison's
// E-values have the scale of its genomes): `losat_web_run_pair` returns the bytes of the shared
// vector (`tests/fixtures/web_v1_search_space/`), which `tests/web_v1_search_space.rs` holds to
// the CLI; BLASTP's handles return them too.
#[test]
fn search_space_and_database_size_give_the_cli_bytes_of_the_shared_vector() {
    let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/web_v1_search_space");
    let read = |name: &str| {
        let path = dir.join(name);
        std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
    };
    let cases = String::from_utf8(read("cases.tsv")).expect("UTF-8 cases");
    let mut count = 0;
    for line in cases.lines().filter(|line| !line.starts_with('#')) {
        let [name, program, query, subject, outfmt, extra] =
            line.split('\t').collect::<Vec<_>>()[..]
        else {
            panic!("six columns: {line:?}");
        };
        let (query, subject) = (read(query), read(subject));
        let extra = extra.split_whitespace().collect::<Vec<_>>().join("\0");
        let expected = Ok(read(&format!("{name}.out")));
        assert_eq!(
            run_pair(
                program,
                std::str::from_utf8(&query).expect("UTF-8 query"),
                std::str::from_utf8(&subject).expect("UTF-8 subject"),
                outfmt,
                &extra
            ),
            expected,
            "{name} (run_pair)"
        );
        if program == "blastp" {
            let store = |bytes: &Vec<u8>, label: &str| {
                fasta_store()
                    .lock()
                    .expect("store")
                    .insert(
                        count_stored_records(bytes).expect("records"),
                        label.to_string(),
                        bytes.clone(),
                    )
                    .expect("handle")
            };
            let (query, subject) = (store(&query, "query"), store(&subject, "subject"));
            assert_eq!(
                run_pair_handles("blastp", query, subject, outfmt, &extra),
                expected,
                "{name} (run_pair_handles)"
            );
            fasta_store().lock().expect("store").clear();
        }
        count += 1;
    }
    assert_eq!(count, 12);
}

// The two values are NCBI's `CArg_Int8` (`-dbsize`) and `CArg_Int8` of at least 0
// (`-searchsp`), read by the CLI's value parsers; an invalid one is refused with the parser's
// message, in the words of the other options of ABI v1 (`<flag> parse error: <detail>`).
#[test]
fn invalid_search_space_and_database_size_are_refused_as_the_cli_refuses_them() {
    use crate::blastinput::value_parsers::{ncbi_int8, nonnegative_ncbi_int8};
    let no_path = || (PathBuf::new(), PathBuf::new(), PathBuf::new());
    let parse = |program: &str, words: &[&str]| -> Result<(Option<i64>, Option<i64>), String> {
        let (query, subject, out) = no_path();
        match program {
            "blastn" => parse_blastn_args(words, query, subject, out)
                .map(|args| (args.searchsp, args.dbsize)),
            "tblastx" => parse_tblastx_args(words, query, subject, out)
                .map(|args| (args.searchsp, args.dbsize)),
            _ => parse_blastp_args(words, query, subject, out)
                .map(|args| (args.searchsp, args.dbsize)),
        }
    };
    for program in ["blastn", "tblastx", "blastp"] {
        for (flag, parser) in [
            (
                "-searchsp",
                nonnegative_ncbi_int8 as fn(&str) -> Result<i64, String>,
            ),
            ("-dbsize", ncbi_int8 as fn(&str) -> Result<i64, String>),
        ] {
            for value in ["-1", "abc", "", "1.5", "0x", "9223372036854775808"] {
                let joined = format!("{flag}={value}");
                for words in [vec![flag, value], vec![joined.as_str()]] {
                    let expected = parser(value)
                        .err()
                        .map(|detail| format!("{flag} parse error: {detail}"));
                    let got = parse(program, &words).err();
                    assert_eq!(got, expected, "{program} {words:?}");
                }
            }
            let error = parse(program, &[flag]).unwrap_err();
            assert_eq!(error, format!("{flag} requires a value"), "{program}");
        }
        // The values that NCBI reads: a sign, hexadecimal, 64 bits; 0 is no search space.
        assert_eq!(
            parse(program, &["-searchsp", "0x10", "-dbsize=-3000000000"]),
            Ok((Some(16), Some(-3_000_000_000))),
            "{program}"
        );
        assert_eq!(
            parse(program, &["-searchsp=9223372036854775807"]),
            Ok((Some(i64::MAX), None)),
            "{program}"
        );
        assert_eq!(parse(program, &[]), Ok((None, None)), "{program}");
    }
}
