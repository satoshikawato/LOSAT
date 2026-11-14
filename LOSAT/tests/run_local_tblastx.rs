//! The TBLASTX search of `losat::api` (which calls the shared entry `run_local`) writes several
//! output formats from one search. Every format must be byte-identical to the CLI run that
//! requests only that format, the warnings must equal that run's standard error, and the hit
//! records must be the final HSP list. With `-num_threads 2` the native search reduces its
//! results in a collector thread; the result must not change. (Identity of outfmt 0 and 7 with
//! NCBI is checked by docs/evidence/losat_web_e2a/check_losat.py over the frozen fixtures;
//! identity with the previous release over the regression manifests by
//! docs/evidence/losat_web_e1a/capture_outputs.py.) The entry's own contract (its observer, its
//! refusals, its errors) is checked in the crate, `LOSAT/src/api/local_blast_tests/tblastx.rs`.

mod run_local_support;

use std::process::Command;

use losat::api::{Outcome, Program};
use run_local_support::{cli_out_path, fixture_sequence, Files, TempFasta};

/// The last three are the custom tabular reports of the shared writer (line CLI C2, PR 2).
const FORMATS: [&str; 13] = [
    "0",
    "6",
    "7",
    "8",
    "9",
    "12",
    "5",
    "15",
    "16",
    "11",
    "10 delim=| std qseq sseq sframe",
    "7 std qseq sseq btop positive ppos frames qcovs qcovhsp",
    "20 qseqid sseqid evalue sstrand",
];

/// Two queries against three subjects cut from two genomes; the first query overlaps
/// two subjects.
struct Inputs {
    query: TempFasta,
    subject: TempFasta,
}

impl Inputs {
    fn new() -> Self {
        let a = fixture_sequence("LC738884.fasta");
        let b = fixture_sequence("LC741431.fasta");
        Self {
            query: TempFasta::new(
                "tblastx_query.fna",
                &[("q_a0", &a[0..3_000]), ("q_a100k", &a[100_000..103_000])],
            ),
            subject: TempFasta::new(
                "tblastx_subject.fna",
                &[
                    ("s_a1k", &a[1_000..6_000]),
                    ("s_a99k", &a[99_000..105_000]),
                    ("s_a0", &a[0..2_000]),
                    ("s_b", &b[0..5_000]),
                ],
            ),
        }
    }

    /// The library's search of `formats` with the options `extra`, with the files' paths as the
    /// names that the reports show, as the CLI shows them.
    fn collect(&self, formats: &[&str], extra: &[&str]) -> Outcome {
        Files::new(Program::Tblastx, &self.query.0, &self.subject.0).collect(formats, extra)
    }

    /// The CLI output file and standard error of one single-format run.
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = cli_out_path(&self.query.0, outfmt);
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("tblastx")
            .arg("-query")
            .arg(&self.query.0)
            .arg("-subject")
            .arg(&self.subject.0)
            .args(extra)
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-out")
            .arg(&out)
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI");
        assert!(output.status.success(), "CLI failed: {output:?}");
        let bytes = std::fs::read(&out).expect("CLI output");
        let _ = std::fs::remove_file(&out);
        (bytes, output.stderr)
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new();
    let serial = inputs.collect(&FORMATS, &[]);
    let mut row_counts = Vec::new();
    for extra in [
        &[][..],
        &["-num_threads", "2"][..],
        &["-max_target_seqs", "1"][..],
    ] {
        let together = inputs.collect(&FORMATS, extra);
        let rows = together.reports[1].spans.iter().flatten().count();
        assert_eq!(
            together.hsps.len(),
            rows,
            "one hit record per outfmt 6 row {extra:?}"
        );
        row_counts.push(rows);
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.collect(&[outfmt], extra);
            assert!(
                together.reports[index].bytes == alone.reports[0].bytes,
                "run_collect outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert!(
                together.reports[index].bytes == cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            assert_eq!(
                together.diagnostics,
                cli_stderr,
                "warnings of {} formats equal those of the single-format CLI run {extra:?}",
                FORMATS.len()
            );
        }
        if extra.first() == Some(&"-num_threads") {
            assert!(
                together
                    .reports
                    .iter()
                    .map(|report| &report.bytes)
                    .eq(serial.reports.iter().map(|report| &report.bytes)),
                "the collector thread must not change the result"
            );
        }
    }
    assert!(row_counts[0] > 0, "the fixture must produce hits");
    assert!(
        row_counts[2] < row_counts[0],
        "-max_target_seqs 1 must drop a subject of this fixture: {row_counts:?}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91 (the query
// batch size, read where `CBlastInput` is made, before the formatter: tblastx_app.cpp:136-137)
// ```c++
//     char* batch_sz_str = getenv("BATCH_SIZE");
//     if (batch_sz_str) {
//         retval = NStr::StringToInt(batch_sz_str);
// ```
// A value that NCBI cannot convert ends the program with its `CStringException` (`CATCH_ALL`,
// exit 255) before anything is written (env parity probe h_probe of 2026-11-12, cases
// `BATCH_SIZE.*.tblastx-*`), also before the records without a title are searched.
#[test]
fn a_batch_size_that_ncbi_cannot_convert_stops_before_anything() {
    let inputs = Inputs::new();
    for (value, line, pos) in [("abc", 861, 0), (" ", 861, 0), ("12x", 870, 2)] {
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("tblastx")
            .arg("-query")
            .arg(&inputs.query.0)
            .arg("-subject")
            .arg(&inputs.subject.0)
            .args(["-outfmt", "0"])
            .env("BATCH_SIZE", value)
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI");
        assert_eq!(output.status.code(), Some(255), "{output:?}");
        assert!(output.stdout.is_empty(), "no prolog: {output:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            run_local_support::string_to_int_exception(value, line, pos)
        );
    }
    let sequence = fixture_sequence("LC738884.fasta");
    let empty_defline = TempFasta::new("tblastx_empty_defline.fna", &[("", &sequence[0..900])]);
    for (query, subject) in [
        (&empty_defline.0, &inputs.subject.0),
        (&inputs.query.0, &empty_defline.0),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("tblastx")
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .args(["-outfmt", "6"])
            .env("BATCH_SIZE", "abc")
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI");
        assert_eq!(output.status.code(), Some(255), "{output:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            run_local_support::string_to_int_exception("abc", 861, 0)
        );
    }
}

// NCBI reference: c++/src/algo/blast/api/local_blast.cpp:54-62,
// c++/src/algo/blast/api/split_query_aux_priv.cpp:53-60 and
// c++/src/app/blast/blast_app_util.cpp:206-210,732-737: NCBI converts CHUNK_SIZE and
// OVERLAP_CHUNK_SIZE in the first batch's query splitter and PRE_FETCH_SEQS_LIMIT before the
// first batch's report with NStr::StringToInt, whose `CStringException` ends the program (exit
// 255; outfmt 6 writes nothing before it), and searches each subject on its own with
// BL2SEQ_LEGACY, which LOSAT refuses (env parity probe h_probe of 2026-11-12).
#[test]
fn environment_that_ncbi_cannot_convert_or_reports_otherwise_stops_the_run() {
    let inputs = Inputs::new();
    let run = |variable: &str, value: &str| {
        Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("tblastx")
            .arg("-query")
            .arg(&inputs.query.0)
            .arg("-subject")
            .arg(&inputs.subject.0)
            .args(["-outfmt", "6"])
            .env(variable, value)
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI")
    };
    for (variable, value, line, pos) in [
        ("CHUNK_SIZE", "abc", 861, 0),
        ("OVERLAP_CHUNK_SIZE", "1e3", 870, 1),
        ("PRE_FETCH_SEQS_LIMIT", "", 861, 0),
    ] {
        let output = run(variable, value);
        assert_eq!(output.status.code(), Some(255), "{variable}: {output:?}");
        assert!(output.stdout.is_empty(), "{variable}: {output:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            run_local_support::string_to_int_exception(value, line, pos),
            "{variable}"
        );
    }
    let output = run("BL2SEQ_LEGACY", "0");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(
        output.stdout.is_empty()
            && stderr.contains("the environment variable BL2SEQ_LEGACY")
            && stderr.contains("not supported by LOSAT's TBLASTX"),
        "{stderr}"
    );
}

// A standard output that the caller opened on /dev/null, read and write (as Python's
// subprocess.DEVNULL and Node's `stdio: 'ignore'`), is written as any file (S08 audit (c),
// round 2, N2). A standard output closed at the start becomes the same /dev/null in Rust's
// runtime (`cli::report_standard_output`).
#[cfg(unix)]
#[test]
fn a_standard_output_on_dev_null_is_written() {
    let inputs = Inputs::new();
    for outfmt in ["0", "6", "7"] {
        let null = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open("/dev/null")
            .expect("open /dev/null");
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("tblastx")
            .arg("-query")
            .arg(&inputs.query.0)
            .arg("-subject")
            .arg(&inputs.subject.0)
            .args(["-outfmt", outfmt])
            .stdout(std::process::Stdio::from(null))
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI");
        assert_eq!(output.status.code(), Some(0), "outfmt {outfmt}: {output:?}");
        assert!(output.stderr.is_empty(), "outfmt {outfmt}: {output:?}");
    }
}

/// Runs the CLI with `query` and `subject` and returns its exit code, standard output and
/// standard error.
fn cli_run(
    query: &std::path::Path,
    subject: &std::path::Path,
    extra: &[&str],
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(extra)
        .env_remove("LOSAT_TIMING");
    for (variable, value) in env {
        command.env(variable, value);
    }
    let output = command.output().expect("run LOSAT CLI");
    (
        output.status.code(),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

// NCBI reference: c++/src/objtools/readers/fasta.cpp:375-384
// ```c
//         if (line.empty()) {
//             continue; // ignore lines containing only whitespace
//         }
//         c = line[0];
//
//         if (c == '!'  ||  c == '#' || c == ';') {
//             // no content, just a comment or blank line
//             continue;
// ```
// NCBI reads the subjects before it finds the query empty (tblastx_app.cpp:118-132), with
// its reader: after the lines it skips, a line that is not plausible data stops the
// reading (`BLAST query error`), and a file of such lines only has no subject.
#[test]
fn subjects_are_read_before_the_query_is_found_empty() {
    let empty = TempFasta::new("tblastx_no_query.fna", &[]);
    let subject = TempFasta::new("tblastx_skipped_lines.fna", &[]);
    let run = |text: &[u8]| {
        std::fs::write(&subject.0, text).expect("write subject");
        cli_run(&empty.0, &subject.0, &["-outfmt", "6"], &[])
    };
    let (code, stdout, stderr) = run(b"\n; comment\n>s1\nACGTAC%%%%\n");
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stdout.is_empty(), "{stderr}");
    assert!(
        stderr.starts_with("BLAST query error: CFastaReader: Near line 4, there's a line that doesn't look like plausible data"),
        "{stderr}"
    );
    let (code, stdout, stderr) = run(b"\n; comment\n>s1\nACGTACGTAC\n");
    assert_eq!(code, Some(0), "{stderr}");
    assert!(stdout.is_empty(), "{stderr}");
    assert_eq!(stderr, "Warning: [tblastx] Query is Empty!\n");
    let (code, stdout, stderr) = run(b"; only comment\n!x\n");
    assert_eq!(code, Some(3), "{stderr}");
    assert!(stdout.is_empty(), "{stderr}");
    assert_eq!(stderr, "BLAST engine error: Empty CBlastQueryVector\n");
}

// NCBI reference: c++/src/objtools/readers/fasta.cpp:919-935 (CFastaReader: a residue that is
// not an IUPAC nucleotide letter is removed with a warning).
// This is the command line's half: the warning for the invalid residue and the rows. The
// entry's half (NCBI's check of the hit saving options, blast_options.c:1518-1523, before
// `Query is Empty!`) is in the crate, `LOSAT/src/api/local_blast_tests/tblastx.rs`.
#[test]
fn residues_are_read_as_ncbis_reader_reads_them() {
    let inputs = Inputs::new();
    let sequence = fixture_sequence("LC738884.fasta");
    let mut with_digit = sequence[1_000..3_000].to_vec();
    with_digit[700] = b'7';
    let digit = TempFasta::new("tblastx_digit.fna", &[("digit", &with_digit)]);
    let (code, stdout, stderr) = cli_run(&digit.0, &inputs.subject.0, &["-outfmt", "6"], &[]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        "FASTA-Reader: Ignoring invalid residues at position(s): On line 13: 41\n"
    );
    assert!(stdout.starts_with(b"digit\t"), "{stderr}");
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-788
// ```c
//         catch(CBlastException & e ) {
//         	// Skip bad subject sequence
//         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
//         		seqblk_vec->push_back(subj);
//         ...
//         		warning += "Subject sequence contains no data";
//         		ERR_POST(Warning << warning);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-636
// ```c
//         } catch (const CException& e) {
//             // FIXME: is index this the right value for the 2nd arg? Also, how
//             // to determine whether the message should contain a warning or
//             // error?
//             CRef<CSearchMessage> m
// ```
// A record without residues stays in its role: a subject gets its warning before the
// search, counts in the database statistics and is never searched; a query gets its
// warning before its report, which has no hits and no valid context, and is counted.
#[test]
fn records_without_residues_are_kept_with_ncbis_warnings() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "tblastx_empty_records_query.fna",
            &[
                ("q1", &genome[0..2_000]),
                ("q2 empty two", b""),
                ("q3", &genome[10_000..12_000]),
            ],
        ),
        subject: TempFasta::new(
            "tblastx_empty_records_subject.fna",
            &[
                ("s1", &genome[0..5_000]),
                ("s2 no letters", b""),
                ("s3", &genome[9_000..13_000]),
            ],
        ),
    };
    let run = inputs.collect(&FORMATS, &[]);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [tblastx] Subject_2 s2 no letters: Subject sequence contains no data\n\
         Warning: [tblastx] Query_2 q2 empty two: Sequence contains no data \n"
    );
    let pairwise = String::from_utf8_lossy(&run.reports[0].bytes);
    assert!(
        pairwise.contains("3 sequences; 9,000 total letters"),
        "{pairwise}"
    );
    assert!(
        pairwise.contains("Query= q2 empty two\n\nLength=0\n\n\n***** No hits found *****\n"),
        "{pairwise}"
    );
    let tabular = String::from_utf8_lossy(&run.reports[1].bytes);
    assert!(
        tabular.lines().any(|row| row.starts_with("q1\ts1\t"))
            && tabular.lines().any(|row| row.starts_with("q3\ts3\t")),
        "{tabular}"
    );
    assert!(
        !tabular.contains("q2") && !tabular.contains("s2"),
        "{tabular}"
    );
    let commented = String::from_utf8_lossy(&run.reports[2].bytes);
    assert!(
        commented.contains("# Query: q2 empty two\n# Database: ")
            && commented.contains("# BLAST processed 3 queries\n"),
        "{commented}"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eAsnText
    //       ...
    //       || m_FormatType == CFormattingArgs::eJsonSeqalign
    //       || m_FormatType == CFormattingArgs::eSAM)
    //     {
    //         x_PrintStructuredReport(results, queries);
    //         return;
    //     }
    //     ...
    //     if (results.HasWarnings()) {
    //         ERR_POST(Warning << results.GetWarningStrings());
    //     }
    // ```
    // The warnings go with the first format (outfmt 0 here); outfmt 5, 8, 9, 11, 12, 15 and 16 write
    // no query result warnings, so their CLI runs equal the library's run with that format alone.
    for (index, outfmt) in FORMATS.into_iter().enumerate() {
        let (cli_output, cli_stderr) = inputs.cli(outfmt, &[]);
        assert!(run.reports[index].bytes == cli_output, "outfmt {outfmt}");
        if matches!(outfmt, "5" | "8" | "9" | "11" | "12" | "15" | "16") {
            let alone = inputs.collect(&[outfmt], &[]);
            assert_eq!(
                String::from_utf8_lossy(&alone.diagnostics),
                "Warning: [tblastx] Subject_2 s2 no letters: Subject sequence contains no data\n"
            );
            assert_eq!(alone.diagnostics, cli_stderr, "outfmt {outfmt}");
        } else {
            assert_eq!(run.diagnostics, cli_stderr, "outfmt {outfmt}");
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:146-152
// ```c
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
// ```
// What ends the reading of the queries comes in the batch where NCBI reads it: a reader
// error in a later batch after the reports of the batches before (without the epilog), in
// the first batch before any report. Blank and comment lines after the last record end
// that record, and the run ends as usual.
#[test]
fn what_ends_the_query_reading_comes_with_its_batch() {
    let genome = fixture_sequence("LC738884.fasta");
    let subject = TempFasta::new("tblastx_batch_subject.fna", &[("s1", &genome[0..5_000])]);
    let query = TempFasta::new("tblastx_batch_query.fna", &[]);
    let mut text = b">q1\n".to_vec();
    text.extend_from_slice(&genome[0..600]);
    text.extend_from_slice(b"\n>q2\n");
    text.extend_from_slice(&genome[1_000..1_600]);
    text.extend_from_slice(b"\n");
    let mut bad = text.clone();
    bad.extend_from_slice(b">q3\nACGT%%%%\n");
    std::fs::write(&query.0, &bad).expect("write query");
    let (code, stdout, stderr) = cli_run(
        &query.0,
        &subject.0,
        &["-outfmt", "7"],
        &[("BATCH_SIZE", "600")],
    );
    assert_eq!(code, Some(1), "{stderr}");
    let report = String::from_utf8_lossy(&stdout);
    assert!(
        report.contains("# Query: q1\n")
            && report.contains("# Query: q2\n")
            && !report.contains("# BLAST processed"),
        "{report}"
    );
    assert!(
        stderr.starts_with("BLAST query error: CFastaReader: "),
        "{stderr}"
    );
    let (code, stdout, stderr) = cli_run(&query.0, &subject.0, &["-outfmt", "7"], &[]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stdout.is_empty(), "{stderr}");
    assert!(
        stderr.starts_with("BLAST query error: CFastaReader: "),
        "{stderr}"
    );
    let mut blank = text.clone();
    blank.extend_from_slice(b"\n; comment\n\n");
    std::fs::write(&query.0, &blank).expect("write query");
    let (code, stdout, stderr) = cli_run(
        &query.0,
        &subject.0,
        &["-outfmt", "7"],
        &[("BATCH_SIZE", "600")],
    );
    assert_eq!(code, Some(0), "{stderr}");
    assert!(String::from_utf8_lossy(&stdout).ends_with("# BLAST processed 2 queries\n"));
}
