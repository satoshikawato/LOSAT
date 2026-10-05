#![allow(warnings, clippy::all)]

//! The shared TBLASTX entry `run_local` writes several output formats from one search.
//! Every format must be byte-identical to the CLI run that requests only that format, the
//! warnings must equal that run's standard error, the hit records must be the final HSP
//! list, and the observer must report the exact row and pairwise section of every HSP.
//! With `-num_threads 2` the native search reduces its results in a collector thread; the
//! result must not change. (Identity of outfmt 0 and 7 with NCBI is checked by
//! docs/evidence/losat_web_e2a/check_losat.py over the frozen fixtures; identity with the
//! previous release over the regression manifests by
//! docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::process::Command;

use run_local_support::{
    assert_observer_ranges, fixture_sequence, read_records, run_formats, Run, TempFasta,
};
use LOSAT::algorithm::tblastx::TblastxArgs;
use LOSAT::api::local_blast::{run_local_tblastx, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 3] = ["0", "6", "7"];

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

    fn args(&self, extra: &[&str]) -> TblastxArgs {
        let query = self.query.0.display().to_string();
        let subject = self.subject.0.display().to_string();
        let mut argv = vec![
            "LOSAT", "tblastx", "-query", &query, "-subject", &subject, "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv).expect("tblastx argv");
        match cli.command {
            Commands::Tblastx(args) => args,
            _ => unreachable!("tblastx subcommand"),
        }
    }

    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        let queries = read_records(&self.query.0);
        let subjects = read_records(&self.subject.0);
        run_formats(formats, observe, |outputs| {
            run_local_tblastx(self.args(extra), &queries, &subjects, outputs)
        })
    }

    /// The CLI output file and standard error of one single-format run.
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = self.query.0.with_extension(format!("cli{outfmt}.out"));
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
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
    let serial = inputs.run(&FORMATS, &[], true);
    let mut row_counts = Vec::new();
    for extra in [
        &[][..],
        &["-num_threads", "2"][..],
        &["-max_target_seqs", "1"][..],
    ] {
        let together = inputs.run(&FORMATS, extra, true);
        assert_eq!(
            together.hits.len(),
            together.ranges[1].len(),
            "one hit record per outfmt 6 row {extra:?}"
        );
        row_counts.push(together.ranges[1].len());
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.run(&[outfmt], extra, false);
            assert!(
                together.outputs[index] == alone.outputs[0],
                "run_local outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert!(
                together.outputs[index] == cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            assert_eq!(
                together.diagnostics, cli_stderr,
                "warnings of three formats equal those of the single-format CLI run {extra:?}"
            );
        }
        if extra.first() == Some(&"-num_threads") {
            assert!(
                together.outputs == serial.outputs,
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

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     // Add tab in front of field, except for the first field.
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
#[test]
fn observer_ranges_are_exact_rows_and_sections_of_the_same_hsp() {
    let inputs = Inputs::new();
    for extra in [
        &[][..],
        &["-num_threads", "2"][..],
        &["-max_target_seqs", "1"][..],
    ] {
        let result = inputs.run(&FORMATS, extra, true);
        assert_observer_ranges(&result, [0, 1, 2], &format!("tblastx {extra:?}"));
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
// ```c
// arg_desc.AddDefaultKey(kArgOutputFormat, "format",
//                        kOutputFormatDescription,
//                        CArgDescriptions::eString,
//                        NStr::IntToString(dft_outfmt));
// ```
#[test]
fn unsupported_formats_fail_before_searching() {
    let inputs = Inputs::new();
    let queries = read_records(&inputs.query.0);
    let subjects = read_records(&inputs.subject.0);
    for outfmt in ["5", "6 qseqid", "abc"] {
        let (mut valid, mut invalid, mut diagnostics) = (Vec::new(), Vec::new(), Vec::new());
        let mut outputs = ReportOutputs {
            formats: vec![
                FormatOutput {
                    outfmt: "6",
                    sink: OutputSink::Writer(&mut valid),
                },
                FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(&mut invalid),
                },
            ],
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        let result = run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:86-88
// ```c
//     char* batch_sz_str = getenv("BATCH_SIZE");
//     if (batch_sz_str) {
//         retval = NStr::StringToInt(batch_sz_str);
// ```
// NCBI stops with a CStringException (exit 255, with the path of its build in the message);
// LOSAT rejects the value, before the prolog as NCBI.
#[test]
fn a_batch_size_that_is_not_an_integer_is_rejected() {
    let inputs = Inputs::new();
    for value in ["abc", " ", "12x"] {
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
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
        assert_eq!(output.status.code(), Some(1), "{output:?}");
        assert!(output.stdout.is_empty(), "no prolog: {output:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!(
                "Error: the BATCH_SIZE value '{value}' is not an integer; NCBI stops with a \
                 CStringException for it, which is not supported by LOSAT's TBLASTX\n"
            )
        );
    }
    // NCBI reads the batch size before the queries (tblastx_app.cpp:136-137) and reads
    // the deflines that LOSAT rejects without a message: the batch size fails first.
    let sequence = fixture_sequence("LC738884.fasta");
    let empty_defline = TempFasta::new("tblastx_empty_defline.fna", &[("", &sequence[0..900])]);
    for (query, subject) in [
        (&empty_defline.0, &inputs.subject.0),
        (&inputs.query.0, &empty_defline.0),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
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
        assert_eq!(output.status.code(), Some(1), "{output:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).starts_with("Error: the BATCH_SIZE value"),
            "{output:?}"
        );
    }
}

// NCBI reference: c++/src/algo/blast/api/local_blast.cpp:54-62,
// c++/src/algo/blast/api/split_query_aux_priv.cpp:53-60 and
// c++/src/app/blast/blast_app_util.cpp:206-210,732-737: NCBI converts CHUNK_SIZE,
// OVERLAP_CHUNK_SIZE and PRE_FETCH_SEQS_LIMIT with NStr::StringToInt (its CStringException
// names its build's source files) and searches each subject on its own with BL2SEQ_LEGACY.
#[test]
fn environment_that_ncbi_cannot_convert_or_reports_otherwise_is_rejected() {
    let inputs = Inputs::new();
    for (variable, value, text) in [
        (
            "CHUNK_SIZE",
            "abc",
            "the environment variable CHUNK_SIZE has the value \"abc\"",
        ),
        (
            "OVERLAP_CHUNK_SIZE",
            "1e3",
            "the environment variable OVERLAP_CHUNK_SIZE has the value \"1e3\"",
        ),
        (
            "PRE_FETCH_SEQS_LIMIT",
            "",
            "the environment variable PRE_FETCH_SEQS_LIMIT has the value \"\"",
        ),
        (
            "BL2SEQ_LEGACY",
            "0",
            "the environment variable BL2SEQ_LEGACY",
        ),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
            .arg("tblastx")
            .arg("-query")
            .arg(&inputs.query.0)
            .arg("-subject")
            .arg(&inputs.subject.0)
            .args(["-outfmt", "6"])
            .env(variable, value)
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{variable}: {output:?}");
        assert!(
            output.stdout.is_empty()
                && stderr.contains(text)
                && stderr.contains("not supported by LOSAT's TBLASTX"),
            "{variable}: {stderr}"
        );
    }
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
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
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
// NCBI reads the subjects before it finds the query empty (tblastx_app.cpp:118-132): after
// the lines it skips, it reads the records as those of any file (S08 audit round 3, (a)
// N1-ii and N1-iii, (d) L1). A subject that `bio` cannot read because of such lines has its
// residues checked when it is read, and a file of such lines only has no subject.
#[test]
fn subjects_after_skipped_lines_are_read_before_the_query_is_found_empty() {
    let empty = TempFasta::new("tblastx_no_query.fna", &[]);
    let subject = TempFasta::new("tblastx_skipped_lines.fna", &[]);
    let run = |text: &[u8]| {
        std::fs::write(&subject.0, text).expect("write subject");
        Command::new(env!("CARGO_BIN_EXE_LOSAT"))
            .arg("tblastx")
            .arg("-query")
            .arg(&empty.0)
            .arg("-subject")
            .arg(&subject.0)
            .args(["-outfmt", "6"])
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI")
    };
    let output = run(b"\n; comment\n>s1\nACGTAC%%%%\n");
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("subject record 1 (s1) has '%' at residue 7")
            && stderr.contains("LOSAT's TBLASTX"),
        "{stderr}"
    );
    let output = run(b"; only comment\n!x\n");
    assert_eq!(output.status.code(), Some(3), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "BLAST engine error: Empty CBlastQueryVector\n"
    );
}

// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312 (x_CleanAndCompress) and
// 4066 (`NStr::HtmlDecode`): NCBI reads past the end of the title `, ,` (and crashes when
// such a subject has hits) and decodes `&amp;` in outfmt 0; the tabular formats print
// only the ids. LOSAT writes the title `, ` stopped at the end of the string (approved
// exception 2 of PD-LOSAT-NCBI-DEFECTS) and rejects the decoded title.
#[test]
fn outfmt0_titles_that_ncbi_reads_past_or_decodes() {
    let inputs = Inputs::new();
    let queries = read_records(&inputs.query.0);
    let sequence = fixture_sequence("LC738884.fasta");
    let unknown = vec![b'N'; 900];
    for (defline, reason) in [(", ,", ""), ("s &amp; t", "HTML character reference")] {
        let subject = TempFasta::new(
            "tblastx_title_subject.fna",
            &[(defline, &sequence[1_000..6_000])],
        );
        let subjects = read_records(&subject.0);
        for outfmt in FORMATS {
            let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut report), &mut diagnostics);
            let result = run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs);
            drop(outputs);
            if outfmt == "0" && reason.is_empty() {
                result.unwrap_or_else(|error| panic!("outfmt 0 {defline:?}: {error}"));
                let report = String::from_utf8(report).expect("UTF-8 report");
                assert!(report.contains("\n> , \n"), "{report}");
            } else if outfmt == "0" {
                let error = result
                    .expect_err("outfmt 0 must reject the title")
                    .to_string();
                assert!(
                    error.contains(reason) && error.contains("not supported by LOSAT's TBLASTX"),
                    "{error}"
                );
                // The rejection comes after the prolog, before the first query's report.
                let report = String::from_utf8(report).expect("UTF-8 report");
                assert!(report.starts_with("TBLASTX 2.17.0+"), "{report}");
                assert!(!report.contains("Query=") && diagnostics.is_empty());
            } else {
                result.unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            }
        }
        // NCBI makes the titles of the subjects that the report shows only: a subject
        // without hits keeps its title out of the report.
        let subject = TempFasta::new(
            "tblastx_title_hitless.fna",
            &[("hit", &sequence[1_000..6_000]), (defline, &unknown)],
        );
        let subjects = read_records(&subject.0);
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs =
            ReportOutputs::single("0", OutputSink::Writer(&mut report), &mut diagnostics);
        run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs)
            .unwrap_or_else(|error| panic!("{defline:?} without hits: {error}"));
    }
    // Titles that NCBI's HtmlDecode leaves as they are (not a name of its table, a final
    // `;` trimmed before the decoding) are written.
    for defline in ["R&D; x", "a&foo;b", "a&amp;"] {
        let subject = TempFasta::new(
            "tblastx_title_kept.fna",
            &[(defline, &sequence[1_000..6_000])],
        );
        let subjects = read_records(&subject.0);
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs =
            ReportOutputs::single("0", OutputSink::Writer(&mut report), &mut diagnostics);
        run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs)
            .unwrap_or_else(|error| panic!("{defline:?}: {error}"));
    }
}

// NCBI reference: c++/src/objtools/readers/fasta.cpp:919-935 (CFastaReader: a residue that is
// not an IUPAC nucleotide letter is removed with a warning) and the "Sequence contains no
// data" warning of a record without residues; LOSAT rejects both, as BLASTN does.
#[test]
fn inputs_ncbi_reads_differently_are_rejected() {
    let inputs = Inputs::new();
    let queries = read_records(&inputs.query.0);
    let subjects = read_records(&inputs.subject.0);
    let sequence = fixture_sequence("LC738884.fasta");
    let mut with_digit = sequence[1_000..3_000].to_vec();
    with_digit[700] = b'7';
    let digit = TempFasta::new("tblastx_digit.fna", &[("digit", &with_digit)]);
    let empty = TempFasta::new(
        "tblastx_empty_record.fna",
        &[("empty", &[][..]), ("full", &sequence[0..900])],
    );
    let run = |extra: &[&str], queries: &[_], subjects: &[_]| {
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs =
            ReportOutputs::single("6", OutputSink::Writer(&mut report), &mut diagnostics);
        let result = run_local_tblastx(inputs.args(extra), queries, subjects, &mut outputs);
        drop(outputs);
        assert!(report.is_empty() && diagnostics.is_empty());
        result.expect_err("rejected").to_string()
    };
    let error = run(&[], &read_records(&digit.0), &subjects);
    assert!(
        error.contains("query record 1 (digit) has '7' at residue 701")
            && error.contains("LOSAT's TBLASTX"),
        "{error}"
    );
    let error = run(&[], &queries, &read_records(&empty.0));
    assert!(
        error.contains("subject record 1 (empty) has no residues")
            && error.contains("LOSAT's TBLASTX"),
        "{error}"
    );
    // NCBI's check of the hit saving options (blast_options.c:1518-1523), before
    // `Query is Empty!`.
    for evalue in ["0", "1e-400"] {
        for queries in [&queries[..], &[]] {
            assert_eq!(
                run(&["-evalue", evalue], queries, &subjects),
                "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n"
            );
        }
    }
}
