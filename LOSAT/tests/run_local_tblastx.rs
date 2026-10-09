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
    assert_observer_ranges, fixture_sequence, reader_records, run_formats, Run, TempFasta,
};
use LOSAT::algorithm::tblastx::TblastxArgs;
use LOSAT::api::local_blast::{run_local_tblastx, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
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
        let queries = queries(&self.query.0);
        let subjects = subjects(&self.subject.0);
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

/// The query records of a FASTA file as NCBI's reader reads them for TBLASTX.
fn queries(path: &std::path::Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::query("TBLASTX", false, true))
}

/// The subject records of a FASTA file as NCBI's reader reads them for TBLASTX.
fn subjects(path: &std::path::Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::subject("TBLASTX", false, true))
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
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
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
    // NCBI reads the batch size before the queries are searched (tblastx_app.cpp:136-137):
    // it fails before the records without a title are searched.
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

/// Runs the CLI with `query` and `subject` and returns its exit code, standard output and
/// standard error.
fn cli_run(
    query: &std::path::Path,
    subject: &std::path::Path,
    extra: &[&str],
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

/// Runs `run_local` with one writer per format and returns the outputs, the diagnostics
/// and the result.
fn run_local_outputs(
    inputs: &Inputs,
    formats: &[&str],
    extra: &[&str],
) -> (Vec<Vec<u8>>, Vec<u8>, anyhow::Result<()>) {
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
    let mut sinks: Vec<Vec<u8>> = vec![Vec::new(); formats.len()];
    let mut diagnostics = Vec::new();
    let result = {
        let mut outputs = ReportOutputs {
            formats: formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(outfmt, sink)| FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(sink),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        run_local_tblastx(inputs.args(extra), &queries, &subjects, &mut outputs)
    };
    (sinks, diagnostics, result)
}

/// The exit code and the message of NCBI's error that a run ended with.
fn native_error(result: anyhow::Result<()>) -> (i32, String) {
    let error = result.expect_err("the run fails");
    let native = error
        .downcast_ref::<LOSAT::cli::NativeError>()
        .unwrap_or_else(|| panic!("NCBI's error: {error:#}"));
    (native.exit, native.message.clone())
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

// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312 (x_CleanAndCompress) and
// 4066 (`NStr::HtmlDecode`): NCBI reads past the end of the title `, ,` (and crashes when
// such a subject has hits) and decodes `&amp;` in outfmt 0; the tabular formats print
// only the ids. LOSAT writes the title `, ` stopped at the end of the string (approved
// exception 2 of PD-LOSAT-NCBI-DEFECTS) and the decoded title (`report/defline.rs`).
#[test]
fn outfmt0_titles_are_written_as_ncbi_makes_them() {
    let inputs = Inputs::new();
    let queries = queries(&inputs.query.0);
    let sequence = fixture_sequence("LC738884.fasta");
    for (defline, heading) in [
        (", ,", "\n> , \n"),
        ("s &amp; t", "\n> s & t\n"),
        ("R&D; x", "\n> R&D; x\n"),
        ("a&foo;b", "\n> a&foo;b\n"),
        ("a&amp;", "\n> a&amp\n"),
    ] {
        let subject = TempFasta::new(
            "tblastx_title_subject.fna",
            &[(defline, &sequence[1_000..6_000])],
        );
        let subjects = subjects(&subject.0);
        for outfmt in FORMATS {
            let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut report), &mut diagnostics);
            run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs)
                .unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            drop(outputs);
            if outfmt == "0" {
                let report = String::from_utf8(report).expect("UTF-8 report");
                assert!(report.contains(heading), "{defline:?}: {report}");
            }
        }
    }
}

// NCBI reference: c++/src/objtools/readers/fasta.cpp:919-935 (CFastaReader: a residue that is
// not an IUPAC nucleotide letter is removed with a warning). NCBI's check of the hit saving
// options (blast_options.c:1518-1523) comes before `Query is Empty!`.
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
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
    let run = |extra: &[&str], queries: &[FastaRecord]| {
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs =
            ReportOutputs::single("6", OutputSink::Writer(&mut report), &mut diagnostics);
        let result = run_local_tblastx(inputs.args(extra), queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(report.is_empty() && diagnostics.is_empty());
        result.expect_err("rejected").to_string()
    };
    for evalue in ["0", "1e-400"] {
        for queries in [&queries[..], &[]] {
            assert_eq!(
                run(&["-evalue", evalue], queries),
                "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n"
            );
        }
    }
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
    let run = inputs.run(&FORMATS, &[], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [tblastx] Subject_2 s2 no letters: Subject sequence contains no data\n\
         Warning: [tblastx] Query_2 q2 empty two: Sequence contains no data \n"
    );
    let pairwise = String::from_utf8_lossy(&run.outputs[0]);
    assert!(
        pairwise.contains("3 sequences; 9,000 total letters"),
        "{pairwise}"
    );
    assert!(
        pairwise.contains("Query= q2 empty two\n\nLength=0\n\n\n***** No hits found *****\n"),
        "{pairwise}"
    );
    let tabular = String::from_utf8_lossy(&run.outputs[1]);
    assert!(
        tabular.lines().any(|row| row.starts_with("q1\ts1\t"))
            && tabular.lines().any(|row| row.starts_with("q3\ts3\t")),
        "{tabular}"
    );
    assert!(
        !tabular.contains("q2") && !tabular.contains("s2"),
        "{tabular}"
    );
    let commented = String::from_utf8_lossy(&run.outputs[2]);
    assert!(
        commented.contains("# Query: q2 empty two\n# Database: ")
            && commented.contains("# BLAST processed 3 queries\n"),
        "{commented}"
    );
    for (index, outfmt) in FORMATS.into_iter().enumerate() {
        let (cli_output, cli_stderr) = inputs.cli(outfmt, &[]);
        assert!(run.outputs[index] == cli_output, "outfmt {outfmt}");
        assert_eq!(run.diagnostics, cli_stderr, "outfmt {outfmt}");
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:649-652
// ```c
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// A batch whose queries all have fewer letters than a codon stops NCBI after the prolog
// when one of them has no letters, with one message per query without letters; a batch of
// queries of one or two letters only is not searched.
#[test]
fn a_batch_without_valid_contexts_and_with_queries_without_letters_stops() {
    let genome = fixture_sequence("LC738884.fasta");
    let subject = TempFasta::new("tblastx_no_data_subject.fna", &[("s1", &genome[0..5_000])]);
    for (records, message) in [
        (
            vec![("q1", &b""[..]), ("q2 two", &b""[..])],
            "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n",
        ),
        (
            vec![("q1", &b"AC"[..]), ("q2 two", &b""[..])],
            "BLAST engine error: Warning: Sequence contains no data \n",
        ),
    ] {
        let inputs = Inputs {
            query: TempFasta::new("tblastx_no_data_query.fna", &records),
            subject: TempFasta::new("tblastx_no_data_subject.fna", &[("s1", &genome[0..5_000])]),
        };
        let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["0", "6"], &[]);
        assert_eq!(native_error(result), (3, message.to_string()));
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        let prolog = String::from_utf8_lossy(&outputs[0]);
        assert!(
            prolog.starts_with("TBLASTX 2.17.0+") && !prolog.contains("Query="),
            "{prolog}"
        );
        assert!(outputs[1].is_empty());
    }
    let inputs = Inputs {
        query: TempFasta::new("tblastx_short_query.fna", &[("q1", b"AC"), ("q2", b"A")]),
        subject,
    };
    let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["7"], &[]);
    result.expect("a batch of short queries is reported unsearched");
    let warnings = String::from_utf8_lossy(&diagnostics);
    assert_eq!(
        warnings
            .matches("Could not calculate ungapped Karlin-Altschul parameters")
            .count(),
        2,
        "{warnings}"
    );
    let report = String::from_utf8_lossy(&outputs[0]);
    assert!(
        report.ends_with("# BLAST processed 2 queries\n"),
        "{report}"
    );
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-980
// ```c
//    if (sbp->gbp) {
//        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//        if (Blast_SubjectIsTranslated(program_number)) {
//            min_subject_length/=3;
//        }
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
//
//    if(min_subject_length <=0) {
// 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
//    }
// ```
// Subjects with fewer letters than subjects (every subject without residues; AUTHORITY
// §E3, BI-44) stop NCBI when it searches the first batch with a valid context, after the
// subjects' warnings and the prolog.
#[test]
fn subjects_with_fewer_letters_than_subjects_stop_the_search() {
    let genome = fixture_sequence("LC738884.fasta");
    for (subjects, warnings) in [
        (
            vec![("s1", &b"A"[..]), ("s2", &b""[..])],
            "Warning: [tblastx] Subject_2 s2: Subject sequence contains no data\n",
        ),
        (
            vec![("s1 one", &b""[..]), ("s2", &b""[..])],
            "Warning: [tblastx] Subject_1 s1 one: Subject sequence contains no data\n\
             Warning: [tblastx] Subject_2 s2: Subject sequence contains no data\n",
        ),
    ] {
        let inputs = Inputs {
            query: TempFasta::new(
                "tblastx_short_subjects_query.fna",
                &[("q1", &genome[0..500])],
            ),
            subject: TempFasta::new("tblastx_short_subjects.fna", &subjects),
        };
        let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["0", "7"], &[]);
        assert_eq!(
            native_error(result),
            (
                3,
                "BLAST engine error: The average subject length is too short\n".to_string()
            )
        );
        assert_eq!(String::from_utf8_lossy(&diagnostics), warnings);
        assert!(String::from_utf8_lossy(&outputs[0]).starts_with("TBLASTX 2.17.0+"));
        assert!(outputs[1].is_empty());
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-460
// ```c++
//     // set sequence range
//     retval->SetInt().SetFrom(from);
//     retval->SetInt().SetTo((to > 0 && to < seqlen) ? to : (seqlen-1));
// ```
// A range that starts just past a record's end gives an interval without letters, which
// NCBI sets up as a record without residues.
#[test]
fn intervals_without_letters_are_records_without_residues() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "tblastx_empty_interval_query.fna",
            &[("q1", &genome[0..2_000]), ("q2 short", &genome[0..300])],
        ),
        subject: TempFasta::new(
            "tblastx_empty_interval_subject.fna",
            &[("s1", &genome[0..5_000]), ("s2 short", &genome[0..300])],
        ),
    };
    let run = inputs.run(&["6"], &["-subject_loc", "301-1000"], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [tblastx] Subject_2 s2 short: Subject sequence contains no data\n"
    );
    let rows = String::from_utf8_lossy(&run.outputs[0]);
    assert!(rows.lines().all(|row| row.contains("\ts1\t")), "{rows}");
    assert!(!rows.is_empty());
    let run = inputs.run(&["7"], &["-query_loc", "301-1000"], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [tblastx] Query_2 q2 short: Sequence contains no data \n"
    );
    let report = String::from_utf8_lossy(&run.outputs[0]);
    assert!(
        report.contains("# Query: q2 short\n") && report.contains("# BLAST processed 2 queries\n"),
        "{report}"
    );
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
