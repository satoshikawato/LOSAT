//! The shared TBLASTN entry `run_local` in process: what the entry itself does, not what the
//! command line writes. The saved NCBI BLAST+ 2.17.0 bytes of three searches, the observer's
//! ranges, the entry's refusals and its sinks are checked here. The tests that compare the
//! entry's reports with the `losat` binary go through `losat::api`
//! (`LOSAT/tests/run_local_tblastn.rs`).
//! (Identity with the previous release over the Stage G matrix is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.)

use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use super::support::{
    assert_observer_ranges, assert_tabular_ranges, reader_records, run_formats, temp_path, Run,
};
use crate::algorithm::tblastn::TblastnArgs;
use crate::api::local_blast::{
    run_local_database_tblastn, run_local_tblastn, FormatOutput, OutputSink, ReportOutputs,
};
use crate::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::{try_parse_from, Cli, Commands};

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
/// The formats of the saved NCBI outputs (`docs/evidence/tlosan_stage_g/`).
const SAVED_FORMATS: [&str; 3] = ["0", "6", "7"];

fn repository(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join(path)
}

/// Parses the TBLASTN options. `query` and `subject` are display names only.
fn tblastn_args(query: &str, subject: &str, extra: &[&str]) -> TblastnArgs {
    let mut argv = vec!["losat", "tblastn", "-query", query, "-subject", subject];
    argv.extend_from_slice(extra);
    let cli: Cli = try_parse_from(argv).expect("tblastn argv");
    match cli.command {
        Commands::Tblastn(args) => args,
        _ => unreachable!("tblastn subcommand"),
    }
}

/// The queries as NCBI's protein reader reads them (the CLI's `fasta_reader`).
fn queries(path: &Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::query("TBLASTN", true, true))
}

/// The subjects as NCBI's nucleotide reader reads them (the CLI's `fasta_reader`).
fn subjects(path: &Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::subject("TBLASTN", false, true))
}

fn run(
    query: &Path,
    subject: &Path,
    names: (&str, &str),
    formats: &[&str],
    extra: &[&str],
    observe: bool,
) -> Run {
    let queries = queries(query);
    let subjects = subjects(subject);
    run_formats(formats, observe, |outputs| {
        run_local_tblastn(
            tblastn_args(names.0, names.1, extra),
            &queries,
            &subjects,
            outputs,
        )
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:662-673 (`ThrowIfInvalidTask`,
// quoted at `TaskParser::check`). The command line's parser refuses a task that is not TBLASTN's;
// a caller of `run_local` (or `run_local_database`, before it opens the database) that sets one
// itself gets NCBI's error for a name NCBI does not know, LOSAT's for another program's task,
// instead of a search with the default task.
#[test]
fn a_task_that_is_not_tblastns_is_refused() {
    let query = repository("LOSAT/tests/fasta/CoBV.faa");
    let subject = repository("LOSAT/tests/fasta/LC738874.fasta");
    for (task, database, message) in [
        (
            "tblastx",
            false,
            "the task tblastx is not supported by LOSAT's TBLASTN (use tblastn or tblastn-fast)",
        ),
        ("TBLASTZ", false, "'tblastz' is not a supported task"),
        ("tblastz", true, "'tblastz' is not a supported task"),
    ] {
        let mut args = tblastn_args("q", "s", &[]);
        args.task = task.into();
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        let result = {
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    sink: OutputSink::Writer(&mut sink),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            if database {
                args.db = Some("no_such_database".into());
                run_local_database_tblastn(
                    args,
                    &queries(&query),
                    ApplicationSettings::default(),
                    &mut outputs,
                )
            } else {
                run_local_tblastn(args, &queries(&query), &subjects(&subject), &mut outputs)
            }
        };
        let error = format!("{:#}", result.expect_err(message));
        assert!(error.contains(message), "{task}: {error}");
        assert!(sink.is_empty() && diagnostics.is_empty());
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
// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1443-1451
// ```c
// if (results.HasWarnings()) ERR_POST(Warning << results.GetWarningStrings());
// ```
// The saved NCBI outputs name the inputs by the paths NCBI read, so those paths are
// passed as the display names.
#[test]
fn formats_of_one_search_equal_the_ncbi_bytes_and_warn_once() {
    // Display names only (the library API reads `subject`, never these paths): the prefix is
    // the one that LOSAT/tests/frozen_paths.py names FROZEN_PREFIX for the tools' comparisons.
    let code1 = "/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna";
    let cases = [
        (
            "docs/evidence/tlosan_stage_g/batch_boundary/valid_then_skipped",
            "/tmp/tlosan-stageg-batch-boundary-20260926/valid_then_skipped.faa",
            "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna",
            code1,
            ".ncbi",
            Some(".fmt0.ncbi.stderr"),
        ),
        (
            "docs/evidence/tlosan_stage_g/batch_boundary/skipped_then_valid",
            "/tmp/tlosan-stageg-batch-boundary-20260926/skipped_then_valid.faa",
            "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna",
            code1,
            ".ncbi",
            Some(".fmt0.ncbi.stderr"),
        ),
        (
            "docs/evidence/tlosan_stage_g/native_real/av_code1",
            "/tmp/tlosan-stageg-final-real-nohit-20260926/av_code1.faa",
            "LOSAT/tests/fasta/AvCLPV.fasta",
            // A display name only, as `code1` above.
            "/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests/fasta/AvCLPV.fasta",
            ".oracle",
            None,
        ),
    ];
    // (saved output stem, NCBI query path, subject, NCBI subject path, output suffix,
    //  suffix of the saved NCBI standard error, if NCBI wrote any)
    for (stem, query_name, subject, subject_name, suffix, stderr_suffix) in cases {
        let query = repository(&format!("{stem}.faa"));
        let together = run(
            &query,
            &repository(subject),
            (query_name, subject_name),
            &SAVED_FORMATS,
            &[],
            true,
        );
        for (index, outfmt) in SAVED_FORMATS.iter().enumerate() {
            let expected = std::fs::read(repository(&format!("{stem}.fmt{outfmt}{suffix}")))
                .expect("saved NCBI output");
            assert!(
                together.outputs[index] == expected,
                "{stem} outfmt {outfmt} differs from the saved NCBI bytes"
            );
        }
        let expected_stderr = stderr_suffix.map_or_else(Vec::new, |stderr_suffix| {
            std::fs::read(repository(&format!("{stem}{stderr_suffix}"))).expect("saved stderr")
        });
        assert_eq!(
            String::from_utf8_lossy(&together.diagnostics),
            String::from_utf8_lossy(&expected_stderr),
            "{stem}: the warnings are written once for three formats"
        );
    }
}

/// Temporary FASTA inputs: several queries (one of them invalid) against several
/// subjects, some of which hold more than one HSP.
struct Inputs {
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let concatenate = |name: &str, parts: &[&str]| {
            let path = temp_path(&format!("tblastn_{name}"));
            let mut bytes = Vec::new();
            for part in parts {
                bytes.extend(std::fs::read(repository(part)).expect("fixture"));
            }
            std::fs::write(&path, bytes).expect("write fixture");
            path
        };
        Self {
            query: concatenate(
                "query.faa",
                &[
                    "docs/evidence/tlosan_stage_c/multi_query_20260924/query.faa",
                    "docs/evidence/tlosan_stage_c/multi_hsp_20260924/query.faa",
                ],
            ),
            subject: concatenate(
                "subject.fna",
                &[
                    "docs/evidence/tlosan_stage_c/lowercase_20260923/subjects.fna",
                    "docs/evidence/tlosan_stage_c/multi_query_20260924/subjects.fna",
                    "docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna",
                    "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna",
                ],
            ),
        }
    }

    fn names(&self) -> (String, String) {
        (
            self.query.display().to_string(),
            self.subject.display().to_string(),
        )
    }

    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        let (query_name, subject_name) = self.names();
        run(
            &self.query,
            &self.subject,
            (&query_name, &subject_name),
            formats,
            extra,
            observe,
        )
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.query);
        let _ = std::fs::remove_file(&self.subject);
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:289-301
// ```c
// CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
// results = lcl_blast.Run();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     formatter.PrintOneResultSet(**result, query);
// }
// ```
// This is the engine half of the test of the same name in `LOSAT/tests/run_local_tblastn.rs`,
// which compares every format with the CLI run of that format through `losat::api`. Here the
// reference is the same entry in process: a search that requests every format together, with the
// hit records and the observer, writes for each format the bytes of the search of that format
// alone without them, and the observer and the hit records change neither the reports nor the
// warnings.
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new();
    let mut hit_counts = Vec::new();
    for extra in [&[][..], &["-max_target_seqs", "1"][..]] {
        let together = inputs.run(&FORMATS, extra, true);
        assert!(!together.hits.is_empty(), "the fixture must produce hits");
        hit_counts.push(together.hits.len());
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.run(&[outfmt], extra, false);
            assert!(
                together.outputs[index] == alone.outputs[0],
                "run_local outfmt {outfmt} {extra:?}"
            );
            // The warnings go with the first format (outfmt 0 here); the formats that write no
            // query result warnings have none of the query warnings alone either.
            if !matches!(*outfmt, "5" | "8" | "9" | "11" | "12" | "15" | "16") {
                assert_eq!(
                    together.diagnostics, alone.diagnostics,
                    "warnings of {} formats equal those of the single-format run, outfmt {outfmt} {extra:?}",
                    FORMATS.len()
                );
            }
        }
    }
    assert!(
        hit_counts[1] < hit_counts[0],
        "-max_target_seqs 1 must drop subjects of this fixture: {hit_counts:?}"
    );
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3971-3980
// ```c++
// x_ShowAlnvecInfo(out,aln_vec_info,show_defline);
// ...
// out<<"\n";
// ```
#[test]
fn observer_ranges_are_exact_rows_and_sections_of_the_same_hsp() {
    let inputs = Inputs::new();
    for extra in [&[][..], &["-max_target_seqs", "1"][..]] {
        let result = inputs.run(&FORMATS, extra, true);
        assert_observer_ranges(&result, [0, 1, 2], &format!("tblastn {extra:?}"));
        // One row per HSP in the custom tabular report too (outfmt 10, no comment lines).
        assert_eq!(
            assert_tabular_ranges(&result, 10, None, "tblastn outfmt 10"),
            result.hits.len()
        );
        // A translated HSP section ends with the blank line of x_DisplayAlnvecInfo.
        for &(hsp, start, end) in &result.ranges[0] {
            assert!(
                result.outputs[0][start..end].ends_with(b"\n\n"),
                "section {hsp} {extra:?}"
            );
        }
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
    let queries = queries(&inputs.query);
    let subjects = subjects(&inputs.subject);
    let (query_name, subject_name) = inputs.names();
    // NCBI blast_args.cpp:2845-2851: a custom specification of outfmt 0 is ignored, so
    // "0 qseqid" is outfmt 0. A delimiter that a field separates its entries with is NCBI's
    // error (blast_args.cpp:2753-2792, `s_ValidateCustomDelim`).
    for outfmt in ["2", "6 delim=; qseqid staxids"] {
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
        let args = tblastn_args(&query_name, &subject_name, &[]);
        let result = run_local_tblastn(args, &queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}

// NCBI derives the hit list size per output format (blast_args.cpp:2894-2973): with
// `-num_descriptions 3` it is 250 for outfmt 0 and 500 for outfmt 6. One run searches once,
// so formats that need different hit list sizes are refused before anything is written
// (`blastinput/report_counts.rs`), while formats of one hit list size still share a run.
#[test]
fn formats_needing_different_hit_list_sizes_are_refused() {
    let inputs = Inputs::new();
    let queries = queries(&inputs.query);
    let subjects = subjects(&inputs.subject);
    let (query_name, subject_name) = inputs.names();
    let attempt = |formats: &[&str], extra: &[&str]| {
        let mut sinks: Vec<Vec<u8>> = formats.iter().map(|_| Vec::new()).collect();
        let mut diagnostics = Vec::new();
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
        let result = run_local_tblastn(
            tblastn_args(&query_name, &subject_name, extra),
            &queries,
            &subjects,
            &mut outputs,
        );
        drop(outputs);
        (result, sinks, diagnostics)
    };
    for (formats, extra) in [
        (&["0", "6"][..], &["-num_descriptions", "3"][..]),
        (&["6", "0"][..], &["-num_descriptions", "3"][..]),
        (&["0", "6"][..], &["-num_alignments", "3"][..]),
        (
            &["0", "6", "7"][..],
            &["-num_descriptions", "3", "-num_alignments", "2"][..],
        ),
    ] {
        let (result, sinks, _) = attempt(formats, extra);
        let error = format!("{:#}", result.expect_err("different hit list sizes"));
        assert!(
            error.contains("need different hit list sizes"),
            "{formats:?} {extra:?}: {error}"
        );
        // The option warnings of the first format may precede the refusal (BLASTP's
        // `check_options`, as the CLI writes a warning before an options error), so the
        // diagnostics are not checked; no report bytes are written.
        assert!(
            sinks.iter().all(Vec::is_empty),
            "{formats:?} {extra:?}: nothing is written"
        );
    }
    // The same formats without these options, or with one hit list size, run together.
    for extra in [&[][..], &["-max_target_seqs", "3"][..]] {
        let (result, sinks, _) = attempt(&["0", "6"], extra);
        result.unwrap_or_else(|error| panic!("{extra:?}: {error:#}"));
        assert!(
            !sinks[0].is_empty(),
            "{extra:?}: outfmt 0 writes its header"
        );
    }
}

// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312 (x_CleanAndCompress) and
// 4066 (`NStr::HtmlDecode`): NCBI reads past the end of the title `, ,` (and crashes when
// such a subject has hits) and decodes `&amp;` in outfmt 0; the tabular formats print
// only the ids. LOSAT writes the title `, ` stopped at the end of the string (approved
// exception 2 of PD-LOSAT-NCBI-DEFECTS) and the decoded title (`report/defline.rs`).
#[test]
fn outfmt0_titles_are_written_as_ncbi_makes_them() {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_tblastn_titles_{}_{nanos}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let query = dir.join("query.faa");
    std::fs::write(&query, ">pep\nWCDMTMSVIQVGGQFKPRTASAYFPYCISLCGKNQIVEHV\n").expect("query");
    let many = std::fs::read_to_string(repository(
        "LOSAT/tests/fasta/outfmt0/tblastx_many_subject.fasta",
    ))
    .expect("subject fixture");
    let sequence: String = many
        .lines()
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .collect();
    let queries = queries(&query);
    for (defline, heading) in [
        (", ,", "\n> , \n"),
        ("s &amp; t", "\n> s & t\n"),
        ("R&D; x", "\n> R&D; x\n"),
    ] {
        let subject = dir.join("subject.fna");
        std::fs::write(&subject, format!(">{defline}\n{sequence}\n")).expect("subject");
        let subjects = subjects(&subject);
        for outfmt in FORMATS {
            let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut report), &mut diagnostics);
            let result = run_local_tblastn(
                tblastn_args("q", "s", &[]),
                &queries,
                &subjects,
                &mut outputs,
            );
            drop(outputs);
            // BLAST XML: NCBI crashes on the punctuation-only title of a subject with hits
            // (outfmt 0, 5, 15 and 16 alike), so LOSAT refuses it and writes nothing.
            if matches!(outfmt, "5" | "15" | "16") && defline == ", ," {
                let error = result.expect_err("outfmt 5, 15 and 16 refuse the punctuation title");
                assert!(error.to_string().contains("crashes on"), "{error}");
                assert!(report.is_empty(), "nothing is written");
                continue;
            }
            result.unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            if outfmt == "0" {
                let report = String::from_utf8(report).expect("UTF-8 report");
                assert!(report.contains(heading), "{defline:?}: {report}");
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
