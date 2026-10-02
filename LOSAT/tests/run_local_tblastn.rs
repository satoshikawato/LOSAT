#![allow(warnings, clippy::all)]

//! The shared TBLASTN entry `run_local` writes several output formats from one search.
//! Every format must be byte-identical to the saved NCBI BLAST+ 2.17.0 bytes and to the
//! CLI run that requests only that format, the warnings must be written once, and the
//! observer must report the exact bytes of every HSP in every format.
//! (Identity with the previous release over the Stage G matrix is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use run_local_support::{assert_observer_ranges, read_records, run_formats, Run};
use LOSAT::algorithm::tblastn::TblastnArgs;
use LOSAT::api::local_blast::{run_local_tblastn, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 3] = ["0", "6", "7"];

fn repository(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join(path)
}

/// Parses the TBLASTN options. `query` and `subject` are display names only.
fn tblastn_args(query: &str, subject: &str, extra: &[&str]) -> TblastnArgs {
    let mut argv = vec!["LOSAT", "tblastn", "-query", query, "-subject", subject];
    argv.extend_from_slice(extra);
    let cli: Cli = try_parse_from(argv).expect("tblastn argv");
    match cli.command {
        Commands::Tblastn(args) => args,
        _ => unreachable!("tblastn subcommand"),
    }
}

fn run(
    query: &Path,
    subject: &Path,
    names: (&str, &str),
    formats: &[&str],
    extra: &[&str],
    observe: bool,
) -> Run {
    let queries = read_records(query);
    let subjects = read_records(subject);
    run_formats(formats, observe, |outputs| {
        run_local_tblastn(
            tblastn_args(names.0, names.1, extra),
            &queries,
            &subjects,
            outputs,
        )
    })
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
            &FORMATS,
            &[],
            true,
        );
        for (index, outfmt) in FORMATS.iter().enumerate() {
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
            let nanos = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("clock")
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "losat_run_local_tblastn_{}_{nanos}_{name}",
                std::process::id()
            ));
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

/// The CLI output file and standard error of one single-format run.
fn cli(inputs: &Inputs, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
    let out = inputs.query.with_extension(format!("cli{outfmt}.out"));
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("tblastn")
        .arg("-query")
        .arg(&inputs.query)
        .arg("-subject")
        .arg(&inputs.subject)
        .args(extra)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-out")
        .arg(&out)
        .output()
        .expect("run LOSAT CLI");
    assert!(output.status.success(), "CLI failed: {output:?}");
    let bytes = std::fs::read(&out).expect("CLI output");
    let _ = std::fs::remove_file(&out);
    (bytes, output.stderr)
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
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
// ```c
// if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
//    hitlist_size = args[kArgMaxTargetSequences].AsInteger();
// }
// ```
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
            let (cli_output, cli_stderr) = cli(&inputs, outfmt, extra);
            assert!(
                together.outputs[index] == cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            assert!(!cli_stderr.is_empty(), "the invalid query is warned about");
            assert_eq!(
                together.diagnostics, cli_stderr,
                "warnings of three formats equal those of the single-format CLI run"
            );
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
    let queries = read_records(&inputs.query);
    let subjects = read_records(&inputs.subject);
    let (query_name, subject_name) = inputs.names();
    for outfmt in ["5", "6 qseq sseq", "0 qseqid"] {
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

// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312 (x_CleanAndCompress) and
// 4066 (`NStr::HtmlDecode`): NCBI reads past the end of the title `, ,` (and crashes when
// such a subject has hits) and decodes `&amp;` in outfmt 0; the tabular formats print
// only the ids.
#[test]
fn outfmt0_rejects_subject_titles_that_ncbi_reads_past_or_decodes() {
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
    for (defline, reason) in [
        (", ,", "reads past its end"),
        ("s &amp; t", "HTML character reference"),
    ] {
        let subject = dir.join("subject.fna");
        std::fs::write(&subject, format!(">{defline}\n{sequence}\n")).expect("subject");
        let queries = read_records(&query);
        let subjects = read_records(&subject);
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
            if outfmt == "0" {
                let error = result
                    .expect_err("outfmt 0 must reject the title")
                    .to_string();
                assert!(
                    error.contains(reason) && error.contains("not supported by LOSAT's TBLASTN"),
                    "{error}"
                );
                assert!(report.is_empty() && diagnostics.is_empty());
            } else {
                result.unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
