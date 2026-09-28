#![allow(warnings, clippy::all)]

//! The shared BLASTP entry `run_local` writes several output formats from one search.
//! Every format must be byte-identical to the CLI run that requests only that format,
//! the warnings must equal that run's standard error, and the observer must report the
//! exact bytes of every HSP in every format.
//! (Identity with the previous release over the regression manifests is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use run_local_support::{assert_observer_ranges, field, read_records, run_formats, Run};
use LOSAT::api::local_blast::{run_local_blastp, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::cli::{try_parse_from, Cli, Commands};

const CUSTOM: &str =
    "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq btop";
const FORMATS: [&str; 4] = ["0", "6", "7", CUSTOM];

/// A temporary FASTA file holding the first `count` records of a fixture, byte for byte.
fn fixture_prefix(name: &str, count: usize) -> PathBuf {
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    let text = std::fs::read_to_string(source).expect("fixture FASTA");
    let mut records = 0;
    let mut kept = String::new();
    for line in text.split_inclusive('\n') {
        if line.starts_with('>') {
            records += 1;
            if records > count {
                break;
            }
        }
        kept.push_str(line);
    }
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "losat_run_local_{}_{nanos}_{count}_{name}",
        std::process::id()
    ));
    std::fs::write(&path, kept).expect("write fixture prefix");
    path
}

fn blastp_args(
    query: &Path,
    subject: &Path,
    extra: &[&str],
) -> LOSAT::algorithm::blastp::BlastpArgs {
    let mut argv = vec![
        "LOSAT".to_string(),
        "blastp".to_string(),
        "-query".to_string(),
        query.display().to_string(),
        "-subject".to_string(),
        subject.display().to_string(),
    ];
    argv.extend(extra.iter().map(|arg| arg.to_string()));
    let cli: Cli = try_parse_from(argv).expect("blastp argv");
    match cli.command {
        Commands::Blastp(args) => args,
        _ => unreachable!("blastp subcommand"),
    }
}

struct Inputs {
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new(records: usize) -> Self {
        Self {
            query: fixture_prefix("SicyWSV.faa", records),
            subject: fixture_prefix("PajaWSV.faa", records),
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.query);
        let _ = std::fs::remove_file(&self.subject);
    }
}

fn run(inputs: &Inputs, formats: &[&str], extra: &[&str], observe: bool) -> Run {
    let queries = read_records(&inputs.query);
    let subjects = read_records(&inputs.subject);
    run_formats(formats, observe, |outputs| {
        let args = blastp_args(&inputs.query, &inputs.subject, extra);
        run_local_blastp(args, &queries, &subjects, "", "", outputs)
    })
}

/// The CLI output file and standard error of one single-format run.
fn cli(inputs: &Inputs, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
    let out = inputs
        .query
        .with_extension(format!("cli{}.out", outfmt.len()));
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("blastp")
        .arg("-query")
        .arg(&inputs.query)
        .arg("-subject")
        .arg(&inputs.subject)
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

// NCBI reference: ncbi-blast/c++/src/app/blast/blastp_app.cpp:225-226
// ```c
// CBlastFormat formatter(opt, *db_adapter,
//                        fmt_args->GetFormattedOutputChoice(),
// ```
// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-465
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
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new(8);
    for extra in [&[][..], &["-max_target_seqs", "1"][..]] {
        let together = run(&inputs, &FORMATS, extra, true);
        assert!(!together.hits.is_empty(), "the fixture must produce hits");
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = run(&inputs, &[outfmt], extra, false);
            assert_eq!(
                together.outputs[index], alone.outputs[0],
                "run_local outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = cli(&inputs, outfmt, extra);
            assert_eq!(
                together.outputs[index], cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            assert_eq!(
                together.diagnostics,
                cli_stderr,
                "warnings of {} formats equal those of the single-format CLI run {extra:?}",
                FORMATS.len()
            );
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
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
    let inputs = Inputs::new(8);
    for extra in [&[][..], &["-max_target_seqs", "1"][..]] {
        let result = run(&inputs, &FORMATS, extra, true);
        assert_observer_ranges(&result, [0, 1, 2], &format!("blastp {extra:?}"));

        // The custom list: rows in hit-list order that tile the output, with the same
        // coordinates as the default row of the same index.
        let (custom, rows) = (&result.outputs[3], &result.ranges[3]);
        assert_eq!(rows.len(), result.hits.len(), "custom rows {extra:?}");
        let mut position = 0;
        for (expected_index, &(hsp, start, end)) in rows.iter().enumerate() {
            assert_eq!(hsp, expected_index);
            assert_eq!(start, position);
            let row = &custom[start..end];
            assert_eq!(row.iter().filter(|&&byte| byte == b'\n').count(), 1);
            let (_, s6, e6) = result.ranges[1][hsp];
            let row6 = &result.outputs[1][s6..e6];
            assert_eq!(
                [4, 5, 6, 7].map(|index| field(row, index)),
                [6, 7, 8, 9].map(|index| field(row6, index)),
                "custom HSP {hsp} {extra:?}"
            );
            position = end;
        }
        assert_eq!(position, custom.len());
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
// ```c
// if (args[kArgOutputFormat]) {
//     string fmt_choice =
//         NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
// ```
#[test]
fn unsupported_formats_fail_before_searching() {
    let inputs = Inputs::new(1);
    let queries = read_records(&inputs.query);
    let subjects = read_records(&inputs.subject);
    for outfmt in ["0 qseqid", "9", "6 nosuchfield"] {
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
        let args = blastp_args(&inputs.query, &inputs.subject, &[]);
        let result = run_local_blastp(args, &queries, &subjects, "", "", &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}
