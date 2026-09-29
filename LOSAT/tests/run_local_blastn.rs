#![allow(warnings, clippy::all)]

//! The shared BLASTN entry `run_local` writes several output formats from one search.
//! Every format must be byte-identical to the CLI run that requests only that format,
//! the warnings must equal that run's standard error, and the observer must report the
//! exact row and pairwise section of every HSP. (Identity of outfmt 0 with NCBI is
//! checked by docs/evidence/losat_web_e2a/check_losat.py over the frozen fixtures.) (Identity with the previous release over the regression
//! manifests is checked by docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::process::Command;

use run_local_support::{
    assert_observer_ranges, fixture_sequence, read_records, run_formats, Run, TempFasta,
};
use LOSAT::algorithm::blastn::BlastnArgs;
use LOSAT::api::local_blast::{run_local_blastn, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 3] = ["0", "6", "7"];

/// Three queries against four subjects cut from two genomes. The first query overlaps
/// two subjects, so `-max_target_seqs 1` removes a subject from its hit list.
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
                "blastn_query.fna",
                &[
                    ("q_a0", &a[0..20_000]),
                    ("q_a100k", &a[100_000..120_000]),
                    ("q_b50k", &b[50_000..60_000]),
                ],
            ),
            subject: TempFasta::new(
                "blastn_subject.fna",
                &[
                    ("s_a10k", &a[10_000..40_000]),
                    ("s_a90k", &a[90_000..130_000]),
                    ("s_a5k", &a[5_000..25_000]),
                    ("s_b45k", &b[45_000..70_000]),
                ],
            ),
        }
    }

    fn args(&self, extra: &[&str]) -> BlastnArgs {
        let query = self.query.0.display().to_string();
        let subject = self.subject.0.display().to_string();
        let mut argv = vec![
            "LOSAT", "blastn", "-query", &query, "-subject", &subject, "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv).expect("blastn argv");
        match cli.command {
            Commands::Blastn(args) => args,
            _ => unreachable!("blastn subcommand"),
        }
    }

    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        let queries = read_records(&self.query.0);
        let subjects = read_records(&self.subject.0);
        run_formats(formats, observe, |outputs| {
            run_local_blastn(self.args(extra), &queries, &subjects, outputs)
        })
    }

    /// The CLI output file and standard error of one single-format run.
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = self.query.0.with_extension(format!("cli{outfmt}.out"));
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
            .arg("blastn")
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
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
// ```c
// if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
//    hitlist_size = args[kArgMaxTargetSequences].AsInteger();
// }
// ```
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new();
    let mut row_counts = Vec::new();
    for extra in [
        &[][..],
        &["-max_target_seqs", "1"][..],
        &["-task", "blastn"][..],
        &["-num_threads", "2"][..],
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
                "warnings of two formats equal those of the single-format CLI run {extra:?}"
            );
        }
    }
    assert!(row_counts[0] > 0, "the fixture must produce hits");
    assert!(
        row_counts[1] < row_counts[0],
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
        &["-max_target_seqs", "1"][..],
        &["-task", "blastn"][..],
    ] {
        let result = inputs.run(&FORMATS, extra, true);
        assert_observer_ranges(&result, [0, 1, 2], &format!("blastn {extra:?}"));
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
    let inputs = Inputs::new();
    let queries = read_records(&inputs.query.0);
    let subjects = read_records(&inputs.subject.0);
    for outfmt in ["0 qseqid", "5", "6 qseqid"] {
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
        let result = run_local_blastn(inputs.args(&[]), &queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}
