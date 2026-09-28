#![allow(warnings, clippy::all)]

//! The shared TBLASTX entry `run_local` writes every requested output from one search.
//! TBLASTX implements outfmt 6 only, so two outfmt 6 outputs are requested; each must
//! be byte-identical to the CLI run, the warnings must equal that run's standard error,
//! and the observer must report the exact row of every HSP. With `-num_threads 2` the
//! native search reduces its results in a collector thread; the result must not change.
//! (Identity with the previous release over the regression manifests is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::process::Command;

use run_local_support::{
    assert_tabular_ranges, fixture_sequence, read_records, run_formats, Run, TempFasta,
};
use LOSAT::algorithm::tblastx::TblastxArgs;
use LOSAT::api::local_blast::{run_local_tblastx, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 2] = ["6", "6"];

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

    /// The CLI output file and standard error of one run.
    fn cli(&self, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = self.query.0.with_extension("cli.out");
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
            .arg("tblastx")
            .arg("-query")
            .arg(&self.query.0)
            .arg("-subject")
            .arg(&self.subject.0)
            .args(extra)
            .arg("-outfmt")
            .arg("6")
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
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
#[test]
fn every_output_of_one_search_matches_the_cli_run() {
    let inputs = Inputs::new();
    let serial = inputs.run(&FORMATS, &[], true);
    assert!(
        serial.hits.is_empty(),
        "TBLASTX produces no hit records yet"
    );
    assert!(
        !serial.ranges[0].is_empty(),
        "the fixture must produce hits"
    );
    for extra in [
        &[][..],
        &["-num_threads", "2"][..],
        &["-max_target_seqs", "1"][..],
    ] {
        let together = inputs.run(&FORMATS, extra, true);
        let (cli_output, cli_stderr) = inputs.cli(extra);
        for (index, outfmt) in FORMATS.iter().enumerate() {
            assert!(
                together.outputs[index] == cli_output,
                "CLI outfmt {outfmt} (output {index}) {extra:?}"
            );
        }
        assert_eq!(
            together.diagnostics, cli_stderr,
            "warnings of two outputs equal those of the CLI run {extra:?}"
        );
        if extra.first() == Some(&"-num_threads") {
            assert!(
                together.outputs[0] == serial.outputs[0],
                "the collector thread must not change the result"
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
#[test]
fn observer_ranges_are_exact_rows_of_the_same_hsp() {
    let inputs = Inputs::new();
    for extra in [&[][..], &["-num_threads", "2"][..]] {
        let result = inputs.run(&FORMATS, extra, true);
        let context = format!("tblastx {extra:?}");
        let rows = assert_tabular_ranges(&result, 0, None, &context);
        assert!(rows > 0, "the fixture must produce hits {extra:?}");
        assert_tabular_ranges(&result, 1, None, &context);
        assert_eq!(
            result.ranges[0], result.ranges[1],
            "{context}: both outputs"
        );
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
    for outfmt in ["0", "7", "6 qseqid"] {
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
