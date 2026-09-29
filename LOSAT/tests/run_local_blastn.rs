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
    // NCBI ignores a custom specification with outfmt 0 (blast_args.cpp:2845-2851).
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
        let result = run_local_blastn(inputs.args(&[]), &queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261-274,299-300
// ```c
//         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
//     ...
//                 if (!batch_size)
//                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
// ```
// The report of an invalid query depends on its query batch. The batches after the first
// depend on NCBI's extension counts, but they have at least 100 residues
// (blast_app_util.hpp:69-76, k_MinBatchSize), so invalid queries of fewer than 100
// residues that a valid query follows are searched with it. Otherwise outfmt 0 and 7,
// which print a section per query, fail; outfmt 6, which has none, does not.
#[test]
fn an_invalid_query_after_the_first_batch_fails_only_where_its_batch_matters() {
    let genome = fixture_sequence("LC738884.fasta");
    for (n_residues, sections_succeed) in [(40, true), (150, false)] {
        let all_n = vec![b'N'; n_residues];
        let query = TempFasta::new(
            "blastn_batches.fna",
            &[
                ("first_batch", &genome[0..6_000]),
                ("all_n", &all_n),
                ("valid", &genome[10_000..12_000]),
            ],
        );
        let subject = TempFasta::new("blastn_batches_subject.fna", &[("s", &genome[0..20_000])]);
        let inputs = Inputs { query, subject };
        let queries = read_records(&inputs.query.0);
        let subjects = read_records(&inputs.subject.0);
        for (outfmt, succeeds) in [
            ("0", sections_succeed),
            ("6", true),
            ("7", sections_succeed),
        ] {
            let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(&mut output),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            let result = run_local_blastn(inputs.args(&[]), &queries, &subjects, &mut outputs);
            drop(outputs);
            assert_eq!(
                result.is_ok(),
                succeeds,
                "{n_residues} N, outfmt {outfmt}: {result:?}"
            );
            if !succeeds {
                assert!(
                    output.is_empty(),
                    "outfmt {outfmt} writes nothing when it fails"
                );
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:672-674
// ```c
//     if (hsp->score == 0) {
//         return CRef<CSeq_align>();
//     }
// ```
// The preliminary search reads the N of the subject as random bases, so seeds of the query
// land in them; the traceback scores those HSPs 0, and no format or hit record shows them.
#[test]
fn hsps_of_score_zero_are_not_reported() {
    let inputs = Inputs {
        query: TempFasta::new("blastn_score0_query.fna", &[("q11", b"GTCTGTCACAA")]),
        subject: TempFasta::new("blastn_score0_subject.fna", &[("sN", &[b'N'; 100])]),
    };
    let run = inputs.run(
        &["6", "7"],
        &["-task", "blastn", "-word_size", "4", "-evalue", "1e6"],
        true,
    );
    assert!(run.outputs[0].is_empty(), "{:?}", run.outputs[0]);
    assert!(String::from_utf8_lossy(&run.outputs[1]).contains("# 0 hits found"));
    assert!(run.hits.is_empty());
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:977,989
// ```c
//     sv.SetCoding(CSeq_data::e_Ncbi4na);
//     ...
//     sv.GetStrandData(strand, buffer);
// ```
// Without -lcase_masking, the letter case of the query changes nothing: both strands come
// from the same case-free sequence data.
#[test]
fn a_lowercase_query_finds_the_same_hits_on_both_strands() {
    let genome = fixture_sequence("LC738884.fasta");
    let segment = &genome[20_000..22_000];
    let reverse: Vec<u8> = segment
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect();
    let subject = TempFasta::new(
        "blastn_case_subject.fna",
        &[("plus", segment), ("minus", &reverse)],
    );
    let lower: Vec<u8> = segment.to_ascii_lowercase();
    let mut outputs = Vec::new();
    for (name, sequence) in [("upper", segment), ("lower", &lower[..])] {
        let query = TempFasta::new(&format!("blastn_case_{name}.fna"), &[("q", sequence)]);
        // The subject file belongs to `subject`; the copy of its path in `inputs` is
        // forgotten below so that it does not delete the file.
        let inputs = Inputs {
            query,
            subject: TempFasta(subject.0.clone()),
        };
        for task in ["megablast", "blastn"] {
            let run = inputs.run(&["6"], &["-task", task], false);
            outputs.push((name, task, run.outputs[0].clone()));
        }
        std::mem::forget(inputs.subject);
    }
    for task in ["megablast", "blastn"] {
        let rows: Vec<&Vec<u8>> = outputs
            .iter()
            .filter(|(_, t, _)| *t == task)
            .map(|(_, _, rows)| rows)
            .collect();
        assert_eq!(rows[0], rows[1], "{task}: lowercase and uppercase queries");
        let text = String::from_utf8_lossy(rows[0]);
        assert!(
            text.lines().any(|row| row.contains("\tminus\t")),
            "{task}: the fixture has a minus-strand hit"
        );
    }
}
