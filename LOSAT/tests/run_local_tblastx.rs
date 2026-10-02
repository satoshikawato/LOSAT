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
}

// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312 (x_CleanAndCompress) and
// 4066 (`NStr::HtmlDecode`): NCBI reads past the end of the title `, ,` (and crashes when
// such a subject has hits) and decodes `&amp;` in outfmt 0; the tabular formats print
// only the ids.
#[test]
fn outfmt0_rejects_subject_titles_that_ncbi_reads_past_or_decodes() {
    let inputs = Inputs::new();
    let queries = read_records(&inputs.query.0);
    let sequence = fixture_sequence("LC738884.fasta");
    for (defline, reason) in [
        (", ,", "reads past its end"),
        ("s &amp; t", "HTML character reference"),
    ] {
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
            if outfmt == "0" {
                let error = result
                    .expect_err("outfmt 0 must reject the title")
                    .to_string();
                assert!(
                    error.contains(reason) && error.contains("not supported by LOSAT's TBLASTX"),
                    "{error}"
                );
                assert!(report.is_empty() && diagnostics.is_empty());
            } else {
                result.unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            }
        }
    }
}

// NCBI reference: c++/src/objtools/readers/fasta.cpp:919-935 (CFastaReader: a residue that is
// not an IUPAC nucleotide letter is removed with a warning) and the "Sequence contains no
// data" warning of a record without residues; LOSAT rejects both, as BLASTN does. NCBI culls
// with hspfilter_culling.c, which LOSAT's culling does not reproduce.
#[test]
fn inputs_ncbi_reads_differently_and_culling_are_rejected() {
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
    let error = run(&["-culling_limit", "2"], &queries, &subjects);
    assert!(
        error.contains("-culling_limit 2 is not supported by LOSAT's TBLASTX"),
        "{error}"
    );
}
