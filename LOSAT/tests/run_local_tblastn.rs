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

use run_local_support::{assert_observer_ranges, reader_records, run_formats, Run, TempFasta};
use LOSAT::algorithm::tblastn::TblastnArgs;
use LOSAT::api::local_blast::{run_local_tblastn, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
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
    let queries = queries(&inputs.query);
    let subjects = subjects(&inputs.subject);
    let (query_name, subject_name) = inputs.names();
    // NCBI blast_args.cpp:2845-2851: a custom specification of outfmt 0 is ignored, so
    // "0 qseqid" is outfmt 0; LOSAT's TBLASTN writes no custom tabular fields or delimiter.
    for outfmt in ["5", "6 qseq sseq", "6 delim=,"] {
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
            run_local_tblastn(
                tblastn_args("q", "s", &[]),
                &queries,
                &subjects,
                &mut outputs,
            )
            .unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            drop(outputs);
            if outfmt == "0" {
                let report = String::from_utf8(report).expect("UTF-8 report");
                assert!(report.contains(heading), "{defline:?}: {report}");
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// One CLI run: the exit code, standard output and standard error.
fn cli_run(query: &Path, subject: &Path, extra: &[&str]) -> (Option<i32>, Vec<u8>, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(extra)
        .output()
        .expect("run LOSAT CLI");
    (
        output.status.code(),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
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
// NCBI reference: c++/src/algo/blast/core/blast_setup.c:969-979
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
//	   return BLASTERR_SUBJECT_LENGTH_INVALID;
// ```
// A subject without residues gets its warning before the search (after `Query is Empty!`),
// counts in the database statistics and is never searched; the gapped tblastn raises the
// shortest subject to 10 letters, so subjects without letters only never stop the search.
#[test]
fn subjects_without_residues_are_kept_with_ncbis_warnings() {
    let query = repository("docs/evidence/tlosan_stage_c/run_20260923/query.faa");
    let whole = repository("docs/evidence/tlosan_stage_c/run_20260923/subjects.fna");
    let records = subjects(&whole);
    assert!(records.len() >= 2, "the fixture has two subjects or more");
    let mut written: Vec<(String, Vec<u8>)> = Vec::new();
    for (index, record) in records.iter().enumerate() {
        written.push((
            String::from_utf8_lossy(&record.title).into_owned(),
            record.sequence.clone(),
        ));
        if index == 0 {
            written.push(("s_empty no letters".to_string(), Vec::new()));
        }
    }
    let with_empty = TempFasta::new(
        "tblastn_empty_subject.fna",
        &written
            .iter()
            .map(|(title, sequence)| (title.as_str(), sequence.as_slice()))
            .collect::<Vec<_>>(),
    );
    let without_empty = TempFasta::new(
        "tblastn_no_empty_subject.fna",
        &written
            .iter()
            .filter(|(_, sequence)| !sequence.is_empty())
            .map(|(title, sequence)| (title.as_str(), sequence.as_slice()))
            .collect::<Vec<_>>(),
    );
    let (code, rows, stderr) = cli_run(&query, &with_empty.0, &["-outfmt", "6"]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        "Warning: [tblastn] Subject_2 s_empty no letters: Subject sequence contains no data\n"
    );
    let (_, expected_rows, _) = cli_run(&query, &without_empty.0, &["-outfmt", "6"]);
    assert!(!rows.is_empty() && rows == expected_rows);
    let (_, report, _) = cli_run(&query, &with_empty.0, &["-outfmt", "7"]);
    assert!(String::from_utf8_lossy(&report).contains("hits found"));
    let (_, report, _) = cli_run(&query, &with_empty.0, &[]);
    let report = String::from_utf8_lossy(&report).into_owned();
    assert!(
        report.contains(&format!(" {} sequences; ", written.len()))
            && report.contains(&format!(
                "Number of sequences in database:  {}\n",
                written.len()
            )),
        "{report}"
    );
    // Every subject without letters: NCBI warns about each and reports no hits (exit 0).
    let empty = TempFasta::new("tblastn_all_empty.fna", &[("e1", b""), ("", b"")]);
    for outfmt in ["0", "6", "7"] {
        let (code, report, stderr) = cli_run(&query, &empty.0, &["-outfmt", outfmt]);
        assert_eq!(code, Some(0), "{stderr}");
        assert_eq!(
            stderr,
            "Warning: [tblastn] Subject_1 e1: Subject sequence contains no data\nWarning: [tblastn] Subject_2 : Subject sequence contains no data\n"
        );
        let report = String::from_utf8_lossy(&report).into_owned();
        match outfmt {
            "0" => assert!(report.contains("***** No hits found *****"), "{report}"),
            "6" => assert!(report.is_empty(), "{report}"),
            _ => assert!(report.contains("# 0 hits found"), "{report}"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
// ```c++
//         } catch (const CException& e) {
//             ...
//             CRef<CSearchMessage> m
//                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
//             messages[index].push_back(m);
//             s_InvalidateQueryContexts(qinfo, index);
//         }
//     ...
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// A query without residues is set up without data: its warning comes before its report and
// the other queries' rows do not change; a batch of such queries only stops the run with
// NCBI's set-up error (exit 3), after the outfmt 0 prolog (oracle `tblastn.batch_late_empty.*`,
// `tblastn.both_all_empty_records.*`).
#[test]
fn queries_without_residues_are_kept_and_a_batch_of_them_stops() {
    let query = repository("docs/evidence/tlosan_stage_c/run_20260923/query.faa");
    let subject = repository("docs/evidence/tlosan_stage_c/run_20260923/subjects.fna");
    let first = queries(&query)[0].sequence.clone();
    let with_empty = TempFasta::new(
        "tblastn_empty_query.faa",
        &[("q1 first", first.as_slice()), ("q2 no letters", &b""[..])],
    );
    let (code, rows, stderr) = cli_run(&with_empty.0, &subject, &["-outfmt", "6"]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        "Warning: [tblastn] Query_2 q2 no letters: Sequence contains no data \n"
    );
    let only_first = TempFasta::new("tblastn_first_query.faa", &[("q1 first", first.as_slice())]);
    let (_, expected_rows, _) = cli_run(&only_first.0, &subject, &["-outfmt", "6"]);
    assert!(!rows.is_empty() && rows == expected_rows);
    let all_empty = TempFasta::new(
        "tblastn_all_empty_query.faa",
        &[("e1", &b""[..]), ("e2", &b""[..])],
    );
    for outfmt in ["0", "6", "7"] {
        let (code, stdout, stderr) = cli_run(&all_empty.0, &subject, &["-outfmt", outfmt]);
        assert_eq!(code, Some(3), "{stderr}");
        assert_eq!(
            stderr,
            "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n"
        );
        let report = String::from_utf8_lossy(&stdout);
        if outfmt == "0" {
            assert!(
                report.starts_with("TBLASTN 2.17.0+\n") && report.ends_with(" total letters\n\n"),
                "{report}"
            );
        } else {
            assert!(report.is_empty(), "{report}");
        }
    }
}
