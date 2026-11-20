//! The TBLASTN search of `losat::api` (which calls the shared entry `run_local`) writes several
//! output formats from one search. Every format must be byte-identical to the CLI run that
//! requests only that format, and the warnings must equal that run's standard error. A subject or
//! a query without residues is kept, with NCBI's warnings, by the CLI. (The saved NCBI BLAST+
//! 2.17.0 bytes of three searches, the observer's ranges and the entry's refusals are checked in
//! the crate, `LOSAT/src/api/local_blast_tests/tblastn.rs`.)
//! (Identity with the previous release over the Stage G matrix is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::path::{Path, PathBuf};
use std::process::Command;

use losat::api::{read_fasta, FastaRecord, Outcome, Program, Role};
use run_local_support::{cli_out_path, temp_path, Files, TempFasta};

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

fn repository(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join(path)
}

/// The records of a TBLASTN query file or subject file as the library reads them.
fn records(path: &Path, role: Role) -> Vec<FastaRecord> {
    read_fasta(Program::Tblastn, role, &std::fs::read(path).expect("FASTA")).expect("records")
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

    /// The library's search of `formats` with the options `extra`, with the files' paths as the
    /// names that the reports show, as the CLI shows them.
    fn collect(&self, formats: &[&str], extra: &[&str]) -> Outcome {
        Files::new(Program::Tblastn, &self.query, &self.subject).collect(formats, extra)
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
    let out = cli_out_path(&inputs.query, outfmt);
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
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
// This half compares `losat::api`'s search with the CLI run of every format (reports and warnings);
// the half about the entry itself (the hit records and the observer change nothing) is the test of
// the same name in `LOSAT/src/api/local_blast_tests/tblastn.rs`.
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new();
    let mut hit_counts = Vec::new();
    for extra in [&[][..], &["-max_target_seqs", "1"][..]] {
        let together = inputs.collect(&FORMATS, extra);
        let rows = together.reports[1].spans.iter().flatten().count();
        assert_eq!(
            together.hsps.len(),
            rows,
            "one hit record per outfmt 6 row {extra:?}"
        );
        assert!(!together.hsps.is_empty(), "the fixture must produce hits");
        hit_counts.push(together.hsps.len());
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.collect(&[outfmt], extra);
            assert!(
                together.reports[index].bytes == alone.reports[0].bytes,
                "run_collect outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = cli(&inputs, outfmt, extra);
            assert!(
                together.reports[index].bytes == cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText
            //       ...
            //       || m_FormatType == CFormattingArgs::eJsonSeqalign
            //       || m_FormatType == CFormattingArgs::eSAM)
            //     {
            //         x_PrintStructuredReport(results, queries);
            //         return;
            //     }
            //     ...
            //     if (results.HasWarnings()) {
            //         ERR_POST(Warning << results.GetWarningStrings());
            //     }
            // ```
            // The warnings go with the first format (outfmt 0 here); outfmt 5, 8, 9, 11, 12, 15 and 16 write
            // no query result warnings, so their CLI runs equal `run_local` with that format alone.
            if matches!(*outfmt, "5" | "8" | "9" | "11" | "12" | "15" | "16") {
                // (The `-max_target_seqs 1` warning of the arguments stays.)
                assert!(
                    !String::from_utf8_lossy(&cli_stderr).contains("Query_3"),
                    "outfmt {outfmt}: no query warning: {}",
                    String::from_utf8_lossy(&cli_stderr)
                );
                assert_eq!(alone.diagnostics, cli_stderr, "outfmt {outfmt}");
                continue;
            }
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

/// One CLI run: the exit code, standard output and standard error.
fn cli_run(query: &Path, subject: &Path, extra: &[&str]) -> (Option<i32>, Vec<u8>, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
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
    let records = records(&whole, Role::Subject);
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
    let first = records(&query, Role::Query)[0].sequence.clone();
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980 and
// c++/src/algo/blast/core/blast_message.c:215-217 (quoted at `stage_d_pipeline::SetupEngineError`).
// Without the Gumbel block (`-searchsp` or `OLD_FSC`) a set of subjects whose average is below one letter
// fails the setup of a batch with NCBI's engine error (exit 3) after the start of the outfmt 0
// report, and after the report of a first batch that NCBI does not search (`BATCH_SIZE=1`, an
// invalid query first). NCBI 2.17.0 tblastn with these inputs (env parity probes f_probe/t of
// 2026-11-12, `tn.*`).
#[test]
fn a_short_average_subject_length_without_the_gumbel_block_is_ncbis_engine_error() {
    let query = repository("LOSAT/tests/fasta/psiblast/ps5c_qx_tk.faa");
    let subject = TempFasta::new(
        "tblastn_short_average.fna",
        &[("a", &b"A"[..]), ("b", &b""[..]), ("c", &b""[..])],
    );
    let path = subject.0.display().to_string();
    let subjects_warning = "Warning: [tblastn] Subject_2 b: Subject sequence contains no data\n\
                            Warning: [tblastn] Subject_3 c: Subject sequence contains no data\n";
    let error = "BLAST engine error: The average subject length is too short\n";
    // The prolog's database line wraps a long path as NCBI's does, so only its end is compared.
    let prolog_end = "           3 sequences; 1 total letters\n\n";
    for (args, env, first_batch) in [
        (&["-searchsp", "1000"][..], &[][..], false),
        (&[][..], &[("OLD_FSC", "1"), ("BATCH_SIZE", "1")][..], true),
    ] {
        for outfmt in ["0", "6", "7"] {
            let args: Vec<&str> = args.iter().copied().chain(["-outfmt", outfmt]).collect();
            let (code, stdout, stderr) =
                run_local_support::subject_cli("tblastn", &query, &subject.0, &args, env);
            assert_eq!(code, Some(3), "{args:?} {env:?}: {stderr}");
            let invalid = if first_batch {
                "Warning: [tblastn] Query_1 qx all X: Could not calculate ungapped Karlin-Altschul \
                 parameters due to an invalid query sequence or its translation. Please verify the \
                 query sequence(s) and/or filtering options \n"
            } else {
                ""
            };
            assert_eq!(
                stderr,
                format!("{subjects_warning}{invalid}{error}"),
                "{args:?} {env:?}"
            );
            match (outfmt, first_batch) {
                ("0", false) => assert!(
                    stdout.starts_with("TBLASTN 2.17.0+\n")
                        && stdout.contains("Database: User specified sequence set (Input:")
                        && stdout.ends_with(prolog_end),
                    "{stdout}"
                ),
                ("0", true) => assert!(
                    stdout.starts_with("TBLASTN 2.17.0+\n")
                        && stdout.contains("Database: User specified sequence set (Input:")
                        && stdout.ends_with(&format!(
                            "{prolog_end}\n\nQuery= qx all X\n\nLength=30\n\n\n***** No hits found *****\n\n\n\nLambda      K        H\n   -1.00    -1.00    -1.00 \n\nGapped\nLambda      K        H\n   -1.00    -1.00    -1.00 \n\nEffective search space used: 0\n"
                        )),
                    "{stdout}"
                ),
                ("7", true) => assert_eq!(
                    stdout,
                    format!(
                        "# TBLASTN 2.17.0+\n# Query: qx all X\n# Database: User specified sequence set (Input: {path})\n"
                    )
                ),
                _ => assert!(stdout.is_empty(), "{args:?} {env:?}: {stdout}"),
            }
        }
        // The XML2 and JSON documents: their start alone (`PrintProlog`, before the search), once;
        // a `-subject` run holds the first batch's report back, and the run never writes it
        // (review 1 of line DBID json, H1). NCBI 2.17.0 with these inputs.
        for (outfmt, start) in [
            ("15", "{\n\"BlastOutput2\": [\n".to_string()),
            (
                "16",
                "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastXML2\n    \
                 xmlns=\"http://www.ncbi.nlm.nih.gov\"\n    \
                 xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n    \
                 xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov \
                 http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"\n>"
                    .to_string(),
            ),
        ] {
            let args: Vec<&str> = args.iter().copied().chain(["-outfmt", outfmt]).collect();
            let (code, stdout, stderr) =
                run_local_support::subject_cli("tblastn", &query, &subject.0, &args, env);
            assert_eq!(code, Some(3), "{args:?} {env:?}: {stderr}");
            assert_eq!(
                stderr,
                format!("{subjects_warning}{error}"),
                "{args:?} {env:?}"
            );
            assert_eq!(stdout, start, "{args:?} {env:?}");
        }
    }
    // The boundary: an average of exactly one letter (tblastn `-subject` of "ACG", ``, ``)
    // is searched, exit 0, by NCBI 2.17.0 and LOSAT alike (review of env parity F, review_f/).
    let one = TempFasta::new(
        "tblastn_average_one.fna",
        &[("a", &b"ACG"[..]), ("b", &b""[..]), ("c", &b""[..])],
    );
    let (code, _, stderr) = run_local_support::subject_cli(
        "tblastn",
        &query,
        &one.0,
        &["-searchsp", "1000", "-outfmt", "6"],
        &[],
    );
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        format!(
            "{subjects_warning}Warning: [tblastn] Query_1 qx all X: Could not calculate ungapped \
             Karlin-Altschul parameters due to an invalid query sequence or its translation. \
             Please verify the query sequence(s) and/or filtering options \n"
        )
    );
}
