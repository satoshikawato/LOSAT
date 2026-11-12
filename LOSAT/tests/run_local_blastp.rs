//! The BLASTP search of `LOSAT::api` (which calls the shared entry `run_local`) writes several
//! output formats from one search. Every format must be byte-identical to the CLI run that
//! requests only that format, and the warnings must equal that run's standard error.
//! (Identity with the previous release over the regression manifests is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.) The entry's own contract (its observer, its
//! refusals, its errors) is checked in the crate, `LOSAT/src/api/local_blast_tests/blastp.rs`.

mod run_local_support;

use std::path::{Path, PathBuf};
use std::process::Command;

use run_local_support::{cli_out_path, read_records, temp_path, Files, TempFasta, Written};
use LOSAT::api::{Error, Outcome, Program};

const CUSTOM: &str =
    "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq btop";
const FORMATS: [&str; 11] = ["0", "6", "7", CUSTOM, "8", "9", "12", "5", "15", "16", "11"];

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
    let path = temp_path(&format!("{count}_{name}"));
    std::fs::write(&path, kept).expect("write fixture prefix");
    path
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

    /// The library's search of `formats` with the options `extra`, with the files' paths as the
    /// names that the reports show, as the CLI shows them.
    fn collect(&self, formats: &[&str], extra: &[&str]) -> Outcome {
        Files::new(Program::Blastp, &self.query, &self.subject).collect(formats, extra)
    }

    /// `Files::write` (`api::run`, which gives the entry no hit records and no observer, as the
    /// command line's call) of `formats`.
    fn write(&self, formats: &[&str], extra: &[&str]) -> Written {
        Files::new(Program::Blastp, &self.query, &self.subject).write(formats, extra)
    }

    /// The CLI output file and standard error of one single-format run.
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = cli_out_path(&self.query, outfmt);
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("blastp")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(&self.subject)
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

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.query);
        let _ = std::fs::remove_file(&self.subject);
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blastp_app.cpp:225-226
// ```c
// CBlastFormat formatter(opt, *db_adapter,
//                        fmt_args->GetFormattedOutputChoice(),
// ```
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
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new(8);
    for extra in [&[][..], &["-max_target_seqs", "1"][..]] {
        let together = inputs.collect(&FORMATS, extra);
        assert!(!together.hsps.is_empty(), "the fixture must produce hits");
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.collect(&[outfmt], extra);
            assert_eq!(
                together.reports[index].bytes, alone.reports[0].bytes,
                "run_collect outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert_eq!(
                together.reports[index].bytes, cli_output,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
// ```c++
//     if (m_FormatType == CFormattingArgs::eTabular ||
//         m_FormatType == CFormattingArgs::eTabularWithComments ||
//         m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
//         m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
// ```
// The tabular reports with a field list or a delimiter, and outfmt 10 and 20 (line CLI C2,
// `algorithm::blastp::tabular`), of one search: each equals its run alone and the CLI run of
// that format, beside the default 6.
#[test]
fn custom_tabular_formats_of_one_search_match_the_cli_run_of_that_format() {
    const TABULAR: [&str; 5] = [
        "6",
        "7 std qcovs qcovus qcovhsp sstrand staxid slen stitle",
        "10 qseqid sseqid qseq sseq btop nident positive ppos mismatch gaps",
        "20 delim=| qseqid sallseqid sacc salltitles frames qlen",
        "6 delim=, std",
    ];
    let inputs = Inputs::new(8);
    for extra in [&[][..], &["-num_threads", "2"][..]] {
        let together = inputs.collect(&TABULAR, extra);
        // The CLI's call has neither hit records nor an observer: `run` makes that call.
        let written = inputs.write(&TABULAR, extra);
        written.result.as_ref().expect("run");
        for (index, outfmt) in TABULAR.iter().enumerate() {
            assert!(
                !together.reports[index].bytes.is_empty(),
                "outfmt {outfmt:?}"
            );
            let alone = inputs.collect(&[outfmt], extra);
            assert!(
                together.reports[index].bytes == alone.reports[0].bytes,
                "run_collect outfmt {outfmt:?} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert!(
                written.outputs[index] == cli_output,
                "run outfmt {outfmt:?} {extra:?}"
            );
            assert_eq!(written.diagnostics, cli_stderr, "run {outfmt:?} {extra:?}");
            assert!(
                together.reports[index].bytes == cli_output,
                "CLI outfmt {outfmt:?} {extra:?}"
            );
            assert_eq!(together.diagnostics, cli_stderr, "{outfmt:?} {extra:?}");
        }
    }
}

/// The first `count` protein sequences of a fixture, from `start` residues on.
fn protein_sequences(name: &str, count: usize) -> Vec<Vec<u8>> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    read_records(&path)
        .into_iter()
        .take(count)
        .map(|record| record.seq().to_vec())
        .collect()
}

/// A BLASTP run of the CLI: the exit code, standard output and standard error.
fn cli_run(
    query: &Path,
    subject: &Path,
    extra: &[&str],
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(extra)
        .env_remove("LOSAT_TIMING")
        .env_remove("BATCH_SIZE");
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
// ```c
//         } catch (const CException& e) {
//             ...
//             CRef<CSearchMessage> m
//                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
//             messages[index].push_back(m);
//             s_InvalidateQueryContexts(qinfo, index);
//         }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
// ```c
//         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
//         		seqblk_vec->push_back(subj);
//         ...
//         		warning += "Subject sequence contains no data";
//         		ERR_POST(Warning << warning);
// ```
// A query without residues stays in its batch with NCBI's warning and a report without
// hits; a subject without residues gets its warning before the prolog and stays in the
// database statistics (oracle rows `blastp.rec_empty_mid.*`, `blastp.both_empty_records.*`).
#[test]
fn records_without_residues_are_kept_with_ncbis_warnings() {
    let q = protein_sequences("SicyWSV.faa", 3);
    let s = protein_sequences("PajaWSV.faa", 3);
    let query = TempFasta::new(
        "blastp_empty_records_query.faa",
        &[("q1", &q[0]), ("q2 empty two", b""), ("q3", &q[2])],
    );
    let subject = TempFasta::new(
        "blastp_empty_records_subject.faa",
        &[("s1", &s[0]), ("s2 no letters", b""), ("s3", &s[2])],
    );
    let run = Files::new(Program::Blastp, &query.0, &subject.0).collect(&["0", "6", "7"], &[]);
    let (outputs, diagnostics) = (&run.reports, &run.diagnostics);
    assert_eq!(
        String::from_utf8_lossy(diagnostics),
        "Warning: [blastp] Subject_2 s2 no letters: Subject sequence contains no data\n\
         Warning: [blastp] Query_2 q2 empty two: Sequence contains no data \n"
    );
    let pairwise = String::from_utf8_lossy(&outputs[0].bytes);
    assert!(pairwise.contains(" 3 sequences; "), "{pairwise}");
    assert!(
        pairwise.contains(
            "Query= q2 empty two\n\nLength=0\n\n\n***** No hits found *****\n\n\n\n\n\nEffective search space used: 0\n"
        ),
        "{pairwise}"
    );
    let commented = String::from_utf8_lossy(&outputs[2].bytes);
    assert!(
        commented.contains("# Query: q2 empty two\n# Database: ")
            && commented.contains("# 0 hits found\n")
            && commented.ends_with("# BLAST processed 3 queries\n"),
        "{commented}"
    );
    for (index, outfmt) in ["0", "6", "7"].into_iter().enumerate() {
        let (code, stdout, stderr) = cli_run(&query.0, &subject.0, &["-outfmt", outfmt], &[]);
        assert_eq!(code, Some(0), "{stderr}");
        assert!(stdout == outputs[index].bytes, "outfmt {outfmt}");
        assert_eq!(stderr.as_bytes(), &diagnostics[..], "outfmt {outfmt}");
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:729-732
// ```c
//    if (db_length == 0 &&
//        !BlastEffectiveLengthsOptions_IsSearchSpaceSet(eff_len_options)) {
//       return 0;
//    }
// ```
// Every subject without residues: the warnings, no hit, and a search space of 0; blastp
// is gapped, so the average subject length is not checked (oracle `blastp.rec_empty_all.s.*`).
// The library's search of the same files writes the CLI's report and warnings.
#[test]
fn subjects_without_residues_only_give_reports_without_hits() {
    let q = protein_sequences("SicyWSV.faa", 2);
    let query = TempFasta::new(
        "blastp_all_empty_query.faa",
        &[("q1", &q[0]), ("q2", &q[1])],
    );
    let subject = TempFasta::new(
        "blastp_all_empty_subject.faa",
        &[("s1 first", b""), ("s2", b"")],
    );
    let (code, stdout, stderr) = cli_run(&query.0, &subject.0, &["-outfmt", "0"], &[]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        "Warning: [blastp] Subject_1 s1 first: Subject sequence contains no data\n\
         Warning: [blastp] Subject_2 s2: Subject sequence contains no data\n"
    );
    let report = String::from_utf8_lossy(&stdout);
    assert!(report.contains(" 2 sequences; 0 total letters"), "{report}");
    assert_eq!(report.matches("***** No hits found *****").count(), 2);
    assert_eq!(
        report.matches("Effective search space used: 0\n").count(),
        2
    );
    let library = Files::new(Program::Blastp, &query.0, &subject.0).collect(&["0"], &[]);
    assert!(library.reports[0].bytes == stdout, "run_collect outfmt 0");
    assert_eq!(library.diagnostics, stderr.as_bytes());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
// ```c
//    if (sbp->gbp) {
//        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//        ...
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
//
//    if(min_subject_length <=0) {
// 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
//    }
// ```
// Without the Gumbel block (`-searchsp`, `OLD_FSC`) a subject set of fewer letters than sequences
// stops blastp and psiblast in the first batch's setup, after the subjects' warnings and the
// outfmt 0 prolog (NCBI 2.17.0 blastp and psiblast, `-subject` of `M`, ``, ``: exit 3, probes
// `x.subjshort.{sp,ofsc}.0` of line PS, item PS5c, and their blastp twins; the outfmt 6 and 7
// reports write nothing).
#[test]
fn subjects_of_less_than_one_letter_on_average_stop_the_search_without_the_gumbel_block() {
    let query = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/psiblast/tk.faa");
    let subject = TempFasta::new(
        "blastp_avg_short_subject.faa",
        &[
            ("s1 one letter", b"M"),
            ("s2 empty", b""),
            ("s3 empty", b""),
        ],
    );
    for (program, title, prolog_end) in [
        (
            "blastp",
            "BLASTP 2.17.0+",
            " 3 sequences; 1 total letters\n\n",
        ),
        (
            "psiblast",
            "PSIBLAST 2.17.0+",
            "Database: User specified sequence set.\n           3 sequences; 1 total letters\n\n",
        ),
    ] {
        for (args, env) in [
            (&["-searchsp", "1000", "-outfmt", "0"][..], None),
            (&["-outfmt", "0"][..], Some(("OLD_FSC", "1"))),
            (&["-searchsp", "1000", "-outfmt", "7"][..], None),
        ] {
            let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
            command
                .arg(program)
                .arg("-query")
                .arg(&query)
                .arg("-subject")
                .arg(&subject.0)
                .args(args)
                .env_remove("LOSAT_TIMING")
                .env_remove("BATCH_SIZE")
                .env_remove("OLD_FSC");
            if let Some((name, value)) = env {
                command.env(name, value);
            }
            let output = command.output().expect("run LOSAT CLI");
            let stdout = String::from_utf8_lossy(&output.stdout);
            assert_eq!(
                output.status.code(),
                Some(3),
                "{program} {args:?}: {stdout}"
            );
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                format!(
                    "Warning: [{program}] Subject_2 s2 empty: Subject sequence contains no data\n\
                     Warning: [{program}] Subject_3 s3 empty: Subject sequence contains no data\n\
                     BLAST engine error: The average subject length is too short\n"
                ),
                "{program} {args:?}"
            );
            if args.last() == Some(&"0") {
                assert!(
                    stdout.starts_with(title) && stdout.ends_with(prolog_end),
                    "{program} {args:?}: {stdout}"
                );
            } else {
                assert!(stdout.is_empty(), "{program} {args:?}: {stdout}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
// ```c++
//     int status = m_PrelimSearch->CheckInternalData();
//     if (status != 0)
//     {
//          // Search was not run, but we send back an empty CSearchResultSet.
// ```
// A first query batch of invalid queries only is not searched: NCBI reports it, then the next
// batch's setup fails (NCBI 2.17.0 `BATCH_SIZE=1 blastp` of an all-X query and tk with `-subject`
// of `MK`, ``, `` and `-searchsp 1000`: probes `q.bp.qx_tk.bs1.sp.{0,7}` of the PS5c review, exit
// 3; line PS DECISIONS 28).
#[test]
fn an_unsearched_first_batch_is_reported_before_the_next_batchs_engine_error() {
    let query = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/psiblast/ps5c_qx_tk.faa");
    let subject = TempFasta::new(
        "blastp_unsearched_first_subject.faa",
        &[("s1 two", b"MK"), ("s2 empty", b""), ("s3 empty", b"")],
    );
    for outfmt in ["0", "7"] {
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("blastp")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject.0)
            .args(["-searchsp", "1000", "-outfmt", outfmt])
            .env_remove("LOSAT_TIMING")
            .env_remove("OLD_FSC")
            .env("BATCH_SIZE", "1")
            .output()
            .expect("run LOSAT CLI");
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert_eq!(output.status.code(), Some(3), "outfmt {outfmt}: {stdout}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Warning: [blastp] Subject_2 s2 empty: Subject sequence contains no data\n\
             Warning: [blastp] Subject_3 s3 empty: Subject sequence contains no data\n\
             Warning: [blastp] Query_1 qx all X: Could not calculate ungapped Karlin-Altschul \
             parameters due to an invalid query sequence or its translation. Please verify the \
             query sequence(s) and/or filtering options \n\
             BLAST engine error: The average subject length is too short\n",
            "outfmt {outfmt}"
        );
        if outfmt == "0" {
            assert!(
                stdout.starts_with("BLASTP 2.17.0+")
                    && stdout.ends_with(
                        " 3 sequences; 2 total letters\n\n\n\nQuery= qx all X\n\nLength=30\n\n\n\
                         ***** No hits found *****\n\n\n\nLambda      K        H\n   -1.00    \
                         -1.00    -1.00 \n\nGapped\nLambda      K        H\n   -1.00    -1.00    \
                         -1.00 \n\nEffective search space used: 0\n"
                    ),
                "{stdout}"
            );
        } else {
            assert!(
                stdout.starts_with("# BLASTP 2.17.0+\n# Query: qx all X\n# Database: ")
                    && stdout.lines().count() == 3,
                "{stdout}"
            );
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:649-652
// ```c
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// A batch of queries without residues only stops NCBI after the prolog, with one message
// per query (oracle `blastp.rec_empty_all.q.*`); with BATCH_SIZE, after the reports of the
// batches before, without the epilog (`blastp.batch_late_allempty.*`). This half compares with
// the binary: the library's run (`api::run`) writes the prolog that the CLI writes to its
// standard output and fails with the message that the CLI writes to its standard error (the
// exit status 3 and the entry's own error are checked in the crate,
// `LOSAT/src/api/local_blast_tests/blastp.rs`); the late batch needs the CLI's BATCH_SIZE, which
// the library does not take, so only the CLI runs it.
#[test]
fn a_batch_of_queries_without_residues_stops() {
    let q = protein_sequences("SicyWSV.faa", 2);
    let s = protein_sequences("PajaWSV.faa", 2);
    let subject = TempFasta::new(
        "blastp_no_data_subject.faa",
        &[("s1", &s[0]), ("s2", &s[1])],
    );
    let query = TempFasta::new("blastp_no_data_query.faa", &[("q1", b""), ("q2", b"")]);
    let Written {
        result,
        outputs,
        diagnostics,
    } = Files::new(Program::Blastp, &query.0, &subject.0).write(&["0", "6"], &[]);
    let error = result.expect_err("the batch fails");
    assert!(matches!(error, Error::Search(_)), "{error:?}");
    assert_eq!(
        error.message(),
        "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n"
    );
    assert!(diagnostics.is_empty());
    let prolog = String::from_utf8_lossy(&outputs[0]);
    assert!(
        prolog.starts_with("BLASTP 2.17.0+\n") && prolog.ends_with(" total letters\n\n"),
        "{prolog}"
    );
    assert!(outputs[1].is_empty());
    let (code, stdout, stderr) = cli_run(&query.0, &subject.0, &["-outfmt", "0"], &[]);
    assert_eq!(code, Some(3), "{stderr}");
    assert_eq!(stdout, outputs[0]);
    assert_eq!(stderr, error.message());

    let mut text = format!(
        ">q1\n{}\n>q2\n{}\n",
        String::from_utf8_lossy(&q[0]),
        String::from_utf8_lossy(&q[1])
    );
    text.push_str(">q3\n>q4 empty\n");
    std::fs::write(&query.0, text).expect("write query");
    let batch = q[0].len().to_string();
    let (code, stdout, stderr) = cli_run(
        &query.0,
        &subject.0,
        &["-outfmt", "7"],
        &[("BATCH_SIZE", &batch)],
    );
    assert_eq!(code, Some(3), "{stderr}");
    let report = String::from_utf8_lossy(&stdout);
    assert!(
        report.contains("# Query: q1\n")
            && report.contains("# Query: q2\n")
            && !report.contains("# Query: q3")
            && !report.contains("# BLAST processed"),
        "{report}"
    );
    assert_eq!(
        stderr,
        "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n"
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
// A reader error in a later batch comes after the reports of the batches before (without
// the epilog); in the first batch, before any report. Blank and comment lines after the
// last record end the reading (oracle `blastp.batch_late_bad.*`). The batches are the CLI's
// (BATCH_SIZE), so the late error and the reports before it are checked on the binary alone; the
// library reads all queries before it searches, so for the first-batch error and for the blank
// and comment lines its run (`api::run_collect`) is held to the CLI's standard error and output.
#[test]
fn what_ends_the_query_reading_comes_with_its_batch() {
    let q = protein_sequences("SicyWSV.faa", 2);
    let s = protein_sequences("PajaWSV.faa", 1);
    let subject = TempFasta::new("blastp_batch_subject.faa", &[("s1", &s[0])]);
    let query = TempFasta::new("blastp_batch_query.faa", &[]);
    let text = format!(
        ">q1\n{}\n>q2\n{}\n",
        String::from_utf8_lossy(&q[0]),
        String::from_utf8_lossy(&q[1])
    );
    let batch = q[0].len().to_string();
    let mut bad = text.clone().into_bytes();
    bad.extend_from_slice(b">q3\n$%&()+= <>?@ [] {} ^~\nMKV\n");
    std::fs::write(&query.0, &bad).expect("write query");
    let (code, stdout, stderr) = cli_run(
        &query.0,
        &subject.0,
        &["-outfmt", "7"],
        &[("BATCH_SIZE", &batch)],
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
    let files = Files::new(Program::Blastp, &query.0, &subject.0);
    let error = files
        .try_collect(&["7"], &[])
        .expect_err("the reader fails");
    assert!(matches!(error, Error::Input(_)), "{error:?}");
    assert!(
        error
            .message()
            .starts_with("BLAST query error: CFastaReader: ")
            && stderr.starts_with(error.message()),
        "{error}"
    );
    let mut blank = text.into_bytes();
    blank.extend_from_slice(b"\n; comment\n\n");
    std::fs::write(&query.0, &blank).expect("write query");
    let (code, stdout, stderr) = cli_run(
        &query.0,
        &subject.0,
        &["-outfmt", "7"],
        &[("BATCH_SIZE", &batch)],
    );
    assert_eq!(code, Some(0), "{stderr}");
    assert!(String::from_utf8_lossy(&stdout).ends_with("# BLAST processed 2 queries\n"));
    let library = Files::new(Program::Blastp, &query.0, &subject.0).collect(&["7"], &[]);
    assert!(library.reports[0].bytes == stdout, "run_collect outfmt 7");
    assert_eq!(library.diagnostics, stderr.as_bytes());
}
