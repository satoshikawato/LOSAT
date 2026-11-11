//! The shared BLASTX entry `run_local` in process: what the entry itself does (its observer, its
//! refusals of hand-set arguments and of format sets, its typed errors), not what the command line
//! writes. The tests that compare the entry's reports with the `losat` binary go through
//! `LOSAT::api` (`LOSAT/tests/run_local_blastx.rs`); the tests of the same name here check the
//! halves of those that the api does not reach.

use std::path::{Path, PathBuf};

use super::support::{
    assert_observer_ranges, assert_tabular_ranges, fixture_sequence, reader_records,
    run_formats_with, temp_path, Run,
};
use crate::algorithm::blastx::BlastxArgs;
use crate::api::local_blast::{
    run_local_blastx, with_run_pool, FormatOutput, OutputSink, ReportOutputs,
};
use crate::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use crate::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 3] = ["0", "6", "7"];

/// The LOSAT crate's directory, where the manifests' paths start.
fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Parses the BLASTX options. `query` and `subject` are display names only.
fn blastx_args(query: &str, subject: &str, extra: &[&str]) -> BlastxArgs {
    let mut argv = vec!["losat", "blastx", "-query", query, "-subject", subject];
    argv.extend_from_slice(extra);
    let cli: Cli = try_parse_from(argv).expect("blastx argv");
    match cli.command {
        Commands::Blastx(args) => args,
        _ => unreachable!("blastx subcommand"),
    }
}

/// The queries as NCBI's nucleotide reader reads them (LOSAT Web's `register`).
fn queries(path: &Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::query("BLASTX", false, true))
}

/// The subjects as NCBI's protein reader reads them (LOSAT Web's `register`).
fn subjects(path: &Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::subject("BLASTX", true, true))
}

/// `run_local` over `queries` and `subjects` (shown as `inputs`' names), with a sink of its own for
/// each format: the result, the outputs and the diagnostics.
fn try_run(
    inputs: &Inputs,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    formats: &[&str],
    extra: &[&str],
) -> (anyhow::Result<()>, Vec<Vec<u8>>, Vec<u8>) {
    let mut sinks = vec![Vec::new(); formats.len()];
    let mut diagnostics = Vec::new();
    let result = {
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
        run_local_blastx(
            blastx_args(&inputs.query, &inputs.subject, extra),
            queries,
            subjects,
            &mut outputs,
        )
    };
    (result, sinks, diagnostics)
}

/// The input paths of a run, as given to the CLI (relative to the crate's directory or
/// absolute); `run_local` shows them as its display names.
struct Inputs {
    query: String,
    subject: String,
}

impl Inputs {
    fn new(query: impl AsRef<Path>, subject: impl AsRef<Path>) -> Self {
        Self {
            query: query.as_ref().display().to_string(),
            subject: subject.as_ref().display().to_string(),
        }
    }

    fn path(name: &str) -> PathBuf {
        crate_dir().join(name)
    }

    /// `run_local` with every format of `formats`.
    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        self.run_with(formats, extra, observe, observe)
    }

    /// `run`, with the observer and the hit records requested apart.
    fn run_with(&self, formats: &[&str], extra: &[&str], observe: bool, hits: bool) -> Run {
        let queries = queries(&Self::path(&self.query));
        let subjects = subjects(&Self::path(&self.subject));
        run_formats_with(formats, observe, hits, |outputs| {
            run_local_blastx(
                blastx_args(&self.query, &self.subject, extra),
                &queries,
                &subjects,
                outputs,
            )
        })
    }
}

/// A query of several reader batches (more than 10,002 letters in all): a record without letters
/// and a long record, then records with the reader's warnings, a record with a lower-case island
/// and a query too short to search; the subjects hold a record without letters.
struct Batches {
    query: PathBuf,
    subject: PathBuf,
}

impl Batches {
    fn new() -> Self {
        let genome = fixture_sequence("MjeNMV.fasta");
        // The first batch: a record without letters and one of 10,100 letters.
        let mut query = b">empty record\n".to_vec();
        query.extend_from_slice(b">long record of 10100 letters\n");
        for line in genome[0..10_100].chunks(70) {
            query.extend_from_slice(line);
            query.push(b'\n');
        }
        // The reader's warnings (invalid residues, hyphens) in the second batch.
        query.extend_from_slice(
            &std::fs::read(Inputs::path(
                "tests/fixtures/fasta_input/batch_late_warn.nuc.fa",
            ))
            .expect("fixture"),
        );
        query.extend_from_slice(b">lower case island\n");
        query.extend_from_slice(&genome[20_000..20_300]);
        query.extend(genome[20_300..20_600].iter().map(u8::to_ascii_lowercase));
        query.extend_from_slice(&genome[20_600..21_200]);
        query.extend_from_slice(b"\n>short\nAC\n");
        let query_path = temp_path("blastx_batches_query.fna");
        std::fs::write(&query_path, query).expect("write query");
        let mut subject =
            std::fs::read(Inputs::path("tests/fasta/MjeNMV.faa")).expect("subject fixture");
        subject.extend_from_slice(b">no letters\n");
        subject.extend_from_slice(
            &std::fs::read(Inputs::path("tests/fixtures/fasta_input/base.s.prot.faa"))
                .expect("fixture"),
        );
        let subject_path = temp_path("blastx_batches_subject.faa");
        std::fs::write(&subject_path, subject).expect("write subject");
        Self {
            query: query_path,
            subject: subject_path,
        }
    }

    fn inputs(&self) -> Inputs {
        Inputs::new(&self.query, &self.subject)
    }
}

impl Drop for Batches {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.query);
        let _ = std::fs::remove_file(&self.subject);
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-297
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
// ...
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ...
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
// ...
//         formatter.PrintEpilog(opt);
// ```
// The half of the test of this name in `LOSAT/tests/run_local_blastx.rs` that the api does not
// reach: the observer must report the exact bytes of every HSP in every format (the CLI half
// there compares the reports and the warnings with the CLI runs). The reference is the entry's
// own hit records.
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let batches = Batches::new();
    // (inputs, whether they have hits)
    let cases = [
        (
            Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa"),
            true,
        ),
        (batches.inputs(), true),
        (
            Inputs::new("tests/fasta/small_test.fasta", "tests/fasta/PajaWSV.faa"),
            false,
        ),
    ];
    for (inputs, hits) in &cases {
        for threads in ["1", "2", "4"] {
            for extra in [&[][..], &["-lcase_masking"][..]] {
                let extra = [extra, &["-num_threads", threads]].concat();
                let context = format!("{} {extra:?}", inputs.query);
                assert_observer_ranges_of(inputs, &extra, *hits, &context);
            }
        }
    }
    // The observer without the hit records numbers the HSPs of later batches on from those of
    // the earlier ones, as with them; outfmt 0 need not come first.
    let inputs = batches.inputs();
    let tabular = inputs.run_with(&["6", "7"], &["-num_threads", "2"], true, false);
    assert!(tabular.hits.is_empty());
    let rows = assert_tabular_ranges(&tabular, 0, Some(1), "observer without hits");
    let reordered = inputs.run(&["7", "0", "6"], &["-num_threads", "2"], true);
    assert_observer_ranges(&reordered, [1, 2, 0], "outfmt 7, 0, 6");
    assert_eq!(rows, reordered.hits.len());
    assert!(tabular.ranges == [reordered.ranges[2].clone(), reordered.ranges[0].clone()]);
}

/// One `run_local` call of outfmt 0, 6 and 7: for inputs with hits, its observer and hits agree
/// with the reports; else it has no hits and no ranges.
fn assert_observer_ranges_of(inputs: &Inputs, extra: &[&str], hits: bool, context: &str) {
    let together = inputs.run(&FORMATS, extra, true);
    if hits {
        assert_observer_ranges(&together, [0, 1, 2], context);
    } else {
        assert!(together.hits.is_empty(), "{context}: no hits");
        assert!(
            together.ranges.iter().all(Vec::is_empty),
            "{context}: no ranges"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973
// ```c++
//     if (args.Exist(kArgNumDescriptions) && args[kArgNumDescriptions]) {
//        m_NumDescriptions = args[kArgNumDescriptions].AsInteger();
//     }
// ```
// The half of the test of this name in `LOSAT/tests/run_local_blastx.rs` that the api does not
// reach: the entry refuses formats that need different hit list sizes (the api splits them into
// searches itself), and the formats searched by separate calls share one pool.
#[test]
fn display_counts_of_the_pairwise_report() {
    let inputs = Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa");
    let extra = [
        "-num_descriptions",
        "3",
        "-num_alignments",
        "2",
        "-num_threads",
        "2",
    ];
    let queries = queries(&Inputs::path(&inputs.query));
    let subjects = subjects(&Inputs::path(&inputs.subject));
    let (refused, outputs, diagnostics) = try_run(&inputs, &queries, &subjects, &FORMATS, &extra);
    let error = refused.expect_err("formats with different hit list sizes are refused");
    assert!(
        error.to_string().contains("different hit list sizes"),
        "{error}"
    );
    assert!(outputs.iter().all(Vec::is_empty) && diagnostics.is_empty());
    let (pairwise, tabular) = with_run_pool(2, "BLASTX", || {
        (
            inputs.run(&["0"], &extra, false),
            inputs.run(&["6", "7"], &extra, false),
        )
    });
    assert!(!pairwise.outputs[0].is_empty(), "outfmt 0");
    assert!(
        tabular.outputs.iter().all(|output| !output.is_empty()),
        "outfmt 6 and 7"
    );
    assert!(
        String::from_utf8_lossy(&tabular.diagnostics).contains("-num_descriptions is ignored"),
        "the tabular formats' warning"
    );
}

/// The side files of outfmt 13 and 14 without a file name are refused, as are the ids of the
/// deflines of records that kept none, before anything is written.
#[test]
fn refusals_come_before_any_output() {
    let inputs = Inputs::new("tests/fasta/small_test.fasta", "tests/fasta/PajaWSV.faa");
    for (formats, extra) in [
        (&["14"][..], &[][..]),
        (&["6"][..], &["-parse_deflines"][..]),
    ] {
        let queries = queries(&Inputs::path(&inputs.query));
        let subjects = subjects(&Inputs::path(&inputs.subject));
        let (result, outputs, diagnostics) = try_run(&inputs, &queries, &subjects, formats, extra);
        let error = result.expect_err("refused");
        assert!(
            error.to_string().contains("not supported in"),
            "{formats:?} {extra:?}: {error}"
        );
        assert!(outputs.iter().all(Vec::is_empty) && diagnostics.is_empty());
    }
}

/// `check_options` (LOSAT Web's `validate`) fails where a run of the same formats fails before it
/// searches, with the run's error: LOSAT Web's refusals, LOSAT's BLASTX limits, formats that need
/// different hit lists, and NCBI's errors of the options, among them the hit list size of each
/// format (`-num_alignments 0` empties the tabular formats' list only).
#[test]
fn check_options_fails_as_a_run_does() {
    let inputs = Inputs::new("tests/fasta/small_test.fasta", "tests/fasta/PajaWSV.faa");
    let queries = queries(&Inputs::path(&inputs.query));
    let subjects = subjects(&Inputs::path(&inputs.subject));
    let mut failed = 0;
    for (formats, extra) in [
        (&FORMATS[..], &[][..]),
        (&["11"][..], &[][..]),
        // NCBI's formatting handler: a custom delimiter that a field uses (W10).
        (&["6 delim=; qseqid staxids", "0"][..], &[][..]),
        (&["6"][..], &["-parse_deflines"][..]),
        // NCBI's error of a matrix BLASTX does not take (IDENTITY, `std_matrix_only`).
        (&FORMATS[..], &["-matrix", "IDENTITY"][..]),
        (&FORMATS[..], &["-evalue", "0"][..]),
        (&FORMATS[..], &["-ungapped"][..]),
        (&["6", "7"][..], &["-num_alignments", "0"][..]),
        (&["0"][..], &["-num_alignments", "0"][..]),
        (&["0", "6"][..], &["-num_alignments", "3"][..]),
        // The ranges' syntax where the command line reads them: the subject range before the
        // query range, both after LOSAT's limits and before NCBI's errors of the options.
        (&FORMATS[..], &["-query_loc", "2-1"][..]),
        (&FORMATS[..], &["-subject_loc", "abc"][..]),
        (
            &FORMATS[..],
            &["-query_loc", "2-1", "-subject_loc", "5-1"][..],
        ),
        (
            &FORMATS[..],
            &["-subject_loc", "2-1", "-matrix", "IDENTITY"][..],
        ),
        (&FORMATS[..], &["-query_loc", "2-1", "-evalue", "0"][..]),
        (
            &FORMATS[..],
            &["-query_loc", "1-100", "-subject_loc", "2-50"][..],
        ),
    ] {
        let args = blastx_args(&inputs.query, &inputs.subject, extra);
        let checked = crate::algorithm::blastx::check_options(&args, formats)
            .map_err(|error| format!("{error:#}"));
        let (result, _, _) = try_run(&inputs, &queries, &subjects, formats, extra);
        let run = result.map_err(|error| format!("{error:#}"));
        assert_eq!(checked, run, "{formats:?} {extra:?}");
        failed += usize::from(run.is_err());
    }
    // Every case but the three searches and the valid ranges fails.
    assert_eq!(failed, 12);
}

/// Values that the command line's parser refuses and a caller of `run_local` can set itself get
/// NCBI's errors of the options where NCBI's core checks them (blast_options.c:1312-1348 for the
/// word size, :1510-1516 for the hit list size; blast_aux.cpp:662-673 for a task NCBI does not
/// know, read in lower case), else LOSAT's refusal of what it does not run, instead of a search.
#[test]
fn hand_set_values_are_refused_as_ncbis_api_or_losat_refuses_them() {
    let inputs = Inputs::new("tests/fasta/small_test.fasta", "tests/fasta/PajaWSV.faa");
    let queries = queries(&Inputs::path(&inputs.query));
    let subjects = subjects(&Inputs::path(&inputs.subject));
    let cases: [(&str, fn(&mut BlastxArgs)); 9] = [
        ("Word-size must be greater than zero", |a| {
            a.word_size = Some(0)
        }),
        (
            "Word-size must be less than 8 for a tblastn, blastp or blastx search",
            |a| a.word_size = Some(8),
        ),
        ("unsupported BLASTX word_size: use 3 or 5", |a| {
            a.word_size = Some(4)
        }),
        ("No hits are being saved", |a| a.max_target_seqs = Some(-1)),
        (
            "the task blastp is not supported by LOSAT's BLASTX (use blastx or blastx-fast)",
            |a| a.task = "blastp".into(),
        ),
        ("'blastz' is not a supported task", |a| {
            a.task = "BLASTZ".into()
        }),
        ("Argument \"strand\". Illegal value", |a| {
            a.strand = "up".into()
        }),
        ("-evalue requires a finite real value", |a| {
            a.evalue = f64::INFINITY
        }),
        ("SEG requires finite cutoffs", |a| a.seg.locut = f64::NAN),
    ];
    for (message, set) in cases {
        let mut args = blastx_args(&inputs.query, &inputs.subject, &[]);
        set(&mut args);
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
            run_local_blastx(args, &queries, &subjects, &mut outputs)
        };
        let error = format!("{:#}", result.expect_err(message));
        assert!(error.contains(message), "{message}: {error}");
        assert!(sink.is_empty(), "{message}");
    }
    // The database entry checks the task first too (before the database is opened).
    let mut args = blastx_args(&inputs.query, &inputs.subject, &[]);
    args.db = Some("no_such_database".into());
    args.task = "blastz".into();
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
        crate::api::local_blast::run_local_database_blastx(
            args,
            &queries,
            crate::blastinput::ncbi_environment::ApplicationSettings::default(),
            &mut outputs,
        )
    };
    let error = format!("{:#}", result.expect_err("refused"));
    assert!(
        error.contains("'blastz' is not a supported task"),
        "{error}"
    );
}

/// A query whose first batch ends on its first record only when the batch counts the whole record:
/// 10,050 letters (9,950 in the range `101-20000`), then 200 (100 in the range) and 50 (skipped).
/// Counting the cut letters would put the first two in one batch and leave the third alone, an
/// empty last batch.
fn whole_record_batch_query() -> PathBuf {
    let genome = fixture_sequence("MjeNMV.fasta");
    let mut query = Vec::new();
    for (name, letters) in [
        ("long", &genome[0..10_050]),
        ("middle", &genome[20_000..20_200]),
        ("short", &genome[30_000..30_050]),
    ] {
        query.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in letters.chunks(70) {
            query.extend_from_slice(line);
            query.push(b'\n');
        }
    }
    let path = temp_path("blastx_whole_record_batch.fna");
    std::fs::write(&path, query).expect("write query");
    path
}

/// A subject file of records with the reader's warnings (`batch_late_warn.prot.faa`) and a record
/// of 10 letters between them, after the first warning.
fn short_subject_between_warnings() -> PathBuf {
    let text = std::fs::read_to_string(Inputs::path(
        "tests/fixtures/fasta_input/batch_late_warn.prot.faa",
    ))
    .expect("fixture");
    let at = text.find(">q5").expect("q5");
    let subject = format!("{}>short subject\nMLCRIFKDHL\n{}", &text[..at], &text[at..]);
    let path = temp_path("blastx_short_subject.faa");
    std::fs::write(&path, subject).expect("write subject");
    path
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-170
// ```c++
//     TSeqPos size_read = 0;
//
//     while (size_read < GetBatchSize()) {
// ...
//         catch (const exception&) {
//             continue; //SB-2307. ignore well formed, not found accession
//         }
// ...
//             size_read += sequence::GetLength(loc->GetInt().GetId(),
//                                              q->GetScope());
// ```
/// The half of the test of this name in `LOSAT/tests/run_local_blastx.rs` that the api does not
/// reach: the observer's ranges are exact for the records cut by `-query_loc` and `-subject_loc`
/// (queries skipped or without letters in the range, subjects without letters, lower-case masks,
/// the minus strand and queries of several batches, at 1, 2 and 4 threads). The CLI half there
/// compares the reports and the warnings with the CLI runs.
#[test]
fn ranges_cut_the_records_as_the_command_line_does() {
    let batches = Batches::new();
    let whole_record_batch = whole_record_batch_query();
    let bslice = || Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa");
    let lowercase_subjects = || {
        Inputs::new(
            "tests/fasta/outfmt0/sml_bx_subj_q.fna",
            "tests/fasta/outfmt0/sml_bx_subj_s.faa",
        )
    };
    // (inputs, ranges and options, whether they have hits)
    let cases = [
        (bslice(), &["-query_loc", "1001-3000"][..], true),
        // The first record ends before the range: it is not searched.
        (bslice(), &["-query_loc", "5200-5600"][..], true),
        // The first record has no letters in the range.
        (bslice(), &["-query_loc", "5001-5600"][..], true),
        (
            bslice(),
            &["-query_loc", "1002-3000", "-strand", "minus"][..],
            true,
        ),
        (bslice(), &["-subject_loc", "10-70"][..], true),
        // Subjects without letters in the range.
        (bslice(), &["-subject_loc", "77-300"][..], true),
        (
            bslice(),
            &["-query_loc", "1001-3000", "-subject_loc", "10-70"][..],
            true,
        ),
        (
            Inputs::new(
                "tests/fasta/outfmt0/x5_bx_q_lc_mix.fna",
                "tests/db_search/bslice.faa",
            ),
            &["-lcase_masking", "-query_loc", "501-2800"][..],
            true,
        ),
        // Several reader batches: the records without letters and the short query are skipped,
        // the batch counts the long record whole.
        (batches.inputs(), &["-query_loc", "50-10050"][..], true),
        // The batch counts the whole record: the first record fills the first batch.
        (
            Inputs::new(&whole_record_batch, "tests/fasta/MjeNMV.faa"),
            &["-query_loc", "101-20000"][..],
            true,
        ),
        // Lower-case runs in the subjects, cut by the subject range.
        (
            lowercase_subjects(),
            &["-lcase_masking", "-subject_loc", "50-150"][..],
            true,
        ),
        (
            lowercase_subjects(),
            &[
                "-lcase_masking",
                "-subject_loc",
                "100-300",
                "-query_loc",
                "1-900",
            ][..],
            true,
        ),
        (
            batches.inputs(),
            &["-query_loc", "50-10050", "-lcase_masking"][..],
            true,
        ),
    ];
    for (inputs, ranges, hits) in &cases {
        for threads in ["1", "2", "4"] {
            let extra = [ranges, &["-num_threads", threads][..]].concat();
            let context = format!("{} {extra:?}", inputs.query);
            assert_observer_ranges_of(inputs, &extra, *hits, &context);
        }
    }
    let _ = std::fs::remove_file(&whole_record_batch);
}

/// The half of the test of this name in `LOSAT/tests/run_local_blastx.rs` that the api does not
/// reach: the entry's error for the ranges and the records is NCBI's typed `NativeError` with the
/// command line's message (the other half compares the output, the warnings and the message with
/// the CLI runs).
#[test]
fn errors_of_the_ranges_follow_the_command_line() {
    let batches = Batches::new();
    let short_subject = short_subject_between_warnings();
    let cases = [
        (
            Inputs::new("tests/db_search/bslice_q.fna", &short_subject),
            &["-subject_loc", "30-50"][..],
            "Invalid from coordinate",
        ),
        (
            Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa"),
            &["-query_loc", "6000-7000"][..],
            "Empty CBlastQueryVector",
        ),
        // The long record fills the first batch; the later records are all skipped.
        (
            batches.inputs(),
            &["-query_loc", "10001-10050"][..],
            "Empty CBlastQueryVector",
        ),
    ];
    for (inputs, ranges, message) in &cases {
        for threads in ["1", "2", "4"] {
            let extra = [ranges, &["-num_threads", threads][..]].concat();
            let context = format!("{} {extra:?}", inputs.query);
            let queries = queries(&Inputs::path(&inputs.query));
            let subjects = subjects(&Inputs::path(&inputs.subject));
            let (result, _, _) = try_run(inputs, &queries, &subjects, &FORMATS, &extra);
            let error = result.expect_err(&context);
            let native = error
                .downcast_ref::<crate::cli::NativeError>()
                .unwrap_or_else(|| panic!("{context}: {error:#}"));
            assert!(native.message.contains(message), "{context}: {error:#}");
        }
    }
    let _ = std::fs::remove_file(&short_subject);
}
