//! The shared BLASTP entry `run_local` in process: what the entry itself does, not what the
//! command line writes. The tests that compare the entry's reports with the `losat` binary go
//! through `LOSAT::api` (`LOSAT/tests/run_local_blastp.rs`, `LOSAT/tests/blastp_row_windows.rs`).

use std::path::{Path, PathBuf};

use super::support::{
    assert_observer_ranges, field, read_records, reader_records, run_formats, run_formats_with,
    temp_path, Run, TempFasta,
};
use crate::algorithm::blastp::BlastpArgs;
use crate::api::local_blast::{run_local_blastp, FormatOutput, OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use crate::cli::{try_parse_from, Cli, Commands};

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

fn blastp_args(query: &Path, subject: &Path, extra: &[&str]) -> BlastpArgs {
    let mut argv = vec![
        "losat".to_string(),
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

/// The query records of a FASTA file as NCBI's reader reads them for BLASTP.
fn queries(path: &Path) -> Vec<FastaRecord> {
    reader_records(
        path,
        ReaderConfig::query("BLASTP", true, true).with_defline(true),
    )
}

/// The subject records of a FASTA file as NCBI's reader reads them for BLASTP.
fn subjects(path: &Path) -> Vec<FastaRecord> {
    reader_records(
        path,
        ReaderConfig::subject("BLASTP", true, true).with_defline(true),
    )
}

fn run(inputs: &Inputs, formats: &[&str], extra: &[&str], observe: bool) -> Run {
    let queries = queries(&inputs.query);
    let subjects = subjects(&inputs.subject);
    run_formats(formats, observe, |outputs| {
        let args = blastp_args(&inputs.query, &inputs.subject, extra);
        run_local_blastp(args, &queries, &subjects, "", "", outputs)
    })
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
    let queries = queries(&inputs.query);
    let subjects = subjects(&inputs.subject);
    // NCBI blast_args.cpp:2845-2851 and tabular.cpp:70-99: a custom specification of
    // outfmt 0 and a token that is not a field name are ignored ("0 qseqid", "6 nosuchfield"
    // run as "0" and "6"); a format that LOSAT's BLASTP does not write is rejected, and a
    // delimiter that a field uses between its entries is NCBI's error (blast_args.cpp:2753-2792).
    for outfmt in ["1", "6 delim=; qseqid staxids", "abc"] {
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

// NCBI derives the hit list size per output format (blast_args.cpp:2894-2973): with
// `-num_descriptions 3` it is 250 for outfmt 0 and 500 for outfmt 6. One run searches once,
// so formats that need different hit list sizes are refused before anything is written
// (`blastinput/report_counts.rs`), while formats of one hit list size still share a run.
#[test]
fn formats_needing_different_hit_list_sizes_are_refused() {
    let inputs = Inputs::new(8);
    let queries = queries(&inputs.query);
    let subjects = subjects(&inputs.subject);
    let attempt = |formats: &[&str], extra: &[&str]| {
        let mut sinks: Vec<Vec<u8>> = formats.iter().map(|_| Vec::new()).collect();
        let mut diagnostics = Vec::new();
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
        let result = run_local_blastp(
            blastp_args(&inputs.query, &inputs.subject, extra),
            &queries,
            &subjects,
            "",
            "",
            &mut outputs,
        );
        drop(outputs);
        (result, sinks, diagnostics)
    };
    for (formats, extra) in [
        (&["0", "6"][..], &["-num_descriptions", "3"][..]),
        (&["6", "0"][..], &["-num_descriptions", "3"][..]),
        (&["0", "6"][..], &["-num_alignments", "3"][..]),
        (
            &["0", "6", "7"][..],
            &["-num_descriptions", "3", "-num_alignments", "2"][..],
        ),
    ] {
        let (result, sinks, _) = attempt(formats, extra);
        let error = format!("{:#}", result.expect_err("different hit list sizes"));
        assert!(
            error.contains("need different hit list sizes"),
            "{formats:?} {extra:?}: {error}"
        );
        // The option warnings of the first format may precede the refusal (BLASTP's
        // `check_options`, as the CLI writes a warning before an options error), so the
        // diagnostics are not checked; no report bytes are written.
        assert!(
            sinks.iter().all(Vec::is_empty),
            "{formats:?} {extra:?}: nothing is written"
        );
    }
    // The same formats without these options, or with one hit list size, run together.
    for extra in [&[][..], &["-max_target_seqs", "3"][..]] {
        let (result, sinks, _) = attempt(&["0", "6"], extra);
        result.unwrap_or_else(|error| panic!("{extra:?}: {error:#}"));
        assert!(
            !sinks[0].is_empty(),
            "{extra:?}: outfmt 0 writes its header"
        );
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

/// `run_local` with one writer per format: the outputs, the diagnostics and the result.
fn run_local_outputs(
    query: &Path,
    subject: &Path,
    formats: &[&str],
) -> (Vec<Vec<u8>>, Vec<u8>, anyhow::Result<()>) {
    let queries = queries(query);
    let subjects = subjects(subject);
    let mut sinks: Vec<Vec<u8>> = vec![Vec::new(); formats.len()];
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
        run_local_blastp(
            blastp_args(query, subject, &[]),
            &queries,
            &subjects,
            "",
            "",
            &mut outputs,
        )
    };
    (sinks, diagnostics, result)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:649-652
// ```c
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// A batch of queries without residues only stops NCBI after the prolog, with one message
// per query (oracle `blastp.rec_empty_all.q.*`). This half checks the entry's own error (NCBI's
// error with exit code 3 and its message), the prolog it wrote and that no other report has a
// byte; the command line's exit status, standard output and standard error of the same run, and
// the stop after the reports of the batches before (BATCH_SIZE), are checked against the binary in
// `LOSAT/tests/run_local_blastp.rs`.
#[test]
fn a_batch_of_queries_without_residues_stops() {
    let s = protein_sequences("PajaWSV.faa", 2);
    let subject = TempFasta::new(
        "blastp_no_data_subject.faa",
        &[("s1", &s[0]), ("s2", &s[1])],
    );
    let query = TempFasta::new("blastp_no_data_query.faa", &[("q1", b""), ("q2", b"")]);
    let (outputs, diagnostics, result) = run_local_outputs(&query.0, &subject.0, &["0", "6"]);
    let error = result.expect_err("the batch fails");
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .unwrap_or_else(|| panic!("NCBI's error: {error:#}"));
    assert_eq!(native.exit, 3);
    assert_eq!(
        native.message,
        "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n"
    );
    assert!(diagnostics.is_empty());
    let prolog = String::from_utf8_lossy(&outputs[0]);
    assert!(
        prolog.starts_with("BLASTP 2.17.0+\n") && prolog.ends_with(" total letters\n\n"),
        "{prolog}"
    );
    assert!(outputs[1].is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
// ```c++
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
// ```
// The rendered rows of blastp's final hits made a window of queries at a time
// (`src/algorithm/blastp/blast_engine/query_rows.rs`): outfmt 0, BLAST XML and the tabular fields
// that show the alignment write the bytes of the rows made for every query up front.
//
// The rows of every query are made up front only for the caller's hit records and for a format
// observer of outfmt 0; every other call makes them per window (`query_rows::row_windows`, one
// thread's window budget of 1 MiB of estimated report, two threads' 2 MiB). The input holds
// 15 queries, 12 of which each hit 180 subjects of about 450 letters, about 5 MB of estimated
// report, so a search makes four windows at one thread and two at two threads; three more queries
// (the letters reversed) have no hits, one of them the last. This half is the entry's own contract,
// in process: the same `run_local_blastp` call with the hit records (rows up front), with
// neither hit records nor observer, with the observer alone, and with the observer without
// outfmt 0 write the same bytes. The CLI run of each format against these bytes, through
// `LOSAT::api`, and the CLI's search in blocks of one query batch each against its search of all
// queries, are checked in `LOSAT/tests/blastp_row_windows.rs`.
// Two queries hold a low-complexity run that SEG masks (`-seg yes`: the lowercase query rows of
// outfmt 0), and some subjects hold lower-case letters (`-lcase_masking`: the lowercase subject
// rows).
#[test]
fn rows_made_per_window_write_the_bytes_of_rows_made_up_front() {
    let inputs = WindowInputs::new();
    let queries = reader_records(&inputs.query.0, ReaderConfig::query("BLASTP", true, true));
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTP", true, true),
    );
    for threads in [1, 2] {
        let args = || {
            let options = window_options(threads);
            let options: Vec<&str> = options.iter().map(String::as_str).collect();
            blastp_args(&inputs.query.0, &inputs.subject.0, &options)
        };
        let run = |observe: bool, hits: bool| {
            run_formats_with(&WINDOW_FORMATS, observe, hits, |outputs| {
                run_local_blastp(args(), &queries, &subjects, "", "", outputs)
            })
        };
        let up_front = run(false, true);
        // The windows: the estimated report is more than twice the budget, so the rows come in
        // two windows or more.
        let mut by_query = vec![Vec::new(); queries.len()];
        for hit in &up_front.hits {
            by_query[hit.hit.q_idx as usize].push(hit.hit.length);
        }
        let estimate: usize = by_query
            .into_iter()
            .map(|lengths| query_estimate(lengths.into_iter()))
            .sum();
        assert!(
            estimate > 2 * window_budget(threads),
            "n{threads}: estimated report {estimate} B makes one window"
        );
        let text = String::from_utf8_lossy(&up_front.outputs[0]).into_owned();
        assert!(
            text.contains("Query= query14 copy 14"),
            "n{threads}: outfmt 0 reaches the last query"
        );
        assert_eq!(
            text.matches("***** No hits found *****").count(),
            NO_HITS.len(),
            "n{threads}: the queries without hits"
        );
        // The letters of the alignment rows that start with `label`.
        let row_letters = |label: &str| -> Vec<String> {
            text.lines()
                .filter(|line| line.starts_with(label))
                .filter_map(|line| line.split_whitespace().nth(2).map(str::to_string))
                .collect()
        };
        assert!(
            row_letters("Query ")
                .iter()
                .any(|row| row.contains("qqqqq")),
            "n{threads}: a SEG mask shows in a query row"
        );
        assert!(
            row_letters("Sbjct ")
                .iter()
                .any(|row| row.bytes().any(|byte| byte.is_ascii_lowercase())),
            "n{threads}: a lower-case subject run shows in a subject row"
        );
        for (name, other) in [
            ("windows", run(false, false)),
            ("observer", run(true, false)),
        ] {
            for (index, outfmt) in WINDOW_FORMATS.iter().enumerate() {
                assert!(
                    other.outputs[index] == up_front.outputs[index],
                    "n{threads} {name}: outfmt {outfmt:?} differs from the rows made up front"
                );
            }
            assert_eq!(other.diagnostics, up_front.diagnostics, "n{threads} {name}");
        }
        // An observer without outfmt 0: the rows of the other formats per window, the observer's
        // HSP numbers over the whole report.
        let observed = run_formats_with(&WINDOW_FORMATS[1..], true, false, |outputs| {
            run_local_blastp(args(), &queries, &subjects, "", "", outputs)
        });
        for (index, outfmt) in WINDOW_FORMATS[1..].iter().enumerate() {
            assert!(
                observed.outputs[index] == up_front.outputs[index + 1],
                "n{threads} observer without outfmt 0: outfmt {outfmt:?} differs"
            );
        }
    }
}

const WINDOW_FORMATS: [&str; 5] = ["0", "5", CUSTOM, "6", "16"];

/// The window budget of `x_pairwise_par::window_bytes` at `threads` threads.
fn window_budget(threads: usize) -> usize {
    (1 << 20) * threads
}

/// `x_pairwise_par::query_estimate_of` of the rows of one query's HSPs of lengths `lengths`.
fn query_estimate(lengths: impl Iterator<Item = usize>) -> usize {
    lengths.fold(512, |bytes, length| bytes + 100 + 400 + 4 * length)
}

/// `sequence` with the letter at `position` replaced by `letter`.
fn with_letter(sequence: &[u8], position: usize, letter: u8) -> Vec<u8> {
    let mut changed = sequence.to_vec();
    changed[position % sequence.len()] = letter;
    changed
}

/// The query and subject files of the windows' search.
struct WindowInputs {
    query: TempFasta,
    subject: TempFasta,
}

impl WindowInputs {
    fn new() -> Self {
        // BDT63507.1 of the fixtures, 449 letters.
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PajaWSV.faa");
        let base = read_records(&fixture)
            .into_iter()
            .find(|record| record.title.starts_with(b"BDT63507.1"))
            .expect("BDT63507.1 in PajaWSV.faa")
            .seq()
            .to_vec();
        assert_eq!(base.len(), 449);
        let subjects: Vec<(String, Vec<u8>)> = (0..180)
            .map(|k| {
                let mut letters = with_letter(&base, 7 * k + 3, b'W');
                // Every fifth subject holds a lower-case run (`-lcase_masking`).
                if k % 5 == 0 {
                    for letter in &mut letters[100..120] {
                        *letter = letter.to_ascii_lowercase();
                    }
                }
                (format!("subj{k}"), letters)
            })
            .collect();
        let queries: Vec<(String, Vec<u8>)> = (0..15)
            .map(|k| {
                let mut letters = with_letter(&base, 11 * k + 5, b'C');
                // Two queries hold a run that SEG masks (`-seg yes`).
                if k == 4 || k == 9 {
                    letters[200..230].fill(b'Q');
                }
                // Three queries without hits (`NO_HITS`).
                if NO_HITS.contains(&k) {
                    letters.reverse();
                }
                (format!("query{k} copy {k}"), letters)
            })
            .collect();
        let file = |name: &str, records: &[(String, Vec<u8>)]| {
            TempFasta::new(
                name,
                &records
                    .iter()
                    .map(|(id, letters)| (id.as_str(), letters.as_slice()))
                    .collect::<Vec<_>>(),
            )
        };
        WindowInputs {
            query: file("rowwin_query.faa", &queries),
            subject: file("rowwin_subject.faa", &subjects),
        }
    }
}

/// The queries whose letters are reversed: no hits at `-evalue 1e-10`.
const NO_HITS: [usize; 3] = [3, 8, 14];

/// The options of every run besides the input files and the format.
fn window_options(threads: usize) -> Vec<String> {
    [
        "-seg",
        "yes",
        "-lcase_masking",
        "-evalue",
        "1e-10",
        "-num_threads",
        &threads.to_string(),
    ]
    .iter()
    .map(|option| option.to_string())
    .collect()
}
