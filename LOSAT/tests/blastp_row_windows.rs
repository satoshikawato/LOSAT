//! The rendered rows of blastp's final hits made a window of queries at a time
//! (`src/algorithm/blastp/blast_engine/query_rows.rs`): outfmt 0, BLAST XML and the tabular fields
//! that show the alignment write the bytes of the rows made for every query up front.
//!
//! The rows of every query are made up front only for the caller's hit records and for a format
//! observer of outfmt 0; every other call makes them per window (`query_rows::row_windows`, one
//! thread's window budget of 1 MiB of estimated report, two threads' 2 MiB). The input here holds
//! 12 queries that each hit 180 subjects of about 450 letters, about 5 MB of estimated report, so
//! a search makes four windows at one thread and two at two threads; three more queries (the
//! letters reversed) have no hits, one of them the last. The same `run_local_blastp` call with the
//! hit records (rows up front), with the observer alone (with and without outfmt 0, which makes the
//! rows up front), and with neither, and the CLI run of each format write the same bytes; so do the
//! CLI's search in blocks of one query batch each (`LOSAT_X_BPQBLOCK`) and its search of all queries
//! with the same batches. Two queries hold a low-complexity run that
//! SEG masks (`-seg yes`: the lowercase query rows of outfmt 0), and some subjects hold lower-case
//! letters (`-lcase_masking`: the lowercase subject rows).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
//! ```c++
//! CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
//! ```

mod run_local_support;

use std::path::Path;
use std::process::Command;

use run_local_support::{run_formats_with, temp_path, TempFasta};
use LOSAT::api::local_blast::run_local_blastp;
use LOSAT::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use LOSAT::cli::{try_parse_from, Cli, Commands};

const CUSTOM: &str =
    "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq btop";
const FORMATS: [&str; 5] = ["0", "5", CUSTOM, "6", "16"];

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
struct Inputs {
    query: TempFasta,
    subject: TempFasta,
}

impl Inputs {
    fn new() -> Self {
        // BDT63507.1 of the fixtures, 449 letters.
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PajaWSV.faa");
        let base = run_local_support::read_records(&fixture)
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
        Inputs {
            query: file("rowwin_query.faa", &queries),
            subject: file("rowwin_subject.faa", &subjects),
        }
    }
}

/// The queries whose letters are reversed: no hits at `-evalue 1e-10`.
const NO_HITS: [usize; 3] = [3, 8, 14];

/// The options of every run besides the input files and the format.
fn options(threads: usize) -> Vec<String> {
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

/// The records of `path` as `run_local_blastp` takes them.
fn records(path: &Path, config: ReaderConfig) -> Vec<FastaRecord> {
    run_local_support::reader_records(path, config)
}

/// The CLI run of `outfmt` with the variables `vars`: its output file.
fn cli(inputs: &Inputs, outfmt: &str, threads: usize, vars: &[(&str, &str)]) -> Vec<u8> {
    let out = temp_path(&format!("rowwin_cli_{}_{threads}.out", outfmt.len()));
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    command.envs(vars.iter().copied());
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(&inputs.query.0)
        .arg("-subject")
        .arg(&inputs.subject.0)
        .args(options(threads))
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-out")
        .arg(&out)
        .output()
        .expect("run LOSAT CLI");
    assert!(output.status.success(), "CLI failed: {output:?}");
    let bytes = std::fs::read(&out).expect("CLI output");
    let _ = std::fs::remove_file(&out);
    bytes
}

#[test]
fn rows_made_per_window_write_the_bytes_of_rows_made_up_front() {
    let inputs = Inputs::new();
    let queries = records(&inputs.query.0, ReaderConfig::query("BLASTP", true, true));
    let subjects = records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTP", true, true),
    );
    for threads in [1, 2] {
        let args = || {
            let query = inputs.query.0.display().to_string();
            let subject = inputs.subject.0.display().to_string();
            let mut argv = vec![
                "LOSAT".to_string(),
                "blastp".to_string(),
                "-query".to_string(),
                query,
                "-subject".to_string(),
                subject,
            ];
            argv.extend(options(threads));
            let cli: Cli = try_parse_from(argv).expect("blastp argv");
            match cli.command {
                Commands::Blastp(args) => args,
                _ => unreachable!("blastp subcommand"),
            }
        };
        let run = |observe: bool, hits: bool| {
            run_formats_with(&FORMATS, observe, hits, |outputs| {
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
            for (index, outfmt) in FORMATS.iter().enumerate() {
                assert!(
                    other.outputs[index] == up_front.outputs[index],
                    "n{threads} {name}: outfmt {outfmt:?} differs from the rows made up front"
                );
            }
            assert_eq!(other.diagnostics, up_front.diagnostics, "n{threads} {name}");
        }
        // An observer without outfmt 0: the rows of the other formats per window, the observer's
        // HSP numbers over the whole report.
        let observed = run_formats_with(&FORMATS[1..], true, false, |outputs| {
            run_local_blastp(args(), &queries, &subjects, "", "", outputs)
        });
        for (index, outfmt) in FORMATS[1..].iter().enumerate() {
            assert!(
                observed.outputs[index] == up_front.outputs[index + 1],
                "n{threads} observer without outfmt 0: outfmt {outfmt:?} differs"
            );
        }
        for (index, outfmt) in FORMATS.iter().enumerate() {
            assert!(
                cli(&inputs, outfmt, threads, &[]) == up_front.outputs[index],
                "n{threads} CLI: outfmt {outfmt:?} differs from the rows made up front"
            );
            // In blocks of one query batch each (three queries of about 450 letters to a batch of
            // 1,500; the composition-based redo reads its batch's queries, so the search of all
            // queries at once with the same batches is the reference).
            let batch = ("BATCH_SIZE", "1500");
            assert!(
                cli(&inputs, outfmt, threads, &[("LOSAT_X_BPQBLOCK", "1"), batch])
                    == cli(&inputs, outfmt, threads, &[("LOSAT_X_BPQBLOCK", "0"), batch]),
                "n{threads} CLI in blocks: outfmt {outfmt:?} differs from the search of all queries"
            );
        }
    }
}
