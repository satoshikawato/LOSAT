//! The rendered rows of blastp's final hits made a window of queries at a time
//! (`src/algorithm/blastp/blast_engine/query_rows.rs`): outfmt 0, BLAST XML and the tabular fields
//! that show the alignment write the bytes of the rows made for every query up front.
//!
//! The rows of every query are made up front only for the caller's hit records and for a format
//! observer of outfmt 0; every other call makes them per window (`query_rows::row_windows`, one
//! thread's window budget of 1 MiB of estimated report, two threads' 2 MiB). The input here holds
//! 15 queries, 12 of which each hit 180 subjects of about 450 letters, about 5 MB of estimated
//! report, so a search makes four windows at one thread and two at two threads; three more queries
//! (the letters reversed) have no hits, one of them the last. This file holds what is compared with
//! the binary: the library's searches (`LOSAT::api::run_collect`, which gives the caller's hit
//! records, and `LOSAT::api::run`, which gives neither) write the bytes of the CLI run of each
//! format, and the CLI's search in blocks of one query batch each (`LOSAT_X_BPQBLOCK`) writes the
//! bytes of its search of all queries with the same batches. The calls of the entry
//! `run_local_blastp` with the hit records, with the observer alone and with the observer without
//! outfmt 0 write the same bytes: that is checked in the crate,
//! `LOSAT/src/api/local_blast_tests/blastp.rs`. Two queries hold a low-complexity run that
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

use run_local_support::{read_records, temp_path, Files, TempFasta};
use LOSAT::api::Program;

const CUSTOM: &str =
    "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq btop";
const FORMATS: [&str; 5] = ["0", "5", CUSTOM, "6", "16"];

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
    let files = Files::new(Program::Blastp, &inputs.query.0, &inputs.subject.0);
    for threads in [1, 2] {
        let options = options(threads);
        let options: Vec<&str> = options.iter().map(String::as_str).collect();
        // The library's searches: `run_collect` (the hit records are the caller's), and `run`
        // (neither hit records nor observer), against the CLI run of each format.
        let collected = files.collect(&FORMATS, &options);
        assert!(
            !collected.hsps.is_empty(),
            "n{threads}: the fixture must produce hits"
        );
        let written = files.write(&FORMATS, &options);
        written
            .result
            .unwrap_or_else(|error| panic!("n{threads}: {error}"));
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let cli_output = cli(&inputs, outfmt, threads, &[]);
            assert!(
                collected.reports[index].bytes == cli_output,
                "n{threads} run_collect: outfmt {outfmt:?} differs from the CLI"
            );
            assert!(
                written.outputs[index] == cli_output,
                "n{threads} run: outfmt {outfmt:?} differs from the CLI"
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
