//! Two options of a psiblast run that name one file, with a `.ncbirc` of the byte 0xEF (line RG
//! item F1w-r3): NCBI's start posts `Error: (110.4) Error reading the registry after line 1:` to
//! the file of `-logfile`, which NCBI opens first, so the log has bytes before the PSSM option
//! opens the file. The run's file is the command line's, which is NCBI BLAST+ 2.17.0's (the item's
//! probes): the PSSM file (version 0) removes the log's file and makes another, a numbered PSSM
//! file empties it, so the start's message is gone. A test process of its own: a registered
//! `.ncbirc` is read by every run of the process.

use losat_web_adapter::{emit, files, run, store, zip};
use sha2::{Digest, Sha256};

fn fasta(name: &str) -> Vec<u8> {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT/tests/fasta");
    std::fs::read(format!("{root}/{name}")).unwrap()
}

fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// The files of the run (stream 8) of psiblast's PSSM probe with `words`, on `query`: each name,
/// the option of its entry's comment, its size and SHA-256.
fn run_files(words: &[&str], query: &[u8]) -> Vec<(String, String, usize, String)> {
    let mut argv = vec![
        "psiblast",
        "-query",
        "psiblast/q1b_gi.faa",
        "-subject",
        "psiblast/q1_subj_ids.faa",
        "-save_pssm_after_last_round",
        "-num_descriptions",
        "3",
        "-num_alignments",
        "2",
    ];
    argv.extend_from_slice(words);
    argv.extend(["-num_threads", "1"]);
    let program = run::Program::parse("psiblast").unwrap();
    let (query, _) = store::register("psiblast", 0, true, false, query).unwrap();
    let (subject, _) = store::register(
        "psiblast",
        1,
        true,
        false,
        &fasta("psiblast/q1_subj_ids.faa"),
    )
    .unwrap();
    emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = store::with_inputs(
        program,
        query,
        subject,
        &argv,
        Some(true),
        |queries, end, subjects| run::run(&argv, queries, end, subjects),
    );
    store::release(query).unwrap();
    store::release(subject).unwrap();
    result
        .unwrap()
        .unwrap_or_else(|error| panic!("{words:?}: {error}"));
    let mut archive = Vec::new();
    emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            if *stream == emit::STREAM_RUN_FILES {
                archive.extend_from_slice(bytes);
            }
        }
    });
    zip::read_with_comments(&archive)
        .into_iter()
        .map(|(name, bytes, comment)| {
            (
                String::from_utf8(name).unwrap(),
                String::from_utf8(comment).unwrap(),
                bytes.len(),
                sha256(&bytes),
            )
        })
        .collect()
}

const PSSM: (usize, &str) = (
    39_579,
    "6f6879f8f1c6bbf7b67dbb38f940fd0c88bafb5932dd83e5f61b8b2bc1d788f3",
);
const ASCII_ROUND_2: (usize, &str) = (
    6_005,
    "641993a472c6fafe5c1a5f561e032262157bf0ae94b9e8050f0f5cc8b4dc8f24",
);
const ASCII_ROUND_1: (usize, &str) = (
    6_005,
    "0f0bf9804c7d1c1aa7d17aa374c7dcf70b83dfb212f792e37fd9b052433a3330",
);
const PSSM_TWO_QUERIES: (usize, &str) = (
    39_523,
    "e6f063a1f650e7e5e6131e24ed7840c9a8d7eb18e1ea27f742efd2c13bcb824d",
);
/// The log of the start's post alone (`Error: (110.4) Error reading the registry after line 1: `).
const START_POST: (usize, &str) = (
    57,
    "eb0e621d721702ff333211f96cf83cf390e719b998a0ba67d3ea574cc77a925c",
);

#[test]
fn a_file_that_the_log_and_a_pssm_option_name_is_the_command_lines_file() {
    let ncbirc = files::begin(files::KIND_REGISTRY, ".ncbirc", 1).unwrap();
    files::chunk(ncbirc, &[0xef]).unwrap();
    files::end(ncbirc).unwrap();
    let row = |name: &str, option: &str, (size, sha): (usize, &str)| {
        (name.to_string(), option.to_string(), size, sha.to_string())
    };
    let one = fasta("psiblast/q1b_gi.faa");
    let mut two = one.clone();
    two.extend_from_slice(b">second query");
    two.extend_from_slice(&one[one.iter().position(|&byte| byte == b'\n').unwrap()..]);
    for (query, words, expected) in [
        // The control: distinct names, the log with the start's post.
        (
            &one,
            &["-num_iterations", "2", "-logfile", "L", "-out_pssm", "P"][..],
            vec![row("P", "out_pssm", PSSM), row("L", "logfile", START_POST)],
        ),
        (
            &one,
            &["-num_iterations", "2", "-logfile", "L", "-out_pssm", "L"],
            vec![row("L", "out_pssm", PSSM)],
        ),
        (
            &one,
            &[
                "-num_iterations",
                "2",
                "-logfile",
                "L",
                "-out_ascii_pssm",
                "L",
            ],
            vec![row("L", "out_ascii_pssm", ASCII_ROUND_2)],
        ),
        (
            &one,
            &[
                "-num_iterations",
                "2",
                "-logfile",
                "E.1",
                "-out_pssm",
                "E",
                "-save_each_pssm",
            ],
            vec![row("E.1", "out_pssm", PSSM), row("E.2", "out_pssm", PSSM)],
        ),
        (
            &one,
            &[
                "-num_iterations",
                "2",
                "-logfile",
                "E.2",
                "-out_ascii_pssm",
                "E",
                "-save_each_pssm",
            ],
            vec![
                row("E.1", "out_ascii_pssm", ASCII_ROUND_1),
                row("E.2", "out_ascii_pssm", ASCII_ROUND_2),
            ],
        ),
        (
            &two,
            &["-num_iterations", "1", "-logfile", "L", "-out_pssm", "L"],
            vec![row("L", "out_pssm", PSSM_TWO_QUERIES)],
        ),
    ] {
        assert_eq!(run_files(words, query), expected, "{words:?}");
    }
    files::release(ncbirc).unwrap();
}
