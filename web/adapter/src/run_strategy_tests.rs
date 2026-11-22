//! The search strategy file of a run (`-export_search_strategy`, line RG item F1w): a file of the
//! run (stream 8) under the name as typed, with the bytes that the command line writes to the file
//! for the argv, which are NCBI BLAST+ 2.17.0's (`LOSAT/tests/fixtures/strategy`, the rows of
//! `LOSAT/tests/export_search_strategy.rs`, in LOSAT's words: `losat_output_bytes`), and the command
//! line's errors and refusals of the name.

use super::*;
use std::collections::BTreeMap;

fn fasta(name: &str) -> Vec<u8> {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT/tests/fasta");
    std::fs::read(format!("{root}/{name}")).unwrap()
}

/// NCBI's strategy file of the row `name` of `LOSAT/tests/fixtures/strategy/cases.tsv`, in LOSAT's
/// words.
fn frozen(name: &str) -> Vec<u8> {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../LOSAT/tests/fixtures/strategy"
    );
    losat::web_adapter::losat_output_bytes(&std::fs::read(format!("{root}/{name}.ss")).unwrap())
        .into_owned()
}

/// The result of `run` of `words` (program name first) on the FASTA files `query` and `subject`
/// (of `LOSAT/tests/fasta`), and the bytes of each stream that it emitted, joined.
fn run_files_of(
    words: &[&str],
    query: &str,
    subject: &str,
) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let program = Program::parse(words[0]).unwrap();
    let (query, _) = crate::store::register(words[0], 0, true, &fasta(query)).unwrap();
    let (subject, _) = crate::store::register(words[0], 1, true, &fasta(subject)).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = crate::store::with_inputs(
        program,
        query,
        subject,
        words,
        Some(true),
        |queries, end, subjects| run(words, queries, end, subjects),
    );
    crate::store::release(query).unwrap();
    crate::store::release(subject).unwrap();
    let mut streams = BTreeMap::new();
    crate::emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            streams
                .entry(*stream)
                .or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        }
    });
    (result.unwrap(), streams)
}

/// The files of the run (stream 8), by name.
fn files(streams: &BTreeMap<u32, Vec<u8>>) -> Vec<(String, Vec<u8>)> {
    streams
        .get(&STREAM_RUN_FILES)
        .map(|archive| {
            crate::zip::read(archive)
                .into_iter()
                .map(|(name, bytes)| (String::from_utf8(name).unwrap(), bytes))
                .collect()
        })
        .unwrap_or_default()
}

/// The rows of `cases.tsv` with `-subject` that a browser's run takes as it is (no environment):
/// the row, the program, its query and subject (of `LOSAT/tests/fasta`) and its other words.
const ROWS: [(&str, &str, &str, &str, &[&str]); 8] = [
    (
        "blastp.sub",
        "blastp",
        "outfmt11/arch_bp_q.faa",
        "outfmt11/arch_bp_s.faa",
        &["-outfmt", "6"],
    ),
    (
        "blastp.sub.nd",
        "blastp",
        "outfmt11/arch_bp_q.faa",
        "outfmt11/arch_bp_s.faa",
        &["-num_descriptions", "3", "-num_alignments", "2"],
    ),
    (
        "blastp.sub.qloc",
        "blastp",
        "outfmt11/arch_bp_q.faa",
        "outfmt11/arch_bp_s.faa",
        &["-query_loc", "10-100", "-outfmt", "6"],
    ),
    (
        "blastn.sub",
        "blastn",
        "outfmt11/arch_q.fna",
        "outfmt11/arch_s.fna",
        &["-outfmt", "6"],
    ),
    (
        "blastn.sub.lcase",
        "blastn",
        "outfmt11/arch_lcout.fna",
        "outfmt11/arch_s.fna",
        &["-lcase_masking", "-outfmt", "6"],
    ),
    (
        "blastx.sub",
        "blastx",
        "outfmt11/arch_bx_q.fna",
        "outfmt11/arch_bx_s.faa",
        &["-outfmt", "6"],
    ),
    (
        "tblastn.sub",
        "tblastn",
        "outfmt11/arch_tbn_q.faa",
        "outfmt11/arch_tbn_s.fna",
        &["-outfmt", "6"],
    ),
    (
        "tblastx.sub",
        "tblastx",
        "outfmt11/arch_tbx_q.fna",
        "outfmt11/arch_tbx_s.fna",
        &["-outfmt", "6"],
    ),
];

/// The argv of a row with `-export_search_strategy NAME`.
fn row_argv<'a>(
    row: &(&'a str, &'a str, &'a str, &'a str, &[&'a str]),
    name: &'a str,
) -> Vec<&'a str> {
    let (_, program, query, subject, words) = *row;
    let mut argv = vec![program, "-query", query, "-subject", subject];
    argv.extend_from_slice(words);
    argv.extend(["-export_search_strategy", name]);
    argv
}

/// Each row's file is a file of the run under the name as typed, with NCBI's bytes, whatever
/// other searches the run makes (`blastp.sub.nd`: outfmt 0's hit list of 3 and outfmt 6's of 500
/// are two searches, and the file is the argv's, outfmt 0's), next to the file of `-logfile`, and
/// the outputs are those of the run without the option. The start's JSON names the file
/// (`strategy_file`).
#[test]
fn the_strategy_file_is_a_file_of_the_run_with_the_command_lines_bytes() {
    for row in &ROWS {
        let (name, _, query, subject, _) = *row;
        for file in ["s.asn", "out/s.asn", "/tmp/s.asn"] {
            let words = row_argv(row, file);
            validate(&words).unwrap();
            assert!(
                start_json(&words)
                    .unwrap()
                    .contains(&format!("\"strategy_file\":\"{file}\"")),
                "{name}"
            );
            let (result, streams) = run_files_of(&words, query, subject);
            result.unwrap();
            assert_eq!(
                files(&streams),
                [(file.to_string(), frozen(name))],
                "{name} {file}"
            );
            let plain: Vec<&str> = words[..words.len() - 2].to_vec();
            let (_, without) = run_files_of(&plain, query, subject);
            for stream in [0, 6, 7, 1, STREAM_DIAGNOSTICS] {
                assert_eq!(
                    streams.get(&stream),
                    without.get(&stream),
                    "{name} {stream}"
                );
            }
        }
        let mut words = row_argv(row, "s.asn");
        words.extend(["-logfile", "run.log"]);
        let (result, streams) = run_files_of(&words, query, subject);
        result.unwrap();
        let named: Vec<String> = files(&streams).into_iter().map(|(name, _)| name).collect();
        assert_eq!(named, ["s.asn", "run.log"], "{name}");
    }
}

/// A run that fails after the file is made has what the command line's file has then: empty for an
/// error of the options (`blastp.sub.ungapped`), the request of the batch whose search fails
/// (`blastp.sub.ee`: a query without letters).
#[test]
fn a_run_that_fails_after_the_file_is_made_keeps_the_command_lines_file() {
    for (name, query, words) in [
        (
            "blastp.sub.ungapped",
            "outfmt11/arch_bp_q.faa",
            &["-ungapped", "-outfmt", "6"][..],
        ),
        (
            "blastp.sub.ee",
            "strategy/ee_blastp.faa",
            &["-outfmt", "6"][..],
        ),
    ] {
        let row = (name, "blastp", query, "outfmt11/arch_bp_s.faa", words);
        let argv = row_argv(&row, "s.asn");
        let (result, streams) = run_files_of(&argv, query, "outfmt11/arch_bp_s.faa");
        assert!(result.is_err(), "{name}");
        assert_eq!(
            files(&streams),
            [("s.asn".to_string(), frozen(name))],
            "{name}"
        );
    }
}

/// The names that the command line cannot open are NCBI's error, at once (no file); `-` and the
/// standard output are the command line's refusals at NCBI's first request, and the outputs hold
/// nothing; the null device gives no file and the run goes on; the standard error and the file of
/// another option of the run are LOSAT Web's refusals at the same point. A run with an empty query
/// ends before that request: NCBI's run, with no file for `-`.
#[test]
fn the_names_the_run_has_no_file_for() {
    let row = &ROWS[0];
    let (_, _, query, subject, _) = *row;
    for name in ["/etc/s.asn", ".", "out/", ""] {
        let words = row_argv(row, name);
        let message = validate(&words).unwrap_err();
        assert_eq!(
            message,
            format!(
                "Command line argument error: Argument \"export_search_strategy\". File is not accessible{}\n",
                match name {
                    "" => String::new(),
                    name => format!(":  `{name}'"),
                }
            ),
            "{name:?}"
        );
        let (result, streams) = run_files_of(&words, query, subject);
        assert_eq!(result.unwrap_err(), message, "{name:?}");
        assert!(files(&streams).is_empty(), "{name:?}");
    }
    for (name, refusal) in [
        ("-", "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's BLASTP; give a file name"),
        ("/dev/stdout", "an -export_search_strategy file that is also the output (-out or the standard output) is not supported by LOSAT's BLASTP: NCBI writes both through separate streams and the file's bytes follow the order of their flushes"),
        ("/dev/stderr", "-export_search_strategy /dev/stderr (a device or a file of the system, where the run has no file) is not supported in the browser; give a file name of the run's own"),
        ("outfmt11/arch_bp_q.faa", "-export_search_strategy outfmt11/arch_bp_q.faa (the file of -query too, whose bytes the command line's strategy file replaces or mixes with its own) is not supported in the browser; give a file name of the run's own"),
    ] {
        let words = row_argv(row, name);
        validate(&words).unwrap();
        let (result, streams) = run_files_of(&words, query, subject);
        assert_eq!(result.unwrap_err(), refusal, "{name}");
        assert!(files(&streams).is_empty(), "{name}");
        for stream in [0, 6, 7] {
            assert!(streams.get(&stream).is_none_or(Vec::is_empty), "{name} {stream}");
        }
    }
    let words = row_argv(row, "/dev/null");
    let (result, streams) = run_files_of(&words, query, subject);
    result.unwrap();
    assert!(files(&streams).is_empty(), "/dev/null");
    assert!(!streams[&6].is_empty());
    let words = row_argv(row, "-");
    let (result, streams) = run_files_of(&words, "outfmt11/arch_bp_empty.faa", subject);
    result.unwrap();
    assert!(files(&streams).is_empty());
    assert_eq!(
        String::from_utf8_lossy(&streams[&STREAM_DIAGNOSTICS]),
        "Warning: [blastp] Query is Empty!\n"
    );
}
