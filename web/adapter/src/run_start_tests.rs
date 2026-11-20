//! The NCBI application start of a run (line RG item E1w-3a): `-logfile`, `-conffile` and
//! `-dryrun` as the command line takes them (`LOSAT/tests/ncbi_toolkit_file_options.rs`, whose
//! command line arms match NCBI BLAST+ 2.17.0: probes `log_warn`, `conf_ok`, `conf_bad`,
//! `conf_missing` of line RG `E/SURVEY_E1.md`).

use super::*;
use crate::emit::{PLACE_NOWHERE, STREAM_MESSAGES};
use std::collections::BTreeMap;

/// The inputs of each program (`LOSAT/tests/fasta/outfmt11`, as the command line's tests).
const SEARCHES: [(&str, &str, &str); 5] = [
    ("blastp", "arch_bp_q.faa", "arch_bp_s.faa"),
    ("blastx", "arch_bx_q.fna", "arch_bx_s.faa"),
    ("tblastn", "arch_tbn_q.faa", "arch_tbn_s.fna"),
    ("tblastx", "arch_tbx_q.fna", "arch_tbx_s.fna"),
    ("blastn", "arch_tbx_q.fna", "arch_tbx_s.fna"),
];

/// NCBI's report of a registry file with the line `[x`, which ends the program (exit status 2).
const SECTION_REPORT: &str = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";

fn input(name: &str) -> Vec<u8> {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../LOSAT/tests/fasta/outfmt11"
    );
    std::fs::read(format!("{root}/{name}")).unwrap()
}

/// The result of `run` of `words` on the inputs of `program` (`SEARCHES`), and the bytes of each
/// stream that it emitted, joined.
fn run_of(words: &[&str]) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let (_, query, subject) = SEARCHES
        .iter()
        .find(|(program, _, _)| *program == words[0])
        .unwrap();
    run_on(words, &input(query), &input(subject))
}

/// `run_of` on the FASTA records `query` and `subject`.
fn run_on(
    words: &[&str],
    query: &[u8],
    subject: &[u8],
) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let program = Program::parse(words[0]).unwrap();
    let (query, _) = crate::store::register(words[0], 0, query).unwrap();
    let (subject, _) = crate::store::register(words[0], 1, subject).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = crate::store::with_inputs(program, query, subject, words, |queries, subjects| {
        run(words, queries, subjects)
    });
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

/// The argv of `program` on its inputs (`SEARCHES`), then `words`.
fn argv<'a>(program: &'a str, words: &[&'a str]) -> Vec<&'a str> {
    let (_, query, subject) = SEARCHES
        .iter()
        .find(|(name, _, _)| *name == program)
        .unwrap();
    let mut argv = vec![program, "-query", query, "-subject", subject];
    argv.extend_from_slice(words);
    argv
}

/// The text of `stream` (empty when it was not emitted).
fn stream_text(streams: &BTreeMap<u32, Vec<u8>>, stream: u32) -> String {
    String::from_utf8(streams.get(&stream).cloned().unwrap_or_default()).unwrap()
}

/// The files of the run (`STREAM_RUN_FILES`), by name, or `None` when the run gave none.
fn run_files(streams: &BTreeMap<u32, Vec<u8>>) -> Option<Vec<(String, String)>> {
    streams.get(&STREAM_RUN_FILES).map(|archive| {
        crate::zip::read(archive)
            .into_iter()
            .map(|(name, bytes)| {
                (
                    String::from_utf8(name).unwrap(),
                    String::from_utf8(bytes).unwrap(),
                )
            })
            .collect()
    })
}

/// The place of the run's messages that its signal of stream 10 gives (`emit::messages`), or `None`
/// where the run sent none (a start that ended the run). One signal per run (line WQA's sweep, D7,
/// 2026-11-20).
fn messages_place(streams: &BTreeMap<u32, Vec<u8>>) -> Option<u32> {
    streams.get(&STREAM_MESSAGES).map(|bytes| {
        assert_eq!(bytes.len(), 4, "one signal of stream 10");
        u32::from_le_bytes(bytes[..].try_into().unwrap())
    })
}

/// Where `text` is in the run's `streams`: on stream 3, in the file of `-logfile` (stream 8), or
/// nowhere (`PLACE_NOWHERE`); never in both.
fn place_of(streams: &BTreeMap<u32, Vec<u8>>, text: &str) -> u32 {
    let on_3 = stream_text(streams, STREAM_DIAGNOSTICS).contains(text);
    let in_file =
        run_files(streams).is_some_and(|files| files.iter().any(|(_, bytes)| bytes.contains(text)));
    match (on_3, in_file) {
        (true, false) => STREAM_DIAGNOSTICS,
        (false, true) => STREAM_RUN_FILES,
        (false, false) => PLACE_NOWHERE,
        (true, true) => panic!("{text:?} is on stream 3 and in the file of -logfile"),
    }
}

/// The signals of stream 9: each failed stream and the place of its error (`emit::discard`).
fn signals(streams: &BTreeMap<u32, Vec<u8>>) -> Vec<(u32, u32)> {
    let bytes = streams
        .get(&crate::emit::STREAM_DISCARD)
        .cloned()
        .unwrap_or_default();
    assert_eq!(bytes.len() % 8, 0);
    bytes
        .chunks(8)
        .map(|signal| {
            (
                u32::from_le_bytes(signal[..4].try_into().unwrap()),
                u32::from_le_bytes(signal[4..].try_into().unwrap()),
            )
        })
        .collect()
}

/// Registers the registry file `name` (kind 2) with `bytes`, for one test: names are unique per
/// test, since the registrations are the process's and the tests run at once.
struct Registered(u32);

impl Registered {
    fn new(name: &str, bytes: &[u8]) -> Self {
        let file = crate::files::begin(crate::files::KIND_REGISTRY, name, bytes.len()).unwrap();
        crate::files::chunk(file, bytes).unwrap();
        let response = crate::files::end(file).unwrap();
        assert!(response.contains("\"kind\":\"registry\""), "{response}");
        Registered(file)
    }
}

impl Drop for Registered {
    fn drop(&mut self) {
        let _ = crate::files::release(self.0);
    }
}

// NCBI BLAST+ 2.17.0 (probes `log_warn`, `log_empty`): the warning of `-max_target_seqs 1` goes to
// the `-logfile` file in TBLASTX and to the standard error in the other programs; the file is made
// in every program and the outputs are the same as without it. In LOSAT Web the file comes in the
// archive of the run's files (stream 8), `""` and `/dev/null` make none, `-` is the standard error
// (stream 3). A name with a directory is a file of the run under its whole name, as the command line
// makes it where the directory exists; so is `/dev/stdout`, whose bytes the command line writes into
// its standard output. The names of the standard error, and of a directory (which NCBI cannot open
// as a file), send the messages to stream 3, and another spelling of `/dev/null` nowhere (NCBI
// BLAST+ 2.17.0 and LOSAT's command line: `fix-logfile-dir/native` of line WQA, D11, 2026-11-20; a
// name with a directory was refused before).
#[test]
fn the_logfile_is_a_file_of_the_run() {
    for (program, _, _) in SEARCHES {
        let to_log = program == "tblastx";
        let warning = format!("Warning: [{program}] Examining 5 or more matches is recommended\n");
        let plain_words = argv(program, &["-outfmt", "6", "-max_target_seqs", "1"]);
        let (result, plain) = run_of(&plain_words);
        result.unwrap();
        assert_eq!(
            stream_text(&plain, STREAM_DIAGNOSTICS),
            warning,
            "{program}"
        );
        assert_eq!(run_files(&plain), None, "{program}");
        // The host is told where the run's messages go: stream 3 without -logfile, the file of
        // -logfile (stream 8) in TBLASTX with one (line WQA's sweep, D7, 2026-11-20).
        assert_eq!(
            messages_place(&plain),
            Some(STREAM_DIAGNOSTICS),
            "{program}"
        );
        let mut words = plain_words.clone();
        words.extend(["-logfile", "run.log"]);
        validate(&words).unwrap();
        let (result, logged) = run_of(&words);
        result.unwrap();
        for stream in [0, 6, 7, 1] {
            assert_eq!(
                logged.get(&stream),
                plain.get(&stream),
                "{program} {stream}"
            );
        }
        let (stderr, file) = match to_log {
            true => (String::new(), warning.clone()),
            false => (warning.clone(), String::new()),
        };
        assert_eq!(
            (stream_text(&logged, STREAM_DIAGNOSTICS), run_files(&logged)),
            (
                stderr.clone(),
                Some(vec![("run.log".to_string(), file.clone())])
            ),
            "{program}"
        );
        let place = match to_log {
            true => STREAM_RUN_FILES,
            false => STREAM_DIAGNOSTICS,
        };
        assert_eq!(messages_place(&logged), Some(place), "{program}");
        // `-logfile=NAME` is the same word.
        let mut joined = plain_words.clone();
        joined.push("-logfile=run.log");
        assert_eq!(run_of(&joined).1, logged, "{program}");
        for name in [
            "logs/run.log",
            "./run.log",
            "a/b/run.log",
            "/tmp/run.log",
            "/dev/stdout",
        ] {
            let mut words = plain_words.clone();
            words.extend(["-logfile", name]);
            validate(&words).unwrap();
            let (result, streams) = run_of(&words);
            result.unwrap();
            for stream in [0, 6, 7, 1] {
                assert_eq!(
                    streams.get(&stream),
                    plain.get(&stream),
                    "{program} {name} {stream}"
                );
            }
            assert_eq!(
                (
                    stream_text(&streams, STREAM_DIAGNOSTICS),
                    run_files(&streams)
                ),
                (stderr.clone(), Some(vec![(name.to_string(), file.clone())])),
                "{program} {name}"
            );
            assert_eq!(messages_place(&streams), Some(place), "{program} {name}");
        }
        for (value, quiet) in [
            ("", to_log),
            ("/dev/null", to_log),
            ("//dev/null", to_log),
            ("-", false),
            ("/dev/stderr", false),
            ("/dev/fd/2", false),
            ("/proc/self/fd/2", false),
            (".", false),
            ("..", false),
            ("logs/", false),
            ("logs/.", false),
        ] {
            let mut words = plain_words.clone();
            words.extend(["-logfile", value]);
            let (result, streams) = run_of(&words);
            result.unwrap();
            let expected = if quiet { "" } else { warning.as_str() };
            assert_eq!(
                (
                    stream_text(&streams, STREAM_DIAGNOSTICS),
                    run_files(&streams)
                ),
                (expected.to_string(), None),
                "{program} {value:?}"
            );
            let place = if quiet {
                PLACE_NOWHERE
            } else {
                STREAM_DIAGNOSTICS
            };
            assert_eq!(messages_place(&streams), Some(place), "{program} {value:?}");
            assert_eq!(streams.get(&6), plain.get(&6), "{program} {value:?}");
        }
    }
}

// NCBI BLAST+ 2.17.0 (probes `conf_ok`, `conf_eq`, `conf_bad`, `conf_missing`): the file of
// `-conffile` is read in place of `<program>.ini`; a malformed one ends the program (exit status
// 2) with NCBI's report, written to the standard error, and where the diagnostics go; a name that
// is not a file is NCBI's critical post and abort (exit status 255). In LOSAT Web the file is a
// registered registry file (kind 2), found by its name, also after `./`; the error is the text
// that the command line writes without `-logfile`, and the file of `-logfile` has its part.
#[test]
fn the_conffile_is_a_registered_registry_file() {
    let good = Registered::new("e1w3a_good.ini", b"[BLAST]\n");
    let bad = Registered::new("e1w3a_bad.ini", b"[x\n");
    for (program, _, _) in SEARCHES {
        let plain_words = argv(program, &["-outfmt", "6"]);
        let (result, plain) = run_of(&plain_words);
        result.unwrap();
        for conffile in [
            &["-conffile", "e1w3a_good.ini"][..],
            &["-conffile=e1w3a_good.ini"],
            &["-conffile", "./e1w3a_good.ini"],
        ] {
            let mut words = plain_words.clone();
            words.extend_from_slice(conffile);
            validate(&words).unwrap();
            let (result, streams) = run_of(&words);
            result.unwrap();
            assert_eq!(streams, plain, "{program} {conffile:?}");
        }
        let mut words = plain_words.clone();
        words.extend(["-conffile", "e1w3a_bad.ini"]);
        let report = SECTION_REPORT.repeat(2).trim_end().to_string();
        assert_eq!(validate(&words).unwrap_err(), report, "{program}");
        let (result, streams) = run_of(&words);
        assert_eq!(result.unwrap_err(), report, "{program}");
        assert_eq!(run_files(&streams), None, "{program}");
        // With `-logfile`: the programs whose application holds its diagnostic handler as a member
        // replay the collected messages to the standard error, the others into the file
        // (`ncbi_environment::replays_to_stderr`).
        words.extend(["-logfile", "run.log"]);
        let (result, streams) = run_of(&words);
        assert_eq!(result.unwrap_err(), report, "{program}");
        let file = match program {
            "tblastx" => SECTION_REPORT.repeat(2),
            _ => SECTION_REPORT.to_string(),
        };
        assert_eq!(
            run_files(&streams),
            Some(vec![("run.log".to_string(), file)]),
            "{program}"
        );
        let mut words = plain_words.clone();
        words.extend(["-conffile", "e1w3a_missing.ini"]);
        let critical =
            "Critical: (106.23) Registry file \"e1w3a_missing.ini\" not found, aborting.";
        assert_eq!(validate(&words).unwrap_err(), critical, "{program}");
        assert_eq!(run_of(&words).0.unwrap_err(), critical, "{program}");
    }
    drop((good, bad));
}

// NCBI BLAST+ 2.17.0 (line RG `E/SURVEY_E1b.md`): `-dryrun` ends the program after its start, its
// arguments and the settings, with exit status 0 and no output; the arguments' errors are still
// errors.
#[test]
fn a_dry_run_ends_after_the_start_and_the_arguments() {
    for (program, _, _) in SEARCHES {
        let words = argv(program, &["-outfmt", "6", "-dryrun"]);
        validate(&words).unwrap();
        let (result, streams) = run_of(&words);
        result.unwrap();
        // No output; the host is told where the run's messages go (stream 10).
        assert_eq!(
            streams.keys().copied().collect::<Vec<u32>>(),
            [STREAM_MESSAGES],
            "{program}"
        );
        assert_eq!(
            messages_place(&streams),
            Some(STREAM_DIAGNOSTICS),
            "{program}"
        );
        let (result, streams) = run_of(&argv(program, &["-dryrun", "-logfile", "run.log"]));
        result.unwrap();
        assert_eq!(
            run_files(&streams),
            Some(vec![("run.log".to_string(), String::new())]),
            "{program}"
        );
        let place = match program {
            "tblastx" => STREAM_RUN_FILES,
            _ => STREAM_DIAGNOSTICS,
        };
        assert_eq!(messages_place(&streams), Some(place), "{program}");
        assert!(validate(&argv(program, &["-dryrun", "-evalue", "x"])).is_err());
    }
}

// The *start* JSON (docs/web/abi_v2.md §9): what NCBI's application reads before the arguments.
#[test]
fn the_start_json_tells_the_files_and_the_dry_run() {
    assert_eq!(
        start_json(&["blastn", "-outfmt", "6"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"dry_run":false,"register_inputs":false}"#
    );
    assert_eq!(
        start_json(&[
            "tblastx",
            "-conffile",
            "a.ini",
            "-conffile=c.ini",
            "-logfile",
            "x.log",
            "-dryrun"
        ])
        .unwrap(),
        r#"{"conffile":"c.ini","log_file":"x.log","dry_run":true,"register_inputs":false}"#
    );
    // After `--` nothing is the start's.
    assert_eq!(
        start_json(&["blastp", "--", "-dryrun"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"dry_run":false,"register_inputs":false}"#
    );
    let version = start_json(&["blastn", "-version"]).unwrap_err();
    assert_eq!(version, parse(&["blastn", "-version"]).unwrap_err());
    assert!(version.contains("-version"), "{version}");
    assert_eq!(
        start_json(&["psiblast"]).unwrap_err(),
        "psiblast is not available in the browser (ABI v2)"
    );
}

// NCBI BLAST+ 2.17.0 (probes `log_warn`, `log_engine_err` of line RG `E/SURVEY_E1.md`): TBLASTX
// writes the messages of its run, and the error that ends it, to the `-logfile` file. In LOSAT Web
// the file of a run that failed comes on stream 8 too. The command line's file of the failed run
// (ws, the gate's build of this branch, NCBI's text): `BLAST query/options error: expect value ...`.
#[test]
fn tblastx_writes_its_run_and_its_error_to_the_logfile() {
    let failed = argv("tblastx", &["-evalue", "0", "-logfile", "run.log"]);
    let (result, streams) = run_of(&failed);
    assert!(result.is_err());
    assert_eq!(
        run_files(&streams),
        Some(vec![(
            "run.log".to_string(),
            "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n".to_string()
        )])
    );
}

// Line WQA's sweep, D7 (2026-11-20): each output that a run could not write is signalled with the
// place where its error is (stream 9), which is where the sink wrote it (`CallerSink::error`), for
// each `-logfile`. The tabular formats' hit list size of 0 fails them alone with NCBI's message,
// which goes where the run's messages go: the CLI's `tblastx -outfmt 6 -num_alignments 0` writes
// it to the file of `-logfile run.log`, nowhere for `""` and `/dev/null`, and to the standard error
// for `-` (LOSAT's command line and NCBI BLAST+ 2.17.0); BLASTN's go to the standard error. A
// subject title that NCBI cannot decode fails them alone with LOSAT's own refusal, which the CLI's
// `tblastx -outfmt 6` writes to the standard error whatever `-logfile` is (`cli::main_entry`, its
// file empty). LOSAT Web writes the run of outfmt 0 and the run of the tabular formats as one run
// (the argv's `-outfmt 0` here): its run's messages hold both searches' messages, where the CLI's
// run of the argv alone writes none of the tabular search's.
#[test]
fn the_error_of_an_output_that_fails_alone_is_where_its_signal_says() {
    const CODING: &str = "ATGAAAACCGCGTATATTGCGAAACAGCGTCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGTCAGCTGGAAGAACGTCTGGGCCTGATTGAAGTGCAGGCGCCGATTCTGAGCCGTGTGGGCGATGGCACCCAGGATAACCTGAGCGGCGCGGAAAAAGCGGTGCAGGTGAAAGTGAAAGCGCTGCCGGATGCGCAG";
    let query = format!(">q\n{CODING}\n").into_bytes();
    let bad = [
        &b">Subject_9 bad x\x81y title\n"[..],
        CODING.as_bytes(),
        b"\n",
    ]
    .concat();
    let saved = "No hits are being saved";
    let refusal = "subject record 1 (Subject_9) has a title";
    for program in ["tblastx", "blastn"] {
        let (result, plain) = run_of(&argv(program, &["-num_alignments", "0"]));
        result.unwrap();
        for logfile in [
            None,
            Some("-"),
            Some("run.log"),
            Some(""),
            Some("/dev/null"),
        ] {
            let extra: Vec<&str> = logfile.map_or(Vec::new(), |name| vec!["-logfile", name]);
            let label = format!("{program} {logfile:?}");
            // NCBI's message.
            let words = argv(program, &[&["-num_alignments", "0"][..], &extra].concat());
            let (result, streams) = run_of(&words);
            result.unwrap();
            assert_eq!(streams.get(&0), plain.get(&0), "{label}");
            assert_eq!(
                run_files(&streams).is_some(),
                logfile == Some("run.log"),
                "{label}"
            );
            let place = place_of(&streams, saved);
            assert_eq!(signals(&streams), [(6, place), (7, place)], "{label}");
            assert_eq!(
                place,
                match (program, logfile) {
                    ("tblastx", Some("run.log")) => STREAM_RUN_FILES,
                    ("tblastx", Some("" | "/dev/null")) => PLACE_NOWHERE,
                    _ => STREAM_DIAGNOSTICS,
                },
                "{label}"
            );
            // The run's messages go there too (stream 10), and stream 3 holds none of them elsewhere.
            assert_eq!(messages_place(&streams), Some(place), "{label}");
            if place != STREAM_DIAGNOSTICS {
                assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), "", "{label}");
            }
            // LOSAT's own refusal.
            let words = [program, "-query", "q", "-subject", "s"];
            let words = [&words[..], &extra].concat();
            let (result, streams) = run_on(&words, &query, &bad);
            result.unwrap();
            assert_eq!(place_of(&streams, refusal), STREAM_DIAGNOSTICS, "{label}");
            assert_eq!(
                signals(&streams),
                [(6, STREAM_DIAGNOSTICS), (7, STREAM_DIAGNOSTICS)],
                "{label}"
            );
            let on_3 = stream_text(&streams, STREAM_DIAGNOSTICS);
            assert_eq!(on_3.matches(refusal).count(), 1, "{label}");
            // Where the run's messages go elsewhere, stream 3 holds the refusal alone.
            if messages_place(&streams) != Some(STREAM_DIAGNOSTICS) {
                assert!(
                    on_3.starts_with("Error: subject record 1 (Subject_9)")
                        && on_3.matches('\n').count() == 1,
                    "{label}: {on_3}"
                );
            }
            if logfile == Some("run.log") {
                assert_eq!(
                    run_files(&streams),
                    Some(vec![("run.log".to_string(), String::new())]),
                    "{label}"
                );
            }
        }
    }
}

// NCBI reads `-outfmt` in `Run`, after a dry run ends (line RG `E/SURVEY_E1b.md`): a value that a run
// rejects passes a dry run, as on the command line.
#[test]
fn a_dry_run_does_not_read_the_outfmt() {
    for outfmt in ["x", "99", "17"] {
        let words = argv("tblastx", &["-outfmt", outfmt, "-dryrun"]);
        validate(&words).unwrap();
        let (result, streams) = run_of(&words);
        result.unwrap();
        assert_eq!(
            streams.keys().copied().collect::<Vec<u32>>(),
            [STREAM_MESSAGES],
            "{outfmt}"
        );
        let words = argv("tblastx", &["-outfmt", outfmt]);
        assert!(validate(&words).is_err(), "{outfmt}");
    }
}

// A registered registry file is the file of its name in the working directory: `./NAME` and `NAME`
// are one file, whichever the host registered or the argv wrote.
#[test]
fn a_leading_dot_slash_names_the_same_registry_file() {
    let _registered = Registered::new("./e1w3a_dot.ini", b"[x\n");
    let report = SECTION_REPORT.repeat(2).trim_end().to_string();
    for name in ["e1w3a_dot.ini", "./e1w3a_dot.ini", "././e1w3a_dot.ini"] {
        let words = argv("blastp", &["-conffile", name]);
        assert_eq!(validate(&words).unwrap_err(), report, "{name}");
    }
}

// The engine's list of the programs whose ABI v2 runs take the toolkit's words names the programs
// that the adapter runs (`Program::parse`).
#[test]
fn the_engines_list_of_web_programs_is_the_adapters() {
    for name in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpsblast",
        "rpstblastn",
        "makeblastdb",
        "rpsbproc",
    ] {
        assert_eq!(
            Program::parse(name).is_ok(),
            losat::web_adapter::WEB_V2_PROGRAMS.contains(&name),
            "{name}"
        );
    }
}
