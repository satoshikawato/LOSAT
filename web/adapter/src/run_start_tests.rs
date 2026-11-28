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

/// Two of [`SECTION_REPORT`] (NCBI's text), in LOSAT's words (`losat::web_adapter::losat_text`, AGENTS.md rule 8),
/// without the last line feed.
fn section_reports() -> String {
    losat::web_adapter::losat_text(&SECTION_REPORT.repeat(2))
        .trim_end()
        .to_string()
}

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
    let (query, _) = crate::store::register(words[0], 0, true, false, query).unwrap();
    let (subject, _) = crate::store::register(words[0], 1, true, false, subject).unwrap();
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

/// The result of `run` of a `-db` search `words` on the program's query fixture (`SEARCHES`,
/// subject handle 0), and the bytes emitted on each stream.
fn run_db_of(words: &[&str]) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let (_, query, _) = SEARCHES
        .iter()
        .find(|(program, _, _)| *program == words[0])
        .unwrap();
    let program = Program::parse(words[0]).unwrap();
    let (query, _) = crate::store::register(words[0], 0, true, false, &input(query)).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = crate::store::with_inputs(
        program,
        query,
        0,
        words,
        Some(true),
        |queries, end, subjects| run(words, queries, end, subjects),
    );
    crate::store::release(query).unwrap();
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
// its standard output, and so are `logs/` and ` run.log`, which NCBI reads as `logs` and `run.log`.
// The names of the standard error, of a directory and of what the system does not open (a file
// descriptor that the run has not opened, `/proc/self/status`, a terminal's `/dev/pts/0` in a run
// without one, a name in a device) send the messages to stream 3, as NCBI sends a log it cannot open
// to `cerr`; another spelling of `/dev/null` and `/dev/stdin` (the null device in a run without a
// standard input) send them nowhere, as a file in a file system that reports no free space does
// (`/dev/fd/2`, `/proc/self/fd/2`, `/proc/self/comm`: NCBI writes nothing to a log file with less
// than 20 KiB free). NCBI BLAST+ 2.17.0 and LOSAT's command line: `fix-logfile-dir/native` and
// `fix-logfile-dir/r2` of line WQA, D11, 2026-11-20; a name with a directory was refused before,
// and round 2 read the last ones as NCBI does. Round 3 (`fix-logfile-dir/r3`): `//run.log` has no
// directory (stream 3), a user's command line makes no file in a directory of root's (`/run.log`,
// `/etc/run.log`, `/var/log/run.log`, `/dev/run.log`) and opens none of its devices
// (`/dev/kmsg`), but a new terminal's master (`/dev/ptmx`), whose messages are lost. Round 4
// (`fix-logfile-dir/r4`): a user makes files only in the working directory, a home and the
// directories where every user does (`/var/lock` is `/run/lock`), every other directory is root's
// (`/data/run.log`), and the system refuses a component of more than 255 bytes.
/// A name of 256 bytes, one more than a component may have.
const LONG_NAME: &str = concat!(
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
);
const _: () = assert!(LONG_NAME.len() == 256);

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
        for (name, expected_name) in [
            ("./run.log", "run.log"),
            ("/tmp/run.log", "/tmp/run.log"),
            ("/dev/stdout", "/dev/stdout"),
            ("logs/", "logs"),
            ("run.log/", "run.log"),
            ("/tmp/abs/run.log/", "/tmp/abs/run.log"),
            (" run.log", "run.log"),
            ("run.log ", "run.log"),
            ("/var/lock/run.log", "/var/lock/run.log"),
            ("/var/crash/run.log", "/var/crash/run.log"),
            ("/home/me/run.log", "/home/me/run.log"),
            ("/etc/default/../../tmp/run.log", "/tmp/run.log"),
            ("/run/user/1000/run.log", "/run/user/1000/run.log"),
            ("/proc/self/cwd/run.log", "run.log"),
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
                (
                    stderr.clone(),
                    Some(vec![(expected_name.to_string(), file.clone())])
                ),
                "{program} {name}"
            );
            assert_eq!(messages_place(&streams), Some(place), "{program} {name}");
        }
        // A relative name whose folder is absent in the run's directory (or is a directory of it)
        // stays on stream 3 with no file (line WQA D19, cases 0 and 1).
        for (value, quiet) in [
            ("logs/run.log", false),
            ("a/b/run.log", false),
            ("a/b/", false),
            ("a/../run.log", false),
            ("src", false),
            ("src/", false),
            ("", to_log),
            ("/dev/null", to_log),
            ("//dev/null", to_log),
            ("/dev/fd/2", to_log),
            ("/dev/fd/1", to_log),
            ("/proc/self/fd/2", to_log),
            ("/proc/thread-self/fd/2", to_log),
            ("/dev/fd/../stderr", to_log),
            ("/dev/stdin", to_log),
            ("/proc/self/comm", to_log),
            ("-", false),
            ("/dev/stderr", false),
            ("//dev/stderr", false),
            (".", false),
            ("..", false),
            ("logs/.", false),
            ("/tmp/", false),
            ("  ", false),
            ("/dev/fd/3", false),
            ("/proc/self/status", false),
            ("/dev/pts/0", false),
            ("/dev/null/x", false),
            ("/dev/ptmx", to_log),
            ("//run.log", false),
            ("/run.log", false),
            ("/etc/run.log", false),
            ("/var/log/run.log", false),
            ("/dev/run.log", false),
            ("/dev/kmsg", false),
            ("/data/run.log", false),
            ("/proc/thread-self/../../../tmp/run.log", false),
            ("/bin/../tmp/run.log", false),
            ("/etc/default/../tmp/run.log", false),
            ("/root/../tmp/run.log", false),
            (LONG_NAME, false),
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
        // The command line opens the file of `-logfile` before it prints its help and its
        // version: a name in a folder present in the run's directory (`-query d/q.fa`) makes that
        // file empty, while one in a missing folder makes none (line WQA's D14 and D19).
        for (word, kind) in [("-h", "help"), ("-version", "version")] {
            let mut with_folder = vec![
                program,
                "-query",
                "d/q.fa",
                "-subject",
                "s.fa",
                "-num_threads",
                "1",
                "-logfile",
                "d/run.log",
                word,
            ];
            assert_eq!(prints(&with_folder), Some(kind), "{program} {word}");
            validate(&with_folder).unwrap();
            let (result, streams) = run_of(&with_folder);
            result.unwrap();
            assert!(
                streams.contains_key(&crate::emit::STREAM_PRINTED),
                "{program} {word}"
            );
            assert_eq!(
                run_files(&streams),
                Some(vec![("d/run.log".to_string(), String::new())]),
                "{program} {word}"
            );
            with_folder[2] = "q.fa";
            validate(&with_folder).unwrap();
            let (result, streams) = run_of(&with_folder);
            result.unwrap();
            assert_eq!(run_files(&streams), None, "{program} {word} missing folder");
        }
    }
}

/// NCBI's frozen file `<program>.<extension>` of `LOSAT/tests/fixtures/syntax_logfile` (the
/// standard error, `err`, and the file of `-logfile`, `log`, of NCBI BLAST+ 2.17.0's
/// `<program> <inputs> -logfile log -bogus 1`), in LOSAT's words.
fn frozen_logfile_case(program: &str, extension: &str) -> String {
    let path = format!(
        "{}/../../LOSAT/tests/fixtures/syntax_logfile/{program}.{extension}",
        env!("CARGO_MANIFEST_DIR")
    );
    losat::web_adapter::losat_text(&std::fs::read_to_string(path).unwrap()).into_owned()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1086-1112 (cited at
// `ncbi_app::Diag::report_arg_exception`)
// An argument error with `-logfile` (line RG item K4i-2b round 2, R-H1): NCBI writes the usage, its
// line that points to the error logs and the rule to the standard error, and its two posts to the
// file (NCBI's case `<program>.log.unk`). A run's error is that standard error, and the file of
// `-logfile` (stream 8) holds the posts; nothing goes to stream 3. `validate`, whose start keeps no
// file, gives the report as the command line writes it without `-logfile`, its posts included.
#[test]
fn an_argument_error_posts_to_the_runs_logfile() {
    for (program, _, _) in SEARCHES {
        let words = argv(program, &["-logfile", "run.log", "-bogus", "1"]);
        let (result, streams) = run_of(&words);
        assert_eq!(
            result.unwrap_err(),
            frozen_logfile_case(program, "err"),
            "{program}"
        );
        assert_eq!(
            run_files(&streams),
            Some(vec![(
                "run.log".to_string(),
                frozen_logfile_case(program, "log")
            )]),
            "{program}"
        );
        assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), "", "{program}");
        let rule = format!("{}\n\n", "=".repeat(72));
        let usage = frozen_logfile_case(program, "err")
            .replace(
                "Error in command-line arguments. See error logs for more details.\n",
                "",
            )
            .replace(&rule, "");
        assert_eq!(
            validate(&words).unwrap_err(),
            format!("{usage}{rule}{}", frozen_logfile_case(program, "log")),
            "{program}: validate"
        );
    }
}

/// Registers the WindowMasker counts fixture `mg1655_u9.obinary` as kind 1 under `name`.
fn register_counts(name: &str) -> u32 {
    let counts_bytes = std::fs::read(format!(
        "{}/../../LOSAT/tests/winmask/counts/mg1655_u9.obinary",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let handle = crate::files::begin(
        crate::files::KIND_WINDOW_MASKER_COUNTS,
        name,
        counts_bytes.len(),
    )
    .unwrap();
    crate::files::chunk(handle, &counts_bytes).unwrap();
    crate::files::end(handle).unwrap();
    handle
}

/// The search's file receiver sees `-query`, `-subject` (unless `-`), registered kind-2 registry
/// files read by the run's start (`.ncbirc`, `-conffile`, `<program>.ini`), and registered kind-1
/// `-window_masker_db` counts files (line WQA D19): their parent folders allow `-logfile`, the
/// folders themselves stay on stream 3, and the files themselves are refused at both `validate`
/// and `run`.
#[test]
fn the_search_receiver_reads_logfile_against_registered_files() {
    let refusal = |program: &str, name: &str| {
        format!(
            "a -logfile file '{name}' that is a file of the run's directory (LOSAT Web reads the \
             run's files and writes none of them) is not supported by LOSAT's {}",
            program.to_ascii_uppercase()
        )
    };
    // `-query in/q.fa` creates folder `in` and file `in/q.fa`; `-subject s.fa` creates `s.fa`.
    let base = [
        "tblastx",
        "-query",
        "in/q.fa",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
    ];
    let warning = "Warning: [tblastx] Examining 5 or more matches is recommended\n";
    let with_log = |name: &'static str| {
        let mut w = base.to_vec();
        w.extend(["-logfile", name]);
        w
    };
    // Case 0 (folder present via `-query in/q.fa`): makes `in/run.log`.
    validate(&with_log("in/run.log")).unwrap();
    let (res, streams) = run_of(&with_log("in/run.log"));
    res.unwrap();
    assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), "");
    assert_eq!(
        run_files(&streams),
        Some(vec![("in/run.log".to_string(), warning.to_string())])
    );
    // Case 1 (folder `in` itself): stays on stream 3, no file.
    validate(&with_log("in")).unwrap();
    let (res, streams) = run_of(&with_log("in"));
    res.unwrap();
    assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), warning, "in");
    assert_eq!(run_files(&streams), None);
    // Case 2 (`in/q.fa`, `s.fa`, `./s.fa`, `in/../s.fa`): refused at validate and run.
    for existing in ["in/q.fa", "s.fa", "./s.fa", "in/../s.fa"] {
        let words = with_log(existing);
        assert_eq!(validate(&words).unwrap_err(), refusal("tblastx", existing));
        assert_eq!(run_of(&words).0.unwrap_err(), refusal("tblastx", existing));
    }
    // `-query -` is stdin, not a file of the run's folder.
    validate(&[
        "tblastx",
        "-query",
        "-",
        "-subject",
        "s.fa",
        "-logfile",
        "-",
        "-num_threads",
        "1",
    ])
    .unwrap();
    // `-conffile sub/d19_conf.ini` (kind 2): creates folder `sub` and file `sub/d19_conf.ini`.
    let _conf = Registered::new("sub/d19_conf.ini", b"[BLAST]\n");
    let conf_base = [
        "tblastx",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-conffile",
        "sub/d19_conf.ini",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
    ];
    let with_conf_log = |name: &'static str| {
        let mut w = conf_base.to_vec();
        w.extend(["-logfile", name]);
        w
    };
    validate(&with_conf_log("sub/run.log")).unwrap();
    let (res, streams) = run_of(&with_conf_log("sub/run.log"));
    res.unwrap();
    assert_eq!(
        run_files(&streams),
        Some(vec![("sub/run.log".to_string(), warning.to_string())])
    );
    let (res, streams) = run_of(&with_conf_log("sub"));
    res.unwrap();
    assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), warning);
    assert_eq!(run_files(&streams), None);
    assert_eq!(
        validate(&with_conf_log("sub/d19_conf.ini")).unwrap_err(),
        refusal("tblastx", "sub/d19_conf.ini")
    );
    // BLASTN `-window_masker_db wm/d19.obinary` (kind 1): creates folder `wm` and file `wm/d19.obinary`.
    let wm_handle = register_counts("wm/d19.obinary");
    let wm_base = [
        "blastn",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-window_masker_db",
        "wm/d19.obinary",
        "-num_threads",
        "1",
    ];
    let with_wm_log = |name: &'static str| {
        let mut w = wm_base.to_vec();
        w.extend(["-logfile", name]);
        w
    };
    validate(&with_wm_log("wm/run.log")).unwrap();
    let (res, streams) = run_of(&with_wm_log("wm/run.log"));
    res.unwrap();
    assert_eq!(
        run_files(&streams),
        Some(vec![("wm/run.log".to_string(), String::new())])
    );
    let (res, streams) = run_of(&with_wm_log("wm"));
    res.unwrap();
    assert_eq!(run_files(&streams), None);
    assert_eq!(
        validate(&with_wm_log("wm/d19.obinary")).unwrap_err(),
        refusal("blastn", "wm/d19.obinary")
    );
    crate::files::release(wm_handle).unwrap();

    // Outside-cwd registered files (line WQA D19 round R1): every registered input whose resolved
    // path `-logfile` opens is refused at both `validate` and `run`.
    let outside_base = |q: &'static str, s: &'static str, log: &'static str| {
        vec![
            "tblastx",
            "-query",
            q,
            "-subject",
            s,
            "-max_target_seqs",
            "1",
            "-num_threads",
            "1",
            "-logfile",
            log,
        ]
    };
    let assert_refused = |label: &str, program: &str, words: &[&str], named: &str| {
        let want = refusal(program, named);
        assert_eq!(validate(words).unwrap_err(), want, "{label}: validate");
        assert_eq!(run_of(words).0.unwrap_err(), want, "{label}: run");
    };
    for (label, words, named) in [
        (
            "abs query",
            outside_base("/tmp/q.fa", "s.fa", "/tmp/q.fa"),
            "/tmp/q.fa",
        ),
        (
            "up query",
            outside_base("../q.fa", "s.fa", "../q.fa"),
            "../q.fa",
        ),
        (
            "up query normalized",
            outside_base("../q.fa", "s.fa", "./../q.fa"),
            "./../q.fa",
        ),
        (
            "proc cwd up query",
            outside_base("/proc/self/cwd/../q.fa", "s.fa", "../q.fa"),
            "../q.fa",
        ),
        (
            "abs subject",
            outside_base("q.fa", "/tmp/s.fa", "/tmp/s.fa"),
            "/tmp/s.fa",
        ),
        (
            "home query",
            outside_base("/home/me/q.fa", "s.fa", "/home/me/q.fa"),
            "/home/me/q.fa",
        ),
        (
            "tmp dotdot query",
            outside_base("/tmp/q.fa", "s.fa", "/tmp/x/../q.fa"),
            "/tmp/x/../q.fa",
        ),
        (
            "data query",
            outside_base("/data/q.fa", "s.fa", "/data/q.fa"),
            "/data/q.fa",
        ),
    ] {
        assert_refused(label, "tblastx", &words, named);
    }
    // Absolute registered `-conffile` (kind 2).
    let _abs_conf = Registered::new("/tmp/d19_r1_abs.ini", b"[BLAST]\n");
    let abs_conf_words = [
        "tblastx",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
        "-conffile",
        "/tmp/d19_r1_abs.ini",
        "-logfile",
        "/tmp/d19_r1_abs.ini",
    ];
    assert_refused(
        "abs conffile",
        "tblastx",
        &abs_conf_words,
        "/tmp/d19_r1_abs.ini",
    );
    // Absolute registered `-window_masker_db` (kind 1).
    let abs_wm_handle = register_counts("/tmp/d19_r1_abs.obinary");
    let abs_wm_words = [
        "blastn",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-window_masker_db",
        "/tmp/d19_r1_abs.obinary",
        "-num_threads",
        "1",
        "-logfile",
        "/tmp/d19_r1_abs.obinary",
    ];
    assert_refused(
        "abs wmdb",
        "blastn",
        &abs_wm_words,
        "/tmp/d19_r1_abs.obinary",
    );
    crate::files::release(abs_wm_handle).unwrap();
    // Controls that must keep their present result:
    // - `/tmp/run.log` with `-query /tmp/q.fa` and `../run.log` with `-query ../q.fa` make files;
    // - `/data/run.log` with `-query /data/q.fa` stays on stream 3 (`/data` is outside
    //   `user_makes`);
    // - `/proc/self/cwd/../q.fa` with `-query ../q.fa` stays on stream 3 because NCBI's
    //   `CDirEntry::ConvertToOSPath` lexically normalizes it to `/proc/self/q.fa` without
    //   following symlinks, which fails to open in procfs.
    for (q, log, want_file, want_diag) in [
        ("/tmp/q.fa", "/tmp/run.log", Some("/tmp/run.log"), ""),
        ("../q.fa", "../run.log", Some("../run.log"), ""),
        ("/data/q.fa", "/data/run.log", None, warning),
        ("../q.fa", "/proc/self/cwd/../q.fa", None, warning),
    ] {
        let ctrl = outside_base(q, "s.fa", log);
        validate(&ctrl).unwrap();
        let (res, streams) = run_of(&ctrl);
        res.unwrap();
        assert_eq!(
            stream_text(&streams, STREAM_DIAGNOSTICS),
            want_diag,
            "{log}"
        );
        assert_eq!(
            run_files(&streams),
            want_file.map(|name| vec![(name.to_string(), warning.to_string())]),
            "{log}"
        );
    }

    // Round R2 finding 1: an option that takes a value (`-export_search_strategy`) whose value
    // starts with `-` (`-query`, `-subject`) must not hide the real `-query`/`-subject` input.
    let strat_query_log = [
        "tblastx",
        "-export_search_strategy",
        "-query",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
        "-logfile",
        "q.fa",
    ];
    assert_refused(
        "strategy=-query -logfile q.fa",
        "tblastx",
        &strat_query_log,
        "q.fa",
    );
    let strat_subject_log = [
        "tblastx",
        "-query",
        "q.fa",
        "-export_search_strategy",
        "-subject",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
        "-logfile",
        "s.fa",
    ];
    assert_refused(
        "strategy=-subject -logfile s.fa",
        "tblastx",
        &strat_subject_log,
        "s.fa",
    );
    // Control (no `-logfile`): `-export_search_strategy -query` still runs and writes `-query`.
    let strat_query_ctrl = [
        "tblastx",
        "-export_search_strategy",
        "-query",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
    ];
    validate(&strat_query_ctrl).unwrap();
    let (res, streams) = run_of(&strat_query_ctrl);
    res.unwrap();
    let strat_files = run_files(&streams).expect("strategy file on stream 8");
    assert_eq!(strat_files.len(), 1);
    assert_eq!(strat_files[0].0, "-query");
    assert!(strat_files[0].1.contains("Blast4-request"));

    // Round R2 finding 2: `-db` list files accepted by LOSAT Web (`-gilist`, `-ipglist`,
    // `-seqidlist`, `-taxidlist` and their negatives) are created and read back by the CLI when
    // named by `-logfile` while missing, so LOSAT Web refuses them at `validate` and `run`,
    // without making missing parent folders present.
    let db_path = format!(
        "{}/../../LOSAT/tests/windowmasker/db/wm",
        env!("CARGO_MANIFEST_DIR")
    );
    for (opt, missing) in [
        ("-seqidlist", "d19-r2-miss-sl.txt"),
        ("-negative_seqidlist", "d19-r2-miss-nsl.txt"),
        ("-taxidlist", "d19-r2-miss-tl.txt"),
        ("-negative_taxidlist", "d19-r2-miss-ntl.txt"),
        ("-gilist", "d19-r2-miss-gl.txt"),
        ("-negative_gilist", "d19-r2-miss-ngl.txt"),
    ] {
        let words = [
            "blastn",
            "-query",
            "q.fa",
            "-db",
            &db_path,
            opt,
            missing,
            "-num_threads",
            "1",
            "-logfile",
            missing,
        ];
        let want = refusal("blastn", missing);
        assert_eq!(validate(&words).unwrap_err(), want, "{opt}: validate");
        assert_eq!(run_db_of(&words).0.unwrap_err(), want, "{opt}: run");
    }
    // `-ipglist` and `-negative_ipglist` are options of the protein searches.
    let protein_db = format!(
        "{}/../../LOSAT/tests/db_search/acbs",
        env!("CARGO_MANIFEST_DIR")
    );
    for (opt, missing) in [
        ("-ipglist", "d19-r2-miss-il.txt"),
        ("-negative_ipglist", "d19-r2-miss-nil.txt"),
    ] {
        let words = [
            "blastp",
            "-query",
            "q.fa",
            "-db",
            &protein_db,
            opt,
            missing,
            "-num_threads",
            "1",
            "-logfile",
            missing,
        ];
        let want = refusal("blastp", missing);
        assert_eq!(validate(&words).unwrap_err(), want, "{opt}: validate");
        assert_eq!(run_db_of(&words).0.unwrap_err(), want, "{opt}: run");
    }
    // Missing parent folder of a missing list file must NOT become present.
    let missing_list_dir = [
        "blastn",
        "-query",
        "q.fa",
        "-db",
        &db_path,
        "-seqidlist",
        "d19-r2-nosub/miss.txt",
        "-num_threads",
        "1",
        "-logfile",
        "d19-r2-nosub/run.log",
    ];
    validate(&missing_list_dir).unwrap();
    let (res, streams) = run_db_of(&missing_list_dir);
    assert!(res.is_err(), "missing seqidlist fails search");
    assert_eq!(run_files(&streams), None, "no log file in missing folder");

    // D19b (item 2): database lookup files (`.nal`/`.nin`, `.pal`/`.pin`, `index.alx`,
    // `-filtering_db`, multi-name `-db`, `taxdb.bti`) are refused at both `validate` and `run`
    // without making missing parent folders present; control `-logfile run.log` returns hits.
    let nucl_db = "../../LOSAT/tests/db_search/arch_tbx_snucl";
    let nucl_nal = format!("{nucl_db}.nal");
    let prot_db = "../../LOSAT/tests/db_search/arch_sprot";
    let prot_pal = format!("{prot_db}.pal");
    let multi_db = format!("{nucl_db} d19b-second-db");
    for (label, program, words, named, is_db) in [
        (
            "nucl db -logfile .nal",
            "blastn",
            vec![
                "blastn",
                "-task",
                "blastn",
                "-query",
                "q.fa",
                "-db",
                nucl_db,
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                &nucl_nal,
            ],
            nucl_nal.as_str(),
            true,
        ),
        (
            "prot db -logfile .pal",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                &prot_pal,
            ],
            prot_pal.as_str(),
            true,
        ),
        (
            "missing prot db -logfile .pin",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                "d19b-miss-prot",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                "d19b-miss-prot.pin",
            ],
            "d19b-miss-prot.pin",
            true,
        ),
        (
            "filtering_db -logfile .nal",
            "blastn",
            vec![
                "blastn",
                "-query",
                "q.fa",
                "-subject",
                "s.fa",
                "-filtering_db",
                nucl_db,
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                &nucl_nal,
            ],
            nucl_nal.as_str(),
            false,
        ),
        (
            "multi -db second -logfile .nal",
            "blastn",
            vec![
                "blastn",
                "-query",
                "q.fa",
                "-db",
                &multi_db,
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                "d19b-second-db.nal",
            ],
            "d19b-second-db.nal",
            true,
        ),
        (
            "db -logfile index.alx",
            "blastn",
            vec![
                "blastn",
                "-query",
                "q.fa",
                "-db",
                "d19b-miss-nucl",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                "index.alx",
            ],
            "index.alx",
            true,
        ),
        (
            "db -outfmt sscinames -logfile taxdb.bti",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "6 qseqid sseqid sscinames",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.bti",
            ],
            "taxdb.bti",
            true,
        ),
        (
            "db -outfmt 15 -logfile taxdb.bti",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "15",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.bti",
            ],
            "taxdb.bti",
            true,
        ),
        (
            "db -outfmt 16 -logfile taxdb.bti",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "16",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.bti",
            ],
            "taxdb.bti",
            true,
        ),
        (
            "db -outfmt 13 -logfile taxdb.bti",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "13",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.bti",
            ],
            "taxdb.bti",
            true,
        ),
        (
            "db -outfmt 14 -logfile taxdb.bti",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "14",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.bti",
            ],
            "taxdb.bti",
            true,
        ),
        (
            "subject -outfmt sscinames -logfile taxdb.bti",
            "blastp",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-subject",
                "s.fa",
                "-outfmt",
                "6 qseqid sseqid sscinames",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.bti",
            ],
            "taxdb.bti",
            false,
        ),
    ] {
        let want = refusal(program, named);
        assert_eq!(validate(&words).unwrap_err(), want, "{label}: validate");
        let run_err = if is_db {
            run_db_of(&words).0.unwrap_err()
        } else {
            run_of(&words).0.unwrap_err()
        };
        assert_eq!(run_err, want, "{label}: run");
    }
    // Line WQA D35 (WQA DECISIONS 62): a search of a database looks for `taxdb.bti` along its
    // search path whatever its report (the native open trace), so `-logfile taxdb.bti` is refused
    // with every format, as on the command line; a `-subject` search looks for it only where its
    // report looks every taxid up (the plural names above, the organism report), so the other
    // formats run with `-logfile taxdb.bti` and give the hits of `-logfile run.log`.
    for fmt in [
        "6",
        "0",
        "5",
        "6 qseqid sseqid staxids",
        "6 qseqid sseqid ssciname",
        "6 qseqid pident",
        "8",
        "11",
        "18",
    ] {
        let words = [
            "blastp",
            "-query",
            "q.fa",
            "-db",
            prot_db,
            "-outfmt",
            fmt,
            "-num_threads",
            "1",
            "-logfile",
            "taxdb.bti",
        ];
        let want = refusal("blastp", "taxdb.bti");
        assert_eq!(validate(&words).unwrap_err(), want, "db outfmt {fmt}: validate");
        assert_eq!(run_db_of(&words).0.unwrap_err(), want, "db outfmt {fmt}: run");
    }
    let subject_18 = [
        "blastp",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-outfmt",
        "18",
        "-num_threads",
        "1",
        "-logfile",
        "taxdb.bti",
    ];
    let want = refusal("blastp", "taxdb.bti");
    assert_eq!(validate(&subject_18).unwrap_err(), want, "subject 18: validate");
    assert_eq!(run_of(&subject_18).0.unwrap_err(), want, "subject 18: run");
    let subject_with = |fmt: &'static str, log: &'static str| {
        vec![
            "blastp",
            "-query",
            "q.fa",
            "-subject",
            "s.fa",
            "-outfmt",
            fmt,
            "-num_threads",
            "1",
            "-logfile",
            log,
        ]
    };
    for fmt in ["6", "0", "5", "6 qseqid sseqid staxids", "6 qseqid sseqid ssciname"] {
        let (ctrl, ctrl_streams) = run_of(&subject_with(fmt, "run.log"));
        ctrl.unwrap();
        let words = subject_with(fmt, "taxdb.bti");
        validate(&words).unwrap_or_else(|e| panic!("subject outfmt {fmt} validate: {e}"));
        let (res, streams) = run_of(&words);
        res.unwrap_or_else(|e| panic!("subject outfmt {fmt} run: {e}"));
        assert_eq!(
            stream_text(&streams, 6),
            stream_text(&ctrl_streams, 6),
            "subject outfmt {fmt} with -logfile taxdb.bti returns the hits of -logfile run.log"
        );
    }
    // Control: `-db` with `-logfile run.log` succeeds and returns hits.
    let db_ctrl = [
        "blastn",
        "-task",
        "blastn",
        "-query",
        "q.fa",
        "-db",
        nucl_db,
        "-outfmt",
        "6",
        "-num_threads",
        "1",
        "-logfile",
        "run.log",
    ];
    validate(&db_ctrl).unwrap();
    let (res, streams) = run_db_of(&db_ctrl);
    res.unwrap();
    assert!(
        !stream_text(&streams, 6).is_empty(),
        "db control returns hits on outfmt 6"
    );
    assert_eq!(
        run_files(&streams),
        Some(vec![("run.log".to_string(), String::new())])
    );
    // Controls that the command line does not read: `-logfile taxdb.btd` without a `taxdb.bti`
    // on the search path (the data file is looked for only beside the index found), and
    // `-logfile taxonomy4blast.sqlite3` of a search that does not widen a taxonomy filter
    // (`-no_taxid_expansion`; line WQA D35: with the widening the run looks for the file along
    // its search path even when there is none, and finds its own log there).
    let sqlite_absent = [
        "blastp",
        "-query",
        "q.fa",
        "-db",
        prot_db,
        "-taxids",
        "9606",
        "-outfmt",
        "6",
        "-num_threads",
        "1",
        "-logfile",
        "taxonomy4blast.sqlite3",
    ];
    let want = refusal("blastp", "taxonomy4blast.sqlite3");
    assert_eq!(validate(&sqlite_absent).unwrap_err(), want, "sqlite absent: validate");
    assert_eq!(run_db_of(&sqlite_absent).0.unwrap_err(), want, "sqlite absent: run");
    for (label, words) in [
        (
            "taxdb.btd without taxdb.bti",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-outfmt",
                "6 qseqid sseqid sscinames",
                "-num_threads",
                "1",
                "-logfile",
                "taxdb.btd",
            ],
        ),
        (
            "taxonomy4blast.sqlite3 with -no_taxid_expansion",
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                prot_db,
                "-taxids",
                "9606",
                "-no_taxid_expansion",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
                "-logfile",
                "taxonomy4blast.sqlite3",
            ],
        ),
    ] {
        validate(&words).unwrap();
        let (res, _) = run_db_of(&words);
        res.unwrap_or_else(|e| panic!("{label}: {e}"));
    }
    // Missing parent folder of a missing `-db` must NOT become present.
    let missing_db_dir = [
        "blastn",
        "-query",
        "q.fa",
        "-db",
        "d19b-nosub/x",
        "-num_threads",
        "1",
        "-logfile",
        "d19b-nosub/run.log",
    ];
    validate(&missing_db_dir).unwrap();
    let (res, streams) = run_db_of(&missing_db_dir);
    assert!(res.is_err(), "missing db fails search");
    assert_eq!(run_files(&streams), None, "no log file in missing folder");

    // Existing alias file (`DBLIST v1`) refuses `-logfile v1.nal`; existing `taxdb.bti` in a
    // registry `BLASTDB` folder refuses its paired `taxdb.btd`; existing `taxonomy4blast.sqlite3`
    // in a registry `BLASTDB` folder refuses `-logfile taxonomy4blast.sqlite3` when `-taxids`
    // expands descendants.
    struct CleanDir(std::path::PathBuf);
    impl Drop for CleanDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let tmp_dir = std::env::temp_dir().join(format!("losat-d19b-search-{}", std::process::id()));
    let _clean = CleanDir(tmp_dir.clone());
    if tmp_dir.exists() {
        std::fs::remove_dir_all(&tmp_dir).unwrap();
    }
    std::fs::create_dir_all(&tmp_dir).unwrap();
    std::fs::write(tmp_dir.join("x.nal"), b"TITLE t\nDBLIST v1\n").unwrap();
    std::fs::write(tmp_dir.join("taxdb.bti"), [0u8; 32]).unwrap();
    let sqlite_fixture = format!(
        "{}/../../LOSAT/tests/db_search/taxfilter/taxonomy4blast.sqlite3",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::copy(&sqlite_fixture, tmp_dir.join("taxonomy4blast.sqlite3")).unwrap();
    let alias_db = tmp_dir.join("x");
    let v1_nal = tmp_dir.join("v1.nal");
    let alias_db_str = alias_db.to_str().unwrap();
    let v1_nal_str = v1_nal.to_str().unwrap();
    let alias_words = [
        "blastn",
        "-query",
        "q.fa",
        "-db",
        alias_db_str,
        "-outfmt",
        "6",
        "-num_threads",
        "1",
        "-logfile",
        v1_nal_str,
    ];
    assert_eq!(
        validate(&alias_words).unwrap_err(),
        refusal("blastn", v1_nal_str),
        "alias sub-entry v1.nal: validate"
    );
    let blastdb_ini = format!("[BLAST]\nBLASTDB={}\n", tmp_dir.display());
    let _tax_conf = Registered::new("d19b_tax_conf.ini", blastdb_ini.as_bytes());
    let btd_path = tmp_dir.join("taxdb.btd");
    let btd_str = btd_path.to_str().unwrap();
    let btd_words = [
        "blastp",
        "-conffile",
        "d19b_tax_conf.ini",
        "-query",
        "q.fa",
        "-db",
        prot_db,
        "-outfmt",
        "6 qseqid sseqid sscinames",
        "-num_threads",
        "1",
        "-logfile",
        btd_str,
    ];
    assert_eq!(
        validate(&btd_words).unwrap_err(),
        refusal("blastp", btd_str),
        "paired taxdb.btd when taxdb.bti exists: validate"
    );
    let sqlite_shadow_words = [
        "blastp",
        "-conffile",
        "d19b_tax_conf.ini",
        "-query",
        "q.fa",
        "-db",
        prot_db,
        "-taxids",
        "9605",
        "-outfmt",
        "6",
        "-num_threads",
        "1",
        "-logfile",
        "taxonomy4blast.sqlite3",
    ];
    assert_eq!(
        validate(&sqlite_shadow_words).unwrap_err(),
        refusal("blastp", "taxonomy4blast.sqlite3"),
        "taxonomy4blast.sqlite3 shadowing BLASTDB: validate"
    );
    // When `taxdb.bti` exists in BLASTDB, every format of a search of a database refuses
    // `-logfile taxdb.bti`, which the run would find first, in place of BLASTDB's file (line WQA
    // D35).
    for fmt in ["18", "6 qseqid sseqid ssciname", "6", "6 qseqid pident"] {
        let shadow_bti = [
            "blastp",
            "-conffile",
            "d19b_tax_conf.ini",
            "-query",
            "q.fa",
            "-db",
            prot_db,
            "-outfmt",
            fmt,
            "-num_threads",
            "1",
            "-logfile",
            "taxdb.bti",
        ];
        assert_eq!(
            validate(&shadow_bti).unwrap_err(),
            refusal("blastp", "taxdb.bti"),
            "shadow taxdb.bti under BLASTDB for outfmt {fmt}: validate"
        );
        assert_eq!(
            run_db_of(&shadow_bti).0.unwrap_err(),
            refusal("blastp", "taxdb.bti"),
            "shadow taxdb.bti under BLASTDB for outfmt {fmt}: run"
        );
    }
    // A `-logfile` named like BLASTDB's own data file is that file: refused; one named like a
    // data file beside no index found is not read (above).
    let blastdb_btd = [
        "blastp",
        "-conffile",
        "d19b_tax_conf.ini",
        "-query",
        "q.fa",
        "-db",
        prot_db,
        "-outfmt",
        "6 qseqid pident",
        "-num_threads",
        "1",
        "-logfile",
        btd_str,
    ];
    assert_eq!(
        validate(&blastdb_btd).unwrap_err(),
        refusal("blastp", btd_str),
        "BLASTDB's taxdb.btd with outfmt 6 qseqid pident: validate"
    );

    // Finding 2 (GXR R1 cases B and C):
    // - Case B: alias in BLASTDB (`B/dbs/a.pal` with `DBLIST v1`, volume `B/dbs/v1.pin`),
    //   `-db a -logfile v1.pin` is refused because `Resolver::expand` checks local copy `./v1.pin`
    //   in cwd whenever `!skip_local`, regardless of whether the DBLIST entry has a directory.
    // - Case C: alias `C/a.pal` with `DBLIST v1`, volume `C/dbs/v1.pin` found via BLASTDB,
    //   `-db C/a -logfile C/dbs/v1.pal` is refused because `Resolver::expand` re-searches `v1`
    //   through `SearchPath::find` (`.pal` before `.pin`).
    std::fs::create_dir_all(tmp_dir.join("B/dbs")).unwrap();
    std::fs::create_dir_all(tmp_dir.join("C/dbs")).unwrap();
    std::fs::write(tmp_dir.join("B/dbs/a.pal"), b"TITLE a\nDBLIST v1\n").unwrap();
    std::fs::write(tmp_dir.join("B/dbs/v1.pin"), b"x").unwrap();
    std::fs::write(tmp_dir.join("C/a.pal"), b"TITLE a\nDBLIST v1\n").unwrap();
    std::fs::write(tmp_dir.join("C/dbs/v1.pin"), b"x").unwrap();
    let ini_b = format!("[BLAST]\nBLASTDB={}\n", tmp_dir.join("B/dbs").display());
    let ini_c = format!("[BLAST]\nBLASTDB={}\n", tmp_dir.join("C/dbs").display());
    let _rb = Registered::new("d19b_case_b.ini", ini_b.as_bytes());
    let _rc = Registered::new("d19b_case_c.ini", ini_c.as_bytes());
    let case_b_words = [
        "blastp",
        "-conffile",
        "d19b_case_b.ini",
        "-query",
        "q.fa",
        "-db",
        "a",
        "-outfmt",
        "6",
        "-num_threads",
        "1",
        "-logfile",
        "v1.pin",
    ];
    assert_eq!(
        validate(&case_b_words).unwrap_err(),
        refusal("blastp", "v1.pin"),
        "case B local copy v1.pin: validate"
    );
    assert_eq!(
        run_db_of(&case_b_words).0.unwrap_err(),
        refusal("blastp", "v1.pin"),
        "case B local copy v1.pin: run"
    );
    let c_db = tmp_dir.join("C/a");
    let c_log = tmp_dir.join("C/dbs/v1.pal");
    let c_db_str = c_db.to_str().unwrap();
    let c_log_str = c_log.to_str().unwrap();
    let case_c_words = [
        "blastp",
        "-conffile",
        "d19b_case_c.ini",
        "-query",
        "q.fa",
        "-db",
        c_db_str,
        "-outfmt",
        "6",
        "-num_threads",
        "1",
        "-logfile",
        c_log_str,
    ];
    assert_eq!(
        validate(&case_c_words).unwrap_err(),
        refusal("blastp", c_log_str),
        "case C search-path re-lookup C/dbs/v1.pal: validate"
    );
    assert_eq!(
        run_db_of(&case_c_words).0.unwrap_err(),
        refusal("blastp", c_log_str),
        "case C search-path re-lookup C/dbs/v1.pal: run"
    );
}

// Line WQA D26 item 2: LOSAT Web refuses the files of a `-db` that the command line refuses (one
// producer, `SeqDb::read_paths`): besides the lookup candidates, a missing component file of a
// volume and the GI-mask files of the alias file's `MASKLIST`.
#[test]
fn the_search_receiver_refuses_the_databases_files_as_the_command_line() {
    struct CleanDir(std::path::PathBuf);
    impl Drop for CleanDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let tmp_dir = std::env::temp_dir().join(format!("losat-d26-search-{}", std::process::id()));
    let _clean = CleanDir(tmp_dir.clone());
    if tmp_dir.exists() {
        std::fs::remove_dir_all(&tmp_dir).unwrap();
    }
    std::fs::create_dir_all(&tmp_dir).unwrap();
    let fixture = format!(
        "{}/../../LOSAT/tests/db_search/arch_tbx_snucl",
        env!("CARGO_MANIFEST_DIR")
    );
    for ext in ["nhr", "nin", "nsq"] {
        std::fs::copy(format!("{fixture}.{ext}"), tmp_dir.join(format!("v.{ext}"))).unwrap();
    }
    std::fs::write(
        tmp_dir.join("m.nal"),
        format!(
            "TITLE m\nDBLIST v\nMASKLIST {}\n",
            tmp_dir.join("mk").display()
        ),
    )
    .unwrap();
    let db = tmp_dir.join("m");
    let db = db.to_str().unwrap();
    let words = |log: &str| -> Vec<String> {
        [
            "blastn",
            "-task",
            "blastn",
            "-query",
            "q.fa",
            "-db",
            db,
            "-outfmt",
            "6",
            "-num_threads",
            "1",
            "-logfile",
            log,
        ]
        .map(str::to_string)
        .to_vec()
    };
    let want = |path: &str| {
        format!(
            "a -logfile file '{path}' that is a file of the run's directory (LOSAT Web reads the \
             run's files and writes none of them) is not supported by LOSAT's BLASTN"
        )
    };
    // The files that a search of the sequences opens, found or not, are refused as the command
    // line refuses them (line WQA D26 round 1: only those the run can open).
    for name in ["v.ndb", "v.nsd", "v.nsi"] {
        let path = tmp_dir.join(name);
        let path = path.to_str().unwrap();
        let owned = words(path);
        let words: Vec<&str> = owned.iter().map(String::as_str).collect();
        assert_eq!(
            validate(&words).unwrap_err(),
            want(path),
            "{name}: validate"
        );
        assert_eq!(run_db_of(&words).0.unwrap_err(), want(path), "{name}: run");
    }
    // The GI index and the GI-mask files are opened with a mask option only.
    for name in ["v.nog", "mk.gni", "mk.gno"] {
        let path = tmp_dir.join(name);
        let path = path.to_str().unwrap();
        let mut owned = words(path);
        owned.extend(["-db_soft_mask".to_string(), "0".to_string()]);
        let masked: Vec<&str> = owned.iter().map(String::as_str).collect();
        assert_eq!(
            validate(&masked).unwrap_err(),
            want(path),
            "{name}: validate"
        );
        assert_eq!(run_db_of(&masked).0.unwrap_err(), want(path), "{name}: run");
        let owned = words(path);
        let plain: Vec<&str> = owned.iter().map(String::as_str).collect();
        validate(&plain).unwrap_or_else(|error| panic!("{name} without a mask option: {error}"));
        run_db_of(&plain)
            .0
            .unwrap_or_else(|error| panic!("{name} without a mask option: {error}"));
    }
    // The id-lookup, header-hash and metadata files are no file of this run.
    for name in ["v.nhd", "v.nnd", "v.nos", "v.njs"] {
        let path = tmp_dir.join(name);
        let path = path.to_str().unwrap();
        let owned = words(path);
        let plain: Vec<&str> = owned.iter().map(String::as_str).collect();
        validate(&plain).unwrap_or_else(|error| panic!("{name}: {error}"));
    }
    // Control: a name that the run does not read is the run's log, with the search's hits.
    let owned = words("run.log");
    let control: Vec<&str> = owned.iter().map(String::as_str).collect();
    validate(&control).unwrap();
    let (res, streams) = run_db_of(&control);
    res.unwrap();
    assert!(!stream_text(&streams, 6).is_empty(), "control returns hits");
}

/// `registered_run_files` uses one reader regardless of whether the argv parses (line WQA D19
/// round R1): the same argv with and without an unknown option and with `-h` yields the same
/// registered names for blastn (with `-window_masker_taxid`), tblastx, and psiblast.
#[test]
fn registered_run_files_uses_one_reader_regardless_of_parse() {
    let taxid_path = format!("{}/9606/wmasker.obinary", crate::files::WINDOW_MASKER_PATH);
    let handle = register_counts(&taxid_path);
    let registry = crate::files::registry_files();

    let blastn_base = [
        "blastn",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-window_masker_taxid",
        "9606",
    ];
    let blastn_bad = [
        "blastn",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-window_masker_taxid",
        "9606",
        "-no_such_option",
    ];
    let blastn_help = [
        "blastn",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-window_masker_taxid",
        "9606",
        "-h",
    ];
    let want_blastn = vec![
        b"q.fa".to_vec(),
        b"s.fa".to_vec(),
        taxid_path.as_bytes().to_vec(),
    ];
    assert_eq!(registered_run_files(&blastn_base, &registry), want_blastn);
    assert_eq!(registered_run_files(&blastn_bad, &registry), want_blastn);
    assert_eq!(registered_run_files(&blastn_help, &registry), want_blastn);

    let tblastx_base = ["tblastx", "-query", "q.fa", "-subject", "s.fa"];
    let tblastx_bad = [
        "tblastx",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-no_such_option",
    ];
    let tblastx_help = ["tblastx", "-query", "q.fa", "-subject", "s.fa", "-h"];
    let want_tblastx = vec![b"q.fa".to_vec(), b"s.fa".to_vec()];
    assert_eq!(registered_run_files(&tblastx_base, &registry), want_tblastx);
    assert_eq!(registered_run_files(&tblastx_bad, &registry), want_tblastx);
    assert_eq!(registered_run_files(&tblastx_help, &registry), want_tblastx);

    // Round R2 finding 1: value of `-export_search_strategy` named `-query` is skipped as its
    // value, yielding `q.fa` and `s.fa` with and without an unknown option and with `-h`.
    let strat_base = [
        "tblastx",
        "-export_search_strategy",
        "-query",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
    ];
    let strat_bad = [
        "tblastx",
        "-export_search_strategy",
        "-query",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-no_such_option",
    ];
    let strat_help = [
        "tblastx",
        "-export_search_strategy",
        "-query",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-h",
    ];
    assert_eq!(registered_run_files(&strat_base, &registry), want_tblastx);
    assert_eq!(registered_run_files(&strat_bad, &registry), want_tblastx);
    assert_eq!(registered_run_files(&strat_help, &registry), want_tblastx);

    let psiblast_base = ["psiblast", "-query", "q.fa", "-subject", "s.fa"];
    let psiblast_bad = [
        "psiblast",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-no_such_option",
    ];
    let psiblast_help = ["psiblast", "-query", "q.fa", "-subject", "s.fa", "-h"];
    let want_psiblast: Vec<Vec<u8>> = vec![b"q.fa".to_vec(), b"s.fa".to_vec()];
    assert_eq!(
        registered_run_files(&psiblast_base, &registry),
        want_psiblast
    );
    assert_eq!(
        registered_run_files(&psiblast_bad, &registry),
        want_psiblast
    );
    assert_eq!(
        registered_run_files(&psiblast_help, &registry),
        want_psiblast
    );
    let rps_base = ["rpsblast", "-query", "q.fa", "-db", "cdd"];
    assert_eq!(
        registered_run_files(&rps_base, &registry),
        Vec::<Vec<u8>>::new()
    );

    crate::files::release(handle).unwrap();
}

/// `-logfile out/x.asn -export_search_strategy out/x.asn` without folder `out` (line WQA D19 round
/// R1): `-logfile` stays on stderr (`Unopened`), so `-export_search_strategy out/x.asn` is not
/// refused as the file of `-logfile` and writes the strategy file on stream 8; with folder `out`
/// present, `-logfile` opens `out/x.asn` and the strategy is refused.
#[test]
fn strategy_target_sees_runs_files_for_logfile() {
    let words = [
        "tblastx",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
        "-logfile",
        "out/x.asn",
        "-export_search_strategy",
        "out/x.asn",
    ];
    validate(&words).unwrap();
    let (res, streams) = run_of(&words);
    res.unwrap();
    let diag = stream_text(&streams, STREAM_DIAGNOSTICS);
    assert_eq!(
        diag,
        "Warning: [tblastx] Examining 5 or more matches is recommended\n"
    );
    let files = run_files(&streams).expect("strategy file on stream 8");
    assert_eq!(files.len(), 1);
    assert_eq!(files[0].0, "out/x.asn");
    assert!(files[0].1.contains("Blast4-request"), "{}", files[0].1);

    // With folder `out` present via a registered `-conffile out/d19_r1_strat.ini`, `-logfile
    // out/x.asn` opens a file, so `-export_search_strategy out/x.asn` is refused.
    let _conf = Registered::new("out/d19_r1_strat.ini", b"[BLAST]\n");
    let with_folder = [
        "tblastx",
        "-query",
        "q.fa",
        "-subject",
        "s.fa",
        "-max_target_seqs",
        "1",
        "-num_threads",
        "1",
        "-conffile",
        "out/d19_r1_strat.ini",
        "-logfile",
        "out/x.asn",
        "-export_search_strategy",
        "out/x.asn",
    ];
    validate(&with_folder).unwrap();
    let (res, _streams) = run_of(&with_folder);
    assert_eq!(
        res.unwrap_err(),
        "-export_search_strategy out/x.asn (the file of -logfile too, whose bytes the command \
         line's strategy file replaces or mixes with its own) is not supported in the browser; \
         give a file name of the run's own"
    );
}

/// NCBI's error for an argv without `-subject` and `-db`, which the command line raises after
/// `-outfmt`'s number and before its other checks (the ranges, the options, the formats that the
/// program does not write); a dry run ends before it. `blastp -query -h` is such an argv: `-h` is
/// the query's name (line WQA's D14, 2026-11-21).
#[test]
fn an_argv_without_subjects_fails_as_on_the_command_line() {
    let error = losat::web_adapter::losat_text(
        "BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified\nPlease refer to the LOSAT user manual.\n",
    );
    for (program, query, _) in SEARCHES {
        for words in [
            vec![program, "-query", "-h"],
            vec![program, "-query", query, "-evalue", "0"],
            vec![program, "-query", query, "-subject_loc", "x"],
            vec![program, "-query", query, "-outfmt", "17"],
        ] {
            assert_eq!(prints(&words), None, "{words:?}");
            assert_eq!(validate(&words).unwrap_err(), *error, "{words:?}");
            assert_eq!(run_of(&words).0.unwrap_err(), *error, "{words:?}");
        }
        assert_eq!(
            validate(&[program, "-query", query, "-outfmt", "x"]).unwrap_err(),
            losat::web_adapter::losat_text(
                "BLAST query/options error: 'x' is not a valid output format\nPlease refer to the LOSAT user manual.\n"
            )
        );
        validate(&[program, "-query", query, "-dryrun"]).unwrap();
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
        let report = section_reports();
        assert_eq!(validate(&words).unwrap_err(), report, "{program}");
        let (result, streams) = run_of(&words);
        assert_eq!(result.unwrap_err(), report, "{program}");
        assert_eq!(run_files(&streams), None, "{program}");
        // The report is the error: nothing came before it on the standard error (stream 3).
        assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), "", "{program}");
        // With `-logfile`: the programs whose application holds its diagnostic handler as a member
        // replay the collected messages to the standard error, the others into the file
        // (`ncbi_environment::replays_to_stderr`).
        words.extend(["-logfile", "run.log"]);
        let (result, streams) = run_of(&words);
        assert_eq!(result.unwrap_err(), report, "{program}");
        assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), "", "{program}");
        // NCBI's text in LOSAT's words (`losat::web_adapter::losat_text`, AGENTS.md rule 8).
        let file = losat::web_adapter::losat_text(&match program {
            "tblastx" => SECTION_REPORT.repeat(2),
            _ => SECTION_REPORT.to_string(),
        })
        .into_owned();
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
        r#"{"conffile":null,"log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
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
        r#"{"conffile":"c.ini","log_file":"x.log","log_file_name":"x.log","strategy_file":null,"dry_run":true,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    // `log_file` keeps NCBI's scan of the argv; `log_file_name` is the normalized host file name
    // under which a run gives the file on stream 8 (line WQA D19, case 4).
    assert_eq!(
        start_json(&["tblastx", "-logfile", "logs/"]).unwrap(),
        r#"{"conffile":null,"log_file":"logs/","log_file_name":"logs","strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    assert_eq!(
        start_json(&["tblastx", "-logfile", "./run.log"]).unwrap(),
        r#"{"conffile":null,"log_file":"./run.log","log_file_name":"run.log","strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    assert_eq!(
        start_json(&["tblastx", "-logfile", "-"]).unwrap(),
        r#"{"conffile":null,"log_file":"-","log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    assert_eq!(
        start_json(&["tblastx", "-logfile", "/dev/stdout"]).unwrap(),
        r#"{"conffile":null,"log_file":"/dev/stdout","log_file_name":"/dev/stdout","strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    // The name of the search strategy file, under which a run gives it (line RG item F1w).
    assert_eq!(
        start_json(&[
            "blastp",
            "-export_search_strategy",
            "out/s.asn",
            "-logfile",
            "run.log"
        ])
        .unwrap(),
        r#"{"conffile":null,"log_file":"run.log","log_file_name":"run.log","strategy_file":"out/s.asn","dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    // After `--` nothing is the start's.
    assert_eq!(
        start_json(&["blastp", "--", "-dryrun"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    // A run that prints the command line's text (line WQA's D14, 2026-11-20): a version word, read
    // before the arguments and the configuration (no `-conffile` to register, not a dry run), or a
    // help word that stops the arguments' reader.
    assert_eq!(
        start_json(&[
            "blastn",
            "-dryrun",
            "-conffile",
            "c.ini",
            "-version",
            "-logfile",
            "v.log"
        ])
        .unwrap(),
        r#"{"conffile":null,"log_file":"v.log","log_file_name":"v.log","strategy_file":null,"dry_run":false,"register_inputs":false,"prints":"version","database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    assert_eq!(
        start_json(&["tblastn", "-conffile", "c.ini", "-dryrun", "-h"]).unwrap(),
        r#"{"conffile":"c.ini","log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":"help","database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    assert_eq!(
        start_json(&["blastx", "-parse_deflines", "-xmlhelp"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":"help","database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    // The words before a help word are read first: an error there is the run's.
    assert_eq!(
        start_json(&["blastp", "-evalue", "x", "-help"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":null,"out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    // PSW1: psiblast PSSM/MSA query and PSSM output files, and tblastn -in_pssm.
    assert_eq!(
        start_json(&[
            "psiblast",
            "-in_pssm",
            "seed.pssm",
            "-subject",
            "s.fa",
            "-out_pssm",
            "./out.asn",
            "-out_ascii_pssm",
            "out.ascii",
            "-save_each_pssm",
        ])
        .unwrap(),
        r#"{"conffile":null,"log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":false,"query_file":"seed.pssm","out_pssm":"out.asn","out_ascii_pssm":"out.ascii","save_each_pssm":true}"#
    );
    assert_eq!(
        start_json(&["tblastn", "-in_pssm", "seed.pssm", "-db", "nt"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"log_file_name":null,"strategy_file":null,"dry_run":false,"register_inputs":false,"prints":null,"database":true,"query_file":"seed.pssm","out_pssm":null,"out_ascii_pssm":null,"save_each_pssm":false}"#
    );
    assert_eq!(
        start_json(&["rpsblast"]).unwrap_err(),
        "rpsblast is not available in the browser (ABI v2)"
    );
}

// NCBI BLAST+ 2.17.0 (probes `log_warn`, `log_engine_err` of line RG `E/SURVEY_E1.md`): TBLASTX
// writes the messages of its run, and the error that ends it, to the `-logfile` file. In LOSAT Web
// the file of a run that failed comes on stream 8 too. The file holds what the CLI's file holds for
// the argv: in a run of two searches what the other search adds goes to stream 3, as without
// `-logfile` (`the_error_of_an_output_that_fails_alone_is_where_its_signal_says`; line WQA's sweep,
// D3 review r2, 2026-11-19: the file had the other search's warning and error). The command line's
// file of the failed run (a Linux host, the gate's build of this branch, NCBI's text):
// `BLAST query/options error: expect value ...`.
#[test]
fn tblastx_writes_its_run_and_its_error_to_the_logfile() {
    let failed = argv("tblastx", &["-evalue", "0", "-logfile", "run.log"]);
    let (result, streams) = run_of(&failed);
    assert!(result.is_err());
    // The file has the error, and the standard error (stream 3) nothing.
    assert_eq!(stream_text(&streams, STREAM_DIAGNOSTICS), "");
    assert_eq!(
        run_files(&streams),
        Some(vec![(
            "run.log".to_string(),
            losat::web_adapter::losat_text(
                "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n"
            )
            .into_owned()
        )])
    );
}

// Line WQA's sweep, D7 (2026-11-20): each output that a run could not write is signalled with the
// place where its error is (stream 9), which is where the sink wrote it (`CallerSink::error`), for
// each `-logfile`. The tabular formats' hit list size of 0 fails them alone with NCBI's message.
// The run of the argv's `-outfmt 0` (here) writes the run's messages where the CLI's run of the
// argv writes them: in TBLASTX the file of `-logfile run.log`, which then holds what the CLI's file
// holds, nothing (the CLI's `tblastx -num_alignments 0 -logfile run.log`); what the tabular formats
// add, their error among it, goes to stream 3, as without `-logfile` (line WQA's sweep, D3 review
// r2, 2026-11-19), and nowhere for `""` and `/dev/null`, where the CLI's `tblastx -outfmt 6
// -num_alignments 0` writes it nowhere (LOSAT's command line and NCBI BLAST+ 2.17.0; D3 review r3);
// BLASTN's go to the standard error. A subject title that NCBI cannot decode fails them alone with
// LOSAT's own refusal, which the CLI's `tblastx -outfmt 6` writes to the standard error whatever
// `-logfile` is (`cli::main_entry`, its file empty).
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
        let stderr = stream_text(&plain, STREAM_DIAGNOSTICS);
        assert!(stderr.contains(saved), "{program}: {stderr}");
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
            let messages = match (program, logfile) {
                ("tblastx", Some("run.log")) => STREAM_RUN_FILES,
                ("tblastx", Some("" | "/dev/null")) => PLACE_NOWHERE,
                _ => STREAM_DIAGNOSTICS,
            };
            assert_eq!(messages_place(&streams), Some(messages), "{label}");
            // The error is where what the tabular formats add goes: stream 3, as without -logfile,
            // also where the run's messages go to the file, which holds the CLI's file of the argv
            // (empty), and nowhere where they go nowhere (line WQA's sweep, D3 review r2 and r3).
            let place = place_of(&streams, saved);
            assert_eq!(signals(&streams), [(6, place), (7, place)], "{label}");
            let nowhere = messages == PLACE_NOWHERE;
            assert_eq!(
                place,
                if nowhere {
                    PLACE_NOWHERE
                } else {
                    STREAM_DIAGNOSTICS
                },
                "{label}"
            );
            assert_eq!(
                stream_text(&streams, STREAM_DIAGNOSTICS),
                if nowhere {
                    String::new()
                } else {
                    stderr.clone()
                },
                "{label}"
            );
            if logfile == Some("run.log") {
                assert_eq!(
                    run_files(&streams),
                    Some(vec![("run.log".to_string(), String::new())]),
                    "{label}"
                );
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
    let report = section_reports();
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

/// The bytes that `run` of `words` emitted on each stream, with no input: a run that prints reads
/// none (`losat_web2_run` passes it none).
fn printing_run(words: &[&str]) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = run(words, &[], &losat::api::QueryEnd::Input, &[]);
    let mut streams = BTreeMap::new();
    crate::emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            streams
                .entry(*stream)
                .or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        }
    });
    (result, streams)
}

// Line WQA's D14 (2026-11-20): the help and version words print what the command line prints to
// its standard output (NCBI's help, the XML of `-xmlhelp`, NCBI's version texts with LOSAT's
// values: `ncbi_app::read` and `ncbi_version_request`, the command line's own producers), on stream
// 11, with no input and no search; `validate` accepts them. `web/adapter/tests/v_abi.js` compares
// the bytes with the native command line.
#[test]
fn help_and_version_words_print_the_command_lines_text() {
    for (program, _, _) in SEARCHES {
        for word in [
            "-h",
            "-help",
            "-help-full",
            "-xmlhelp",
            "-version",
            "-version-full",
            "-version-full-xml",
            "-version-full-json",
        ] {
            let words = argv(program, &["-evalue", "1", word, "-bogus"]);
            validate(&words).unwrap();
            let (result, streams) = printing_run(&words);
            result.unwrap();
            assert_eq!(
                streams.keys().copied().collect::<Vec<u32>>(),
                [STREAM_MESSAGES, crate::emit::STREAM_PRINTED],
                "{program} {word}"
            );
            let text = stream_text(&streams, crate::emit::STREAM_PRINTED);
            let expected = match word {
                // NCBI's forms with LOSAT's values (`ncbi_version::print`; v_abi.js compares the
                // bytes with the command line's).
                "-version" | "-version-full" => {
                    let version = losat::web_adapter::LOSAT_VERSION;
                    assert!(text.starts_with(&format!("{program}: {version}\n Package: LOSAT ")));
                    text.clone()
                }
                "-version-full-xml" => assert_xml(&text),
                "-version-full-json" => {
                    assert!(
                        text.contains(&format!("\"appname\": \"{program}\"")),
                        "{text}"
                    );
                    text.clone()
                }
                "-xmlhelp" => assert_xml(&text),
                _ => {
                    let error = losat::web_adapter::try_parse_cli(
                        std::iter::once("losat").chain(words.iter().copied()),
                        None,
                    )
                    .unwrap_err();
                    match error {
                        losat::web_adapter::ArgvError::Printed(text) => text,
                        error => panic!("{program} {word}: {error:?}"),
                    }
                }
            };
            assert_eq!(text, expected, "{program} {word}");
        }
        // `--help` is no word of NCBI's reader (`Unknown argument: "-help"`, line RG DECISIONS
        // K4i-2 D3), and the words before a help word are read first, as the command line reads
        // them.
        let unknown = argv(program, &["--help"]);
        assert!(validate(&unknown).is_err(), "{program}");
        let refused = argv(program, &["-evalue", "x", "-h"]);
        assert!(validate(&refused).is_err(), "{program}");
        assert!(printing_run(&refused).0.is_err(), "{program}");
    }
}

/// `text` as an XML document (its first line), returned as it is.
fn assert_xml(text: &str) -> String {
    assert!(text.starts_with("<?xml version=\"1.0\"?>\n"), "{text}");
    text.to_string()
}

// NCBI BLAST+ 2.17.0 and the command line (line WQA's D14 probes, 2026-11-20): a version is printed
// after the `-logfile` file is opened and before the configuration is loaded, so a `-conffile` of
// no file does not stop it; the help comes after the load, which the missing file ends (exit
// status 255 on the command line).
#[test]
fn a_version_comes_before_the_configuration_and_the_help_after_it() {
    let version = argv("blastn", &["-version", "-conffile", "d14_missing.ini"]);
    validate(&version).unwrap();
    let (result, streams) = printing_run(&version);
    result.unwrap();
    let losat_version = losat::web_adapter::LOSAT_VERSION;
    assert!(stream_text(&streams, crate::emit::STREAM_PRINTED)
        .starts_with(&format!("blastn: {losat_version}\n")));
    let help = argv("blastn", &["-h", "-conffile", "d14_missing.ini"]);
    let critical = "Critical: (106.23) Registry file \"d14_missing.ini\" not found, aborting.";
    assert_eq!(validate(&help).unwrap_err(), critical);
    let (result, streams) = printing_run(&help);
    assert_eq!(result.unwrap_err(), critical);
    assert!(!streams.contains_key(&crate::emit::STREAM_PRINTED));
    // The file of `-logfile` is made, empty, as the command line makes it.
    let (result, streams) = printing_run(&argv("tblastx", &["-version", "-logfile", "v.log"]));
    result.unwrap();
    assert_eq!(
        run_files(&streams),
        Some(vec![("v.log".to_string(), String::new())])
    );
}

#[test]
fn range_recorder_faults_when_hit_list_count_mismatches_parts() {
    let mut recorder = RangeRecorder::new(vec![6], vec![Arc::new(AtomicU64::new(0))]);
    recorder.round_begin(0, 1);
    recorder.round_begin(0, 2);
    assert_eq!(recorder.parts.len(), 2);

    // Zero hit lists (e.g. search that handed over no hit list) or exact match is valid.
    recorder.check_hit_lists(0);
    assert!(recorder.fault.is_none());
    recorder.check_hit_lists(2);
    assert!(recorder.fault.is_none());

    // Fewer hit lists than parts faults instead of silently truncating later rounds.
    recorder.check_hit_lists(1);
    assert_eq!(
        recorder.fault.as_deref(),
        Some("1 hit lists for 2 parts of the search")
    );

    // More hit lists than parts also faults.
    let mut extra = RangeRecorder::new(vec![6], vec![Arc::new(AtomicU64::new(0))]);
    extra.round_begin(0, 1);
    extra.check_hit_lists(3);
    assert_eq!(
        extra.fault.as_deref(),
        Some("3 hit lists for 1 parts of the search")
    );
}
