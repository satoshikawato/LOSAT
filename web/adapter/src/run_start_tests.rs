//! The NCBI application start of a run (line RG item E1w-3a): `-logfile`, `-conffile` and
//! `-dryrun` as the command line takes them (`LOSAT/tests/ncbi_toolkit_file_options.rs`, whose
//! command line arms match NCBI BLAST+ 2.17.0: probes `log_warn`, `conf_ok`, `conf_bad`,
//! `conf_missing` of line RG `E/SURVEY_E1.md`).

use super::*;
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
    let program = Program::parse(words[0]).unwrap();
    let (query, _) = crate::store::register(words[0], 0, &input(query)).unwrap();
    let (subject, _) = crate::store::register(words[0], 1, &input(subject)).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = crate::store::with_inputs(
        program,
        query,
        subject,
        writes_archive(words),
        |queries, subjects| run(words, queries, subjects),
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
// (stream 3), and a name with a directory is refused.
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
            (stderr, Some(vec![("run.log".to_string(), file)])),
            "{program}"
        );
        // `-logfile=NAME` is the same word.
        let mut joined = plain_words.clone();
        joined.push("-logfile=run.log");
        assert_eq!(run_of(&joined).1, logged, "{program}");
        for (value, quiet) in [("", to_log), ("/dev/null", to_log), ("-", false)] {
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
            assert_eq!(streams.get(&6), plain.get(&6), "{program} {value:?}");
        }
        let mut words = plain_words.clone();
        words.extend(["-logfile", "d/run.log"]);
        let refusal = "the -logfile name d/run.log has a directory, which LOSAT Web does not have: give a file name without one";
        assert_eq!(validate(&words).unwrap_err(), refusal, "{program}");
        let (result, streams) = run_of(&words);
        assert_eq!(result.unwrap_err(), refusal, "{program}");
        assert_eq!(run_files(&streams), None, "{program}");
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
        assert_eq!(
            streams.keys().copied().collect::<Vec<u32>>(),
            Vec::<u32>::new(),
            "{program}"
        );
        let (result, streams) = run_of(&argv(program, &["-dryrun", "-logfile", "run.log"]));
        result.unwrap();
        assert_eq!(
            run_files(&streams),
            Some(vec![("run.log".to_string(), String::new())]),
            "{program}"
        );
        assert!(validate(&argv(program, &["-dryrun", "-evalue", "x"])).is_err());
    }
}

// The *start* JSON (docs/web/abi_v2.md §9): what NCBI's application reads before the arguments.
#[test]
fn the_start_json_tells_the_files_and_the_dry_run() {
    assert_eq!(
        start_json(&["blastn", "-outfmt", "6"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"dry_run":false}"#
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
        r#"{"conffile":"c.ini","log_file":"x.log","dry_run":true}"#
    );
    // After `--` nothing is the start's.
    assert_eq!(
        start_json(&["blastp", "--", "-dryrun"]).unwrap(),
        r#"{"conffile":null,"log_file":null,"dry_run":false}"#
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
// the file of a run that failed comes on stream 8 too, and a run of two searches (here the tabular
// formats' hit list size of 0, which fails them alone) writes into it what stream 3 has without
// `-logfile`. The command line's file of the failed run (ws, the gate's build of this branch, NCBI's
// text): `BLAST query/options error: expect value ...`.
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
    for program in ["tblastx", "blastn"] {
        let plain_words = argv(program, &["-num_alignments", "0"]);
        let (result, plain) = run_of(&plain_words);
        result.unwrap();
        let stderr = stream_text(&plain, STREAM_DIAGNOSTICS);
        assert!(
            stderr.contains("No hits are being saved"),
            "{program}: {stderr}"
        );
        let mut words = plain_words.clone();
        words.extend(["-logfile", "run.log"]);
        let (result, logged) = run_of(&words);
        result.unwrap();
        assert_eq!(logged.get(&0), plain.get(&0), "{program}");
        let expected = match program {
            "tblastx" => (String::new(), stderr),
            _ => (stderr, String::new()),
        };
        assert_eq!(
            (stream_text(&logged, STREAM_DIAGNOSTICS), run_files(&logged)),
            (expected.0, Some(vec![("run.log".to_string(), expected.1)])),
            "{program}"
        );
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
        assert!(streams.is_empty(), "{outfmt}");
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
