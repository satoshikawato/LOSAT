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
    run_inputs_of(words, &fasta(query), &fasta(subject))
}

/// `run_files_of` with the FASTA bytes `query` and `subject`.
fn run_inputs_of(
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

/// Each row's file is a file of the run under the command line's name of it (line RG item
/// F1w-r2: `s.asn` for `./s.asn`, `outfmt11/../s.asn` and `/proc/self/cwd/s.asn`), with NCBI's
/// bytes, whatever other searches the run makes (`blastp.sub.nd`: outfmt 0's hit list of 3 and
/// outfmt 6's of 500 are two searches, and the file is the argv's, outfmt 0's), next to the file of
/// `-logfile`, and the outputs are those of the run without the option, also in the folder of the
/// query (`outfmt11/s.asn` next to `-query outfmt11/...`: a folder of the run's files). The start's
/// JSON names the file as the run does (`strategy_file`).
#[test]
fn the_strategy_file_is_a_file_of_the_run_with_the_command_lines_bytes() {
    for row in &ROWS {
        let (name, _, query, subject, _) = *row;
        for (typed, file) in [
            ("s.asn", "s.asn"),
            ("./s.asn", "s.asn"),
            ("outfmt11/s.asn", "outfmt11/s.asn"),
            ("outfmt11/./s.asn", "outfmt11/s.asn"),
            ("outfmt11/../s.asn", "s.asn"),
            ("/proc/self/cwd/s.asn", "s.asn"),
            ("/tmp/s.asn", "/tmp/s.asn"),
            ("/tmp//s.asn", "/tmp/s.asn"),
        ] {
            let words = row_argv(row, typed);
            validate(&words).unwrap();
            assert!(
                start_json(&words)
                    .unwrap()
                    .contains(&format!("\"strategy_file\":\"{file}\"")),
                "{name} {typed}"
            );
            let (result, streams) = run_files_of(&words, query, subject);
            result.unwrap();
            assert_eq!(
                files(&streams),
                [(file.to_string(), frozen(name))],
                "{name} {typed}"
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

/// The names that the command line cannot open are NCBI's error, at once (no file): also, read
/// against the run's files as `-logfile` is (line RG item F1w-r2), a name in a folder that they do
/// not have (`out/s.asn`, `a/b/../s.asn`) or in a file of them, and a folder of them (the query's
/// `outfmt11`; `src` of the run's directory, this test's working directory); `-` and the standard
/// output are the command line's refusals at NCBI's first request, and the outputs hold nothing;
/// the null device gives no file and the run goes on; the standard error, the file of another
/// option of the run and another file of the run's directory (`Cargo.toml` of this test's) are
/// LOSAT Web's refusals at the same point. A run with an empty query ends before that request:
/// NCBI's run, with no file for `-`.
#[test]
fn the_names_the_run_has_no_file_for() {
    let row = &ROWS[0];
    let (_, _, query, subject, _) = *row;
    for name in [
        "/etc/s.asn",
        ".",
        "out/",
        "",
        "out/s.asn",
        "a/b/../s.asn",
        "outfmt11",
        "outfmt11/arch_bp_q.faa/s.asn",
        "src",
    ] {
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
        ("Cargo.toml", "-export_search_strategy Cargo.toml (a file of the run's directory: LOSAT Web reads the run's files and writes none of them) is not supported in the browser; give a file name of the run's own"),
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
    // No file in any run's files: the start names none.
    assert!(start_json(&words)
        .unwrap()
        .contains("\"strategy_file\":null"));
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

/// The query and subject of each program's rows of `cases.tsv` (of `LOSAT/tests/fasta`).
fn program_inputs(program: &str) -> (&'static str, &'static str) {
    match program {
        "blastp" => ("outfmt11/arch_bp_q.faa", "outfmt11/arch_bp_s.faa"),
        "blastn" => ("outfmt11/arch_q.fna", "outfmt11/arch_s.fna"),
        "blastx" => ("outfmt11/arch_q.fna", "outfmt11/arch_bp_s.faa"),
        "tblastn" => ("outfmt11/arch_tbn_q.faa", "outfmt11/arch_tbn_s.fna"),
        "tblastx" => ("outfmt11/arch_tbx_q.fna", "outfmt11/arch_tbx_s.fna"),
        other => panic!("no rows of {other}"),
    }
}

/// NCBI's standard error (`ext` "err") or output ("out") of the row `name` of `cases.tsv`, in
/// LOSAT's words; empty where NCBI wrote none.
fn frozen_text(name: &str, ext: &str) -> String {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../LOSAT/tests/fixtures/strategy"
    );
    let bytes = std::fs::read(format!("{root}/{name}.{ext}")).unwrap_or_default();
    String::from_utf8(losat::web_adapter::losat_output_bytes(&bytes).into_owned()).unwrap()
}

/// A run that NCBI ends after its prolog and before its first request (`-query_loc` past the
/// query's end: `Empty CBlastQueryVector`), with a name that LOSAT refuses at that request: each of
/// the run's formats (the argv's, then outfmt 6, 7 and 0) keeps the bytes written to it, as the
/// command line's report stream keeps NCBI's writes (line RG T6a), so every stream is that of the
/// run without the option, the argv's own output and its signals too (review 1 B1: one hold for
/// every format put outfmt 0's prolog into the argv's tabular or XML output), for each program
/// (review 2 N-c: BLASTX's report stream and the stops of BLASTN, TBLASTN and TBLASTX had no
/// test). The argv's outfmt 0 and the run's error are
/// NCBI's (`<program>.sub.qlocbig.stdout` of `cases.tsv`, whose argv names the inputs from
/// `LOSAT/`, as the prolog shows the subject's name).
#[test]
fn a_run_that_ends_before_its_first_request_keeps_each_formats_bytes() {
    for program in ["blastp", "blastn", "blastx", "tblastn", "tblastx"] {
        let (query, subject) = program_inputs(program);
        let (query_name, subject_name) = (
            format!("tests/fasta/{query}"),
            format!("tests/fasta/{subject}"),
        );
        for outfmt in ["6", "5", "0", "7"] {
            let plain = [
                program,
                "-query",
                query_name.as_str(),
                "-subject",
                subject_name.as_str(),
                "-query_loc",
                "100000-100100",
                "-outfmt",
                outfmt,
            ];
            let (plain_result, without) = run_files_of(&plain, query, subject);
            let row = format!("{program}.sub.qlocbig.stdout");
            assert_eq!(
                plain_result.as_ref().unwrap_err(),
                &frozen_text(&row, "err"),
                "{program} outfmt {outfmt}"
            );
            if outfmt == "0" {
                assert_eq!(
                    String::from_utf8_lossy(&without[&0]),
                    frozen_text(&row, "out"),
                    "{program}: NCBI's prolog"
                );
            }
            for name in ["-", "/dev/stdout", "/dev/stderr", query_name.as_str()] {
                let mut words = plain.to_vec();
                words.extend(["-export_search_strategy", name]);
                let (result, streams) = run_files_of(&words, query, subject);
                assert_eq!(result, plain_result, "{program} outfmt {outfmt} {name}");
                assert_eq!(streams, without, "{program} outfmt {outfmt} {name}");
            }
        }
    }
}

/// A run that LOSAT refuses at NCBI's first request (`-`) writes to its outfmt 0 output what the
/// command line's report stream writes before that refusal (line RG T6a,
/// `strategy::check_report_stream`): NCBI's prolog, which NCBI writes before its request, and so
/// the start of the run without the option up to its first query's report (`prolog_of`); the
/// refusal is the run's error, for each program.
#[test]
fn a_refused_run_writes_outfmt_0s_prolog_before_the_refusal() {
    for program in ["blastp", "blastn", "blastx", "tblastn", "tblastx"] {
        let (query, subject) = program_inputs(program);
        let (query_name, subject_name) = (
            format!("tests/fasta/{query}"),
            format!("tests/fasta/{subject}"),
        );
        let plain = [
            program,
            "-query",
            query_name.as_str(),
            "-subject",
            subject_name.as_str(),
            "-outfmt",
            "0",
        ];
        let (plain_result, without) = run_files_of(&plain, query, subject);
        plain_result.unwrap();
        let mut words = plain.to_vec();
        words.extend(["-export_search_strategy", "-"]);
        let (result, streams) = run_files_of(&words, query, subject);
        assert_eq!(
            result.unwrap_err(),
            format!(
                "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's {}; give a file name",
                program.to_ascii_uppercase()
            ),
            "{program}"
        );
        assert!(files(&streams).is_empty(), "{program}");
        assert_eq!(
            String::from_utf8_lossy(&streams[&0]),
            String::from_utf8_lossy(prolog_of(&without[&0])),
            "{program}: NCBI's prolog"
        );
    }
}

/// The prolog of the outfmt 0 report `report`: its bytes before the two newlines that start the
/// first query's report (`Query=`), which NCBI writes before its first request (the command line's
/// run of the argv with `-export_search_strategy -` writes these bytes before the refusal).
fn prolog_of(report: &[u8]) -> &[u8] {
    let first_query = report
        .windows(b"\n\nQuery=".len())
        .position(|window| window == b"\n\nQuery=")
        .expect("a query's report");
    &report[..first_query]
}

/// A `-outfmt` value that does not parse, of an argv that names a search strategy file, is read
/// where the command line reads it (`losat::web_adapter::early_format_choice`, the one reading of
/// the command line and the Web): at the program's formatting handler, after the arguments' own
/// errors (`-subject_loc`, `-query_loc`), so the error of `validate` and of the run is NCBI's for
/// the argv (`<row>.err` of `cases.tsv`, in LOSAT's words) and the file is NCBI's, empty; without
/// the option the value's error comes first, as on the command line (the `.plain` rows). Line RG
/// F1w review 2, S-a: the Web read the value first with the option too.
#[test]
fn an_outfmt_that_does_not_parse_fails_where_the_command_line_fails() {
    let mut rows: Vec<(String, &str, Vec<&str>, bool)> = Vec::new();
    for program in ["blastp", "blastn", "blastx", "tblastn", "tblastx"] {
        let words = vec!["-query_loc", "abc", "-outfmt", "abc"];
        rows.push((format!("{program}.sub.qloc.fabc"), program, words, true));
        for (tag, outfmt) in [("fx", "x"), ("f99", "99")] {
            for strategy in [true, false] {
                let plain = if strategy { "" } else { ".plain" };
                let name = format!("{program}.sub.qlocrev.{tag}{plain}");
                let words = vec!["-query_loc", "10-5", "-outfmt", outfmt];
                rows.push((name, program, words, strategy));
            }
        }
    }
    let words = vec!["-subject_loc", "abc", "-outfmt", "99"];
    rows.push(("blastn.sub.sloc.f99".to_string(), "blastn", words, true));
    rows.push((
        "blastp.sub.badfmt".to_string(),
        "blastp",
        vec!["-outfmt", "99"],
        true,
    ));
    for (name, program, words, strategy) in &rows {
        let (query, subject) = program_inputs(program);
        let mut argv = vec![*program, "-query", query, "-subject", subject];
        argv.extend(words);
        if *strategy {
            argv.extend(["-export_search_strategy", "s.asn"]);
        }
        let want = frozen_text(name, "err");
        assert!(!want.is_empty(), "{name}: NCBI's error");
        assert_eq!(validate(&argv).unwrap_err(), want, "{name}: validate");
        let (result, streams) = run_files_of(&argv, query, subject);
        assert_eq!(result.unwrap_err(), want, "{name}: run");
        let file = strategy.then(|| ("s.asn".to_string(), frozen(name)));
        assert_eq!(files(&streams), Vec::from_iter(file), "{name}: files");
        for stream in [0, 6, 7, STREAM_CUSTOM] {
            assert!(
                streams.get(&stream).is_none_or(Vec::is_empty),
                "{name}: stream {stream}"
            );
        }
    }
}

/// A run of two searches (`-num_descriptions 3 -num_alignments 2`: outfmt 0's hit list of 3, the
/// tabular formats' of 500; the argv's outfmt 0 is the second search) with a name that LOSAT
/// refuses at NCBI's first request: each search has the export of the command line's run of its
/// formats (decision F1w-5; the other search's receiver keeps nothing), so the other search, which
/// runs first, stops at its own first request with its outputs holding nothing, and the run fails
/// with the refusal of the argv's own search, its outfmt 0 output holding NCBI's prolog and the
/// tabular outputs nothing, as the CLI's run of the argv writes NCBI's prolog before the refusal
/// (line RG T6a; `ROWS[1]`, whose run with a file the first test checks; without the other
/// search's export the outcome here is the same, and the next test tells them apart).
#[test]
fn a_two_search_run_refuses_a_name_as_the_command_line() {
    let (_, _, query, subject, _) = ROWS[1];
    let words = row_argv(&ROWS[1], "-");
    let groups = {
        let (program, command, requested) = parse_run(&words).unwrap();
        search_groups(&outputs(program, &requested), &command)
    };
    assert_eq!(groups.len(), 2, "two searches");
    assert_eq!(own_search(&groups), 1, "the argv's search is the second");
    let (plain_result, without) = run_files_of(&words[..words.len() - 2], query, subject);
    plain_result.unwrap();
    let (result, streams) = run_files_of(&words, query, subject);
    assert_eq!(result.unwrap_err(), BLASTP_STANDARD_OUTPUT);
    assert!(files(&streams).is_empty());
    assert_eq!(
        String::from_utf8_lossy(&streams[&0]),
        String::from_utf8_lossy(prolog_of(&without[&0])),
        "NCBI's prolog"
    );
    for stream in [6, 7] {
        assert!(streams.get(&stream).is_none_or(Vec::is_empty), "{stream}");
    }
}

/// The other search of a run of two searches reads the names of the run's files against the same
/// run's files as the argv's own search (line RG item F1w-r2, `kept_by_no_one`). For the strategy
/// file, `ROWS[1]` of the first test (`blastp.sub.nd`, `outfmt11/s.asn`) fails without it: the other
/// search read the name against no registered names, where the query's folder is missing, and gave
/// NCBI's error. psiblast's `-out_pssm` takes the same names in the other search: a name in the
/// query's folder is the run's file and the query's own file is LOSAT Web's refusal before any
/// output, as in a run of one search (these hold before the change too: the other search's PSSM
/// writes are not checked, and the refusal comes from the argv's search).
#[test]
fn the_other_search_reads_names_against_the_runs_files() {
    // NCBI writes the PSSM of the last round with `-save_pssm_after_last_round` (V-ABI row
    // psiblast.subject.it2.pssm, web/adapter/tests/v_abi.js).
    let words = [
        "psiblast",
        "-query",
        "psiblast/q1b_gi.faa",
        "-subject",
        "psiblast/q1_subj_ids.faa",
        "-num_iterations",
        "2",
        "-save_pssm_after_last_round",
        "-num_descriptions",
        "3",
        "-num_alignments",
        "2",
        "-out_pssm",
        "psiblast/p.asn",
    ];
    let groups = {
        let (program, command, requested) = parse_run(&words).unwrap();
        search_groups(&outputs(program, &requested), &command)
    };
    assert_eq!(groups.len(), 2, "two searches");
    validate(&words).unwrap();
    let (result, streams) = run_files_of(&words, words[2], words[4]);
    result.unwrap();
    let named: Vec<String> = files(&streams).into_iter().map(|(name, _)| name).collect();
    assert_eq!(named, ["psiblast/p.asn"]);
    // The query's own file is a file of the run (registered) for both searches: LOSAT Web's
    // refusal before either search writes an output.
    let mut refused = words;
    refused[13] = "psiblast/q1b_gi.faa";
    let (result, streams) = run_files_of(&refused, words[2], words[4]);
    assert_eq!(
        result.unwrap_err(),
        "an -out_pssm file 'psiblast/q1b_gi.faa' that is a file of the run's directory (the command \
         line replaces it; LOSAT Web reads the run's files and writes none of them)"
    );
    assert!(files(&streams).is_empty());
    for stream in [0, 6, 7] {
        assert!(streams.get(&stream).is_none_or(Vec::is_empty), "{stream}");
    }
}

/// LOSAT's refusal of `-export_search_strategy -` in BLASTP.
const BLASTP_STANDARD_OUTPUT: &str = "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's BLASTP; give a file name";

/// The two-search run of `a_two_search_run_refuses_a_name_as_the_command_line` on 700 queries
/// whose outfmt 7 report is longer than a chunk, which the adapter sends as the search writes it
/// (`emit::CHUNK`): the run without the option sends that much of stream 7, and the run with `-`
/// sends no byte of the other search's outputs and only NCBI's prolog of outfmt 0 (as the first
/// test), as its other search, which runs first, stops at its first request (decision F1w-5). Without that search's export it would run to its end and send its chunks
/// before the argv's own search fails (review 2 of line RG F1w, N-c, mutation M4).
#[test]
fn a_two_search_run_that_refuses_a_name_sends_nothing_of_the_other_search() {
    let title = "t".repeat(1600);
    let query: Vec<u8> = (0..700)
        .flat_map(|index| format!(">q{index} {title}\nMKTAYIAKQRQISFVKSHFS\n").into_bytes())
        .collect();
    let subject = fasta(ROWS[1].3);
    let refused = row_argv(&ROWS[1], "-");
    let (result, without) = run_inputs_of(&refused[..refused.len() - 2], &query, &subject);
    result.unwrap();
    let sent = without.get(&7).map_or(0, Vec::len);
    assert!(
        sent > crate::emit::CHUNK,
        "stream 7 of the run without the option: {sent} bytes"
    );
    let (result, streams) = run_inputs_of(&refused, &query, &subject);
    assert_eq!(result.unwrap_err(), BLASTP_STANDARD_OUTPUT);
    assert!(files(&streams).is_empty());
    assert_eq!(
        String::from_utf8_lossy(&streams[&0]),
        String::from_utf8_lossy(prolog_of(&without[&0])),
        "NCBI's prolog"
    );
    for stream in [6, 7] {
        assert!(streams.get(&stream).is_none_or(Vec::is_empty), "{stream}");
    }
}

/// The result of `run` of `words` (program name first), a search of `-db` (line WQA D12a), on the
/// query file `query` (of `LOSAT/`), and the bytes of each stream that it emitted, joined. The
/// names of the argv that start with `tests/` are files of `LOSAT/`, given from the adapter's
/// directory, where its tests run (`../../LOSAT/tests/...`): a relative name, as the command line's
/// run from `LOSAT/` gives it to the strategy's request (`db_frozen`).
fn database_run(words: &[&str], query: &str) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let words: Vec<String> = words
        .iter()
        .map(|word| match word.strip_prefix("tests/") {
            Some(name) => format!("{DB_ROOT}/tests/{name}"),
            None => word.to_string(),
        })
        .collect();
    let words: Vec<&str> = words.iter().map(String::as_str).collect();
    let bytes = std::fs::read(format!("{DB_ROOT}/{query}")).unwrap();
    let loaders = data_loaders(&words).unwrap();
    let (handle, _) = crate::store::register(words[0], 0, loaders, false, &bytes).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let program = Program::parse(words[0]).unwrap();
    let result = crate::store::with_inputs(
        program,
        handle,
        0,
        &words,
        Some(loaders),
        |queries, end, subjects| {
            assert!(subjects.is_empty());
            run(&words, queries, end, subjects)
        },
    );
    crate::store::release(handle).unwrap();
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

/// `LOSAT/` from the adapter's directory, the current directory of its tests.
const DB_ROOT: &str = "../../LOSAT";

/// NCBI's strategy file of the `-db` row `name` (`frozen`) as the run from the adapter's directory
/// writes it: the request names the database as the argv does (`subject database "..."`), here
/// `../../LOSAT/tests/...` for NCBI's `tests/...` (`database_run`); nothing else of the file names
/// a file (the lists of the filters are in it as their ids).
fn db_frozen(name: &str) -> Vec<u8> {
    let text = String::from_utf8(frozen(name)).unwrap();
    text.replace("\"tests/", &format!("\"{DB_ROOT}/tests/"))
        .into_bytes()
}

/// The rows of `cases.tsv` with `-db` that a browser's run takes as it is (no environment), one
/// of the lower case of the query's records and one of a filter of the database (a list of ids, a
/// list or a set of taxids) for each program: the row, its query (of `LOSAT/`) and its argv without
/// `-query` and the option.
const DB_ROWS: [(&str, &str, &[&str]); 10] = [
    (
        "blastp.db.lcase",
        "tests/fasta/outfmt11/arch_bp_lcase.faa",
        &[
            "blastp",
            "-db",
            "tests/db_search/arch_sprot",
            "-lcase_masking",
            "-outfmt",
            "6",
        ],
    ),
    (
        "blastp.db.list",
        "tests/db_search/taxfilter/q.faa",
        &[
            "blastp",
            "-db",
            "tests/db_search/taxfilter/tfp",
            "-seqidlist",
            "tests/db_search/seqidlist/sl_tfp.txt",
            "-outfmt",
            "6",
        ],
    ),
    (
        "blastn.db.lcase",
        "tests/fasta/outfmt11/arch_lcout.fna",
        &[
            "blastn",
            "-db",
            "tests/db_search/arch_snucl",
            "-lcase_masking",
            "-outfmt",
            "6",
        ],
    ),
    (
        "blastn.db.list",
        "tests/db_search/taxfilter/q.fna",
        &[
            "blastn",
            "-db",
            "tests/db_search/taxfilter/tfn",
            "-taxids",
            "9606",
            "-no_taxid_expansion",
            "-outfmt",
            "6",
        ],
    ),
    (
        "blastx.db.lcase",
        "tests/fasta/outfmt11/arch_bx_lc.fna",
        &[
            "blastx",
            "-db",
            "tests/db_search/arch_sprot",
            "-lcase_masking",
            "-outfmt",
            "6",
        ],
    ),
    (
        "blastx.db.list",
        "tests/db_search/taxfilter/q.fna",
        &[
            "blastx",
            "-db",
            "tests/db_search/taxfilter/tfp",
            "-negative_taxids",
            "9606",
            "-no_taxid_expansion",
            "-outfmt",
            "6",
        ],
    ),
    (
        "tblastn.db.lcase",
        "tests/fasta/outfmt11/arch_tbn_lc.faa",
        &[
            "tblastn",
            "-db",
            "tests/db_search/arch_tbn_snucl",
            "-lcase_masking",
            "-outfmt",
            "6",
        ],
    ),
    (
        "tblastn.db.list",
        "tests/db_search/taxfilter/q.faa",
        &[
            "tblastn",
            "-db",
            "tests/db_search/taxfilter/tfn",
            "-negative_seqidlist",
            "tests/db_search/seqidlist/sl_tfn.txt",
            "-outfmt",
            "6",
        ],
    ),
    (
        "tblastx.db.lcase",
        "tests/fasta/outfmt11/arch_bx_lc.fna",
        &[
            "tblastx",
            "-db",
            "tests/db_search/arch_tbx_snucl",
            "-lcase_masking",
            "-outfmt",
            "6",
        ],
    ),
    (
        "tblastx.db.list",
        "tests/db_search/taxfilter/q.fna",
        &[
            "tblastx",
            "-db",
            "tests/db_search/taxfilter/tfn",
            "-negative_taxidlist",
            "tests/db_search/taxfilter/tl_neg.txt",
            "-no_taxid_expansion",
            "-outfmt",
            "6",
        ],
    ),
];

/// A search of `-db` (line WQA D12a) writes its strategy file as the command line does, with the
/// database in its request: each row's file is a file of the run with NCBI's bytes, and the outputs
/// are those of the run without the option. A database that is not there is NCBI's error, after
/// which the file stays empty, as the command line leaves it (`*.db.nodb`).
#[test]
fn a_database_search_writes_the_command_lines_strategy_file() {
    for (name, query, words) in DB_ROWS {
        let mut argv = vec![words[0], "-query", "q"];
        argv.extend_from_slice(&words[1..]);
        let plain = argv.clone();
        argv.extend(["-export_search_strategy", "s.asn"]);
        let (result, streams) = database_run(&argv, query);
        result.unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(
            files(&streams),
            [("s.asn".to_string(), db_frozen(name))],
            "{name}"
        );
        let (_, without) = database_run(&plain, query);
        for stream in [0, 6, 7, 1, STREAM_DIAGNOSTICS] {
            assert_eq!(
                streams.get(&stream),
                without.get(&stream),
                "{name} {stream}"
            );
        }
    }
    for (name, program, query) in [
        (
            "blastp.db.nodb",
            "blastp",
            "tests/fasta/outfmt11/arch_bp_q.faa",
        ),
        (
            "blastn.db.nodb",
            "blastn",
            "tests/fasta/outfmt11/arch_q.fna",
        ),
    ] {
        let argv = [
            program,
            "-query",
            "q",
            "-db",
            "tests/db_search/none",
            "-outfmt",
            "6",
            "-export_search_strategy",
            "s.asn",
        ];
        let (result, streams) = database_run(&argv, query);
        let error = result.unwrap_err();
        // NCBI's error in LOSAT's words (the search path is the run's own).
        let start = losat::web_adapter::losat_output_bytes(
            b"BLAST Database error: No alias or index file found",
        );
        assert!(error.as_bytes().starts_with(&start), "{name}: {error}");
        assert_eq!(frozen(name), b"");
        assert_eq!(
            files(&streams),
            [("s.asn".to_string(), Vec::new())],
            "{name}"
        );
    }
}

/// The files that a search reads from the run's working directory (a file of the database of `-db`,
/// a list of its filter, a file of BLASTN's `-filtering_db`) are the files of another option of the
/// run: LOSAT Web's refusal, at NCBI's first request (the command line's strategy file would replace
/// the bytes that the search reads).
#[test]
fn the_files_of_a_database_search_are_no_strategy_file() {
    let query = "tests/db_search/taxfilter/q.faa";
    for (name, option) in [
        ("tests/db_search/taxfilter/tfp.pin", "db"),
        ("tests/db_search/seqidlist/sl_tfp.txt", "seqidlist"),
    ] {
        let argv = [
            "blastp",
            "-query",
            "q",
            "-db",
            "tests/db_search/taxfilter/tfp",
            "-seqidlist",
            "tests/db_search/seqidlist/sl_tfp.txt",
            "-outfmt",
            "6",
            "-export_search_strategy",
            name,
        ];
        let (result, streams) = database_run(&argv, query);
        assert_eq!(
            result.unwrap_err(),
            format!(
                "-export_search_strategy {DB_ROOT}/{name} (the file of -{option} too, whose bytes the command line's strategy file replaces or mixes with its own) is not supported in the browser; give a file name of the run's own"
            ),
            "{name}"
        );
        assert!(files(&streams).is_empty(), "{name}");
        assert!(streams.get(&6).is_none_or(Vec::is_empty), "{name}");
    }
    // BLASTN's `-filtering_db` (line RG item D4) reads its database of repeats from the run's
    // working directory in a search of `-subject` too.
    let name = format!("{DB_ROOT}/tests/db_search/rpt_d4.nin");
    let repeats = format!("{DB_ROOT}/tests/db_search/rpt_d4");
    let words = [
        "blastn",
        "-query",
        "outfmt11/arch_q.fna",
        "-subject",
        "outfmt11/arch_s.fna",
        "-filtering_db",
        &repeats,
        "-outfmt",
        "6",
        "-export_search_strategy",
        &name,
    ];
    let (result, streams) = run_files_of(&words, "outfmt11/arch_q.fna", "outfmt11/arch_s.fna");
    assert_eq!(
        result.unwrap_err(),
        format!(
            "-export_search_strategy {name} (the file of -filtering_db too, whose bytes the command line's strategy file replaces or mixes with its own) is not supported in the browser; give a file name of the run's own"
        )
    );
    assert!(files(&streams).is_empty());
}

/// The run's files go to the host once, after the run (`RunFileArchive::close`): a file handed over
/// after that is refused, so that a strategy file that the engine hands over late is its error
/// (LOSAT's internal error of `strategy::end_run`) and not a file lost without one (review 2 of
/// line RG F1w, N-a).
#[test]
fn a_file_handed_over_after_the_runs_files_went_out_is_refused() {
    use losat::web_adapter::FileReceiver as _;
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let archive = RunFileArchive::default();
    archive.receive(b"s.asn", b"requests").unwrap();
    RunFileArchive::emit(archive.close()).unwrap();
    let late = archive.receive(b"late.asn", b"more").unwrap_err();
    assert_eq!(
        late.to_string(),
        "late.asn came after the run's files were emitted"
    );
    assert!(archive.close().is_empty());
    let mut streams = BTreeMap::new();
    crate::emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            streams
                .entry(*stream)
                .or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        }
    });
    assert_eq!(
        files(&streams),
        [("s.asn".to_string(), b"requests".to_vec())]
    );
}
