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
    let (handle, _) = crate::store::register(words[0], 0, loaders, &bytes).unwrap();
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
        assert!(
            error.starts_with("BLAST Database error: No alias or index file found"),
            "{name}: {error}"
        );
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
