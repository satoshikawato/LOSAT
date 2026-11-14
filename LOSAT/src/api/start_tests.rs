//! The NCBI application start of the library's search calls (`api::start`, line RG item E1w-2):
//! `-dryrun`, `-logfile` and `-conffile` in `run`, `run_collect`, `run_to_file` and
//! `run_batches`, against the command line's behaviour (NCBI BLAST+ 2.17.0, probes of line RG
//! items E1 and E1b, `tests/ncbi_toolkit_file_options.rs` and `tests/ncbi_toolkit_dryrun.rs`,
//! which compare the library's bytes with the command line's through
//! `examples/library_command_line`).
//!
//! These tests change neither the environment nor the current directory of the test process
//! (in-crate tests share it): the start reads the process's environment, which has none of NCBI's
//! settings here, as the other in-crate runs of a start assume (`dustmasker`'s tests). Every file
//! is in a directory of its own test.

use std::path::{Path, PathBuf};

use super::*;
use crate::api::{
    run, run_batches, run_collect, run_to_file, validate, FastaRecord, Inputs, Program,
};

const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

const NUCLEOTIDES: &[u8] =
    b"ATGAAAACCGCGTATATTGCGAAACAGCGCCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGCCAGCTGGAAGAACGCCTGGGC";

/// NCBI's report of a registry file with the line `[x`, which ends the program (exit status 2),
/// and its replay (`tests/ncbi_toolkit_file_options.rs`).
const SECTION_REPORT: &str = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";

/// The warning of `-max_target_seqs 1` (NCBI BLAST+ 2.17.0, probe `log_warn`).
fn few_matches(program: &str) -> String {
    format!("Warning: [{program}] Examining 5 or more matches is recommended\n")
}

fn proteins() -> (Vec<FastaRecord>, Vec<FastaRecord>) {
    (
        vec![
            FastaRecord::new("Query_1", b"q1", PROTEIN),
            FastaRecord::new("Query_2", b"q2", &PROTEIN[5..]),
        ],
        vec![FastaRecord::new("Subject_1", b"s", PROTEIN)],
    )
}

fn nucleotides() -> (Vec<FastaRecord>, Vec<FastaRecord>) {
    (
        vec![FastaRecord::new("Query_1", b"q", NUCLEOTIDES)],
        vec![FastaRecord::new("Subject_1", b"s", NUCLEOTIDES)],
    )
}

/// A directory of the test `name`, empty.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-api-start-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn path_text(path: &Path) -> &str {
    path.to_str().unwrap()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (quoted at
// `ncbi_environment::preparse`): `-dryrun` anywhere before `--`, also as an option's value, but
// not as the value of `-logfile` or `-conffile`, which the start takes with them.
#[test]
fn the_dry_run_and_the_start_are_the_argv_s_words() {
    let words = |argv: &[&str]| {
        let search = validate(Program::Blastp, argv).unwrap();
        (search.dry_run(), search.starts())
    };
    assert_eq!(words(&[]), (false, false));
    assert_eq!(words(&["-evalue", "1"]), (false, false));
    assert_eq!(words(&["-dryrun"]), (true, true));
    assert_eq!(words(&["-evalue", "1", "-dryrun", "-dryrun"]), (true, true));
    assert_eq!(words(&["-logfile", "log.txt"]), (false, true));
    assert_eq!(words(&["-logfile", "-dryrun"]), (false, true));
    assert_eq!(words(&["-conffile", "c.ini"]), (false, true));
    assert_eq!(words(&["-conffile=c.ini"]), (false, true));
    assert_eq!(words(&["-conffile", "-dryrun"]), (false, true));
    // A last `-logfile` or `-conffile` without its value is dropped: nothing to start for.
    assert_eq!(words(&["-evalue", "1", "-logfile"]), (false, false));
    assert_eq!(words(&["-evalue", "1", "-conffile"]), (false, false));
    for program in [
        Program::Blastn,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
        Program::Psiblast,
    ] {
        assert!(
            validate(program, &["-dryrun"]).unwrap().dry_run(),
            "{program}"
        );
    }
}

// No NCBI counterpart: the library's contract (`api`): a search of records without `-logfile`,
// `-conffile` and `-dryrun` reads no NCBI setting, so its call makes no start at its beginning.
#[test]
fn a_records_search_without_the_words_starts_nothing() {
    let mut diagnostics = Vec::<u8>::new();
    let search = validate(Program::Blastp, &["-evalue", "1"]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    assert!(call.settings.is_none());
    assert!(call.log.is_none());
    assert!(!call.dry_run());
    assert!(diagnostics.is_empty());
    // With one of the words the start runs at once, and its settings are kept for the call.
    let dir = scratch("nostart");
    let conffile = dir.join("c.ini");
    std::fs::write(&conffile, "[BLAST]\n").unwrap();
    let search = validate(Program::Blastp, &["-conffile", path_text(&conffile)]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    assert!(matches!(call.settings, Some(Ok(_))));
    // Only tblastx's run messages follow `-logfile` (`run_messages_follow_log`).
    assert!(call.log.is_none());
    assert!(diagnostics.is_empty());
    let log = dir.join("log.txt");
    let search = validate(Program::Tblastx, &["-logfile", path_text(&log)]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    assert!(matches!(call.log, Some(DiagTarget::File(_))));
    assert_eq!(std::fs::read(&log).unwrap(), b"");
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (quoted at `cli::main_entry`):
// a dry run ends after the arguments. The command line's writes nothing and makes no `-out`
// file, whatever its inputs (NCBI BLAST+ 2.17.0, probes `dry`, `dry_noquery`, `dry_old_out`,
// `dry_dbmissing`); `-outfmt` is read by NCBI's `Run`, so a dry run does not look at it.
#[test]
fn a_dry_run_searches_nothing_in_every_call() {
    let dir = scratch("dry");
    let (query, subject) = proteins();
    let search = validate(Program::Blastp, &["-evalue", "1", "-dryrun"]).unwrap();
    let records = Inputs::new(&query, &subject).named("q.faa", "s.faa");
    let database = Inputs::database(&query, "losat-api-start-no-such-db").named("q.faa", "nodb");
    let empty = Inputs::new(&[], &[]);
    for inputs in [records, database, empty] {
        // `run`: nothing is written.
        let (mut report, mut diagnostics) = (Vec::<u8>::new(), Vec::<u8>::new());
        run(
            &search,
            inputs,
            &mut [("6", &mut report), ("0", &mut Vec::<u8>::new())],
            &mut diagnostics,
        )
        .unwrap();
        assert!(report.is_empty() && diagnostics.is_empty());
        // `run_collect`: an empty report per format, in the order asked for, and no HSP; a
        // format that a search refuses is not read.
        let outcome = run_collect(&search, inputs, &["6", "0", "5", "bogus"]).unwrap();
        assert_eq!(
            outcome
                .reports
                .iter()
                .map(|report| report.outfmt.as_str())
                .collect::<Vec<_>>(),
            ["6", "0", "5", "bogus"]
        );
        for report in &outcome.reports {
            assert!(report.bytes.is_empty(), "{}", report.outfmt);
            assert!(report.spans.is_empty() && report.subject_headings.is_empty());
        }
        assert!(outcome.hsps.is_empty());
        assert_eq!(outcome.query_stats, vec![None; inputs.query.len()]);
        assert!(outcome.diagnostics.is_empty());
        // `run_to_file`: no file is made, and one that exists is left as it is.
        let out = dir.join("out.txt");
        run_to_file(&search, inputs, "6", &out, &mut diagnostics).unwrap();
        assert!(!out.exists());
        std::fs::write(&out, "OLD OUT\n").unwrap();
        run_to_file(&search, inputs, "6", &out, &mut diagnostics).unwrap();
        assert_eq!(std::fs::read(&out).unwrap(), b"OLD OUT\n");
        std::fs::remove_file(&out).unwrap();
        assert!(diagnostics.is_empty());
        // `run_batches`: one batch of no query.
        let batches: Vec<_> = run_batches(&search, inputs)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        assert_eq!(batch.query_indices, 0..0);
        assert!(batch.hsps.is_empty() && batch.query_stats.is_empty());
        assert!(batch.diagnostics.is_empty());
    }
    // Without `-dryrun` the same records are searched.
    let search = validate(Program::Blastp, &["-evalue", "1"]).unwrap();
    let outcome = run_collect(&search, records, &["6"]).unwrap();
    assert!(!outcome.reports[0].bytes.is_empty() && !outcome.hsps.is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI BLAST+ 2.17.0 (probe `dry_log`): the `-logfile` file of a dry run is made (appended to),
// and nothing is written to it.
#[test]
fn a_dry_run_makes_the_logfile() {
    let dir = scratch("drylog");
    let (query, subject) = proteins();
    let inputs = Inputs::new(&query, &subject);
    let log = dir.join("log.txt");
    let search = validate(Program::Blastp, &["-dryrun", "-logfile", path_text(&log)]).unwrap();
    let outcome = run_collect(&search, inputs, &["6"]).unwrap();
    assert!(outcome.reports[0].bytes.is_empty() && outcome.diagnostics.is_empty());
    assert_eq!(std::fs::read(&log).unwrap(), b"");
    std::fs::write(&log, "OLD LINE\n").unwrap();
    let batches = run_batches(&search, inputs).unwrap().count();
    assert_eq!(batches, 1);
    assert_eq!(std::fs::read(&log).unwrap(), b"OLD LINE\n");
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (quoted at
// `ncbi_environment::run_messages_follow_log`). NCBI BLAST+ 2.17.0 (probes `log_warn`,
// `log_empty`, `log_nodir`, `log_old`): tblastx's run messages go to the `-logfile` file (appended
// to), the other programs' to the standard error; `-logfile ""` sends them nowhere, a file in a
// missing directory leaves them on the standard error. The library's standard error is the
// call's diagnostics.
#[test]
fn the_logfile_gets_what_the_command_line_writes_there() {
    let dir = scratch("log");
    let log = dir.join("log.txt");
    let nodir = dir.join("nodir").join("log.txt");
    let (proteins_query, proteins_subject) = proteins();
    let (nucleotides_query, nucleotides_subject) = nucleotides();
    for (program, name, (query, subject)) in [
        (
            Program::Tblastx,
            "tblastx",
            (&nucleotides_query, &nucleotides_subject),
        ),
        (
            Program::Blastn,
            "blastn",
            (&nucleotides_query, &nucleotides_subject),
        ),
        (
            Program::Blastp,
            "blastp",
            (&proteins_query, &proteins_subject),
        ),
    ] {
        let to_log = program == Program::Tblastx;
        let warning = few_matches(name);
        let inputs = Inputs::new(query, subject);
        let plain = validate(program, &["-max_target_seqs", "1"]).unwrap();
        let plain = run_collect(&plain, inputs, &["6"]).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&plain.diagnostics),
            warning,
            "{name}"
        );
        let words =
            |value: &str| validate(program, &["-max_target_seqs", "1", "-logfile", value]).unwrap();
        // `run_collect`.
        std::fs::write(&log, "OLD LINE\n").unwrap();
        let logged = run_collect(&words(path_text(&log)), inputs, &["6"]).unwrap();
        assert_eq!(logged.reports, plain.reports, "{name}");
        let (expected_log, expected_diagnostics) = if to_log {
            (format!("OLD LINE\n{warning}"), String::new())
        } else {
            ("OLD LINE\n".to_string(), warning.clone())
        };
        assert_eq!(
            (
                std::fs::read_to_string(&log).unwrap(),
                String::from_utf8_lossy(&logged.diagnostics).into_owned()
            ),
            (expected_log.clone(), expected_diagnostics.clone()),
            "{name} run_collect"
        );
        // `run`.
        std::fs::write(&log, "OLD LINE\n").unwrap();
        let (mut report, mut diagnostics) = (Vec::<u8>::new(), Vec::<u8>::new());
        run(
            &words(path_text(&log)),
            inputs,
            &mut [("6", &mut report)],
            &mut diagnostics,
        )
        .unwrap();
        assert_eq!(report, plain.reports[0].bytes, "{name}");
        assert_eq!(
            (
                std::fs::read_to_string(&log).unwrap(),
                String::from_utf8_lossy(&diagnostics).into_owned()
            ),
            (expected_log.clone(), expected_diagnostics.clone()),
            "{name} run"
        );
        // `run_batches`: the run's messages with the first batch, or in the file.
        std::fs::write(&log, "OLD LINE\n").unwrap();
        let batches: Vec<_> = run_batches(&words(path_text(&log)), inputs)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(
            (
                std::fs::read_to_string(&log).unwrap(),
                String::from_utf8_lossy(&batches[0].diagnostics).into_owned()
            ),
            (expected_log, expected_diagnostics),
            "{name} run_batches"
        );
        // `-logfile ""` and a file in a missing directory.
        for (value, quiet) in [("", to_log), (path_text(&nodir), false)] {
            let outcome = run_collect(&words(value), inputs, &["6"]).unwrap();
            assert_eq!(outcome.reports, plain.reports, "{name} {value:?}");
            let expected = if quiet { "" } else { warning.as_str() };
            assert_eq!(
                String::from_utf8_lossy(&outcome.diagnostics),
                expected,
                "{name} {value:?}"
            );
        }
        assert!(!dir.join("nodir").exists());
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI BLAST+ 2.17.0 (probes `conf_ok`, `conf_bad`, `conf_missing`, `conf_log`): the file of
// `-conffile` is read in place of `<program>.ini`; a malformed one is NCBI's report and its replay
// (exit status 2), one that is not a file NCBI's critical post (exit status 255), where the
// diagnostics go. The library's error is the text that the command line writes to its standard
// error without `-logfile`; the `-logfile` file gets its part, as on the command line.
#[test]
fn the_conffile_is_loaded_as_on_the_command_line() {
    let dir = scratch("conf");
    let (query, subject) = proteins();
    let inputs = Inputs::new(&query, &subject);
    let plain = run_collect(&validate(Program::Blastp, &[]).unwrap(), inputs, &["6"]).unwrap();
    let conffile = dir.join("c.ini");
    let with = |extra: &[&str]| {
        let mut argv = vec!["-conffile", path_text(&conffile)];
        argv.extend_from_slice(extra);
        validate(Program::Blastp, &argv).unwrap()
    };
    // A file that loads, also one that turns NCBI's data loaders off (the records read the same,
    // `CallStart::settings`).
    for registry in ["[BLAST]\n", "[BLAST]\nDATA_LOADERS = none\n"] {
        std::fs::write(&conffile, registry).unwrap();
        let read = run_collect(&with(&[]), inputs, &["6"]).unwrap();
        assert_eq!(read.reports, plain.reports, "{registry:?}");
        assert!(read.diagnostics.is_empty(), "{registry:?}");
    }
    // A malformed file.
    std::fs::write(&conffile, "[x\n").unwrap();
    let report = Error::Search(SECTION_REPORT.repeat(2).trim_end().to_string());
    assert_eq!(run_collect(&with(&[]), inputs, &["6"]), Err(report.clone()));
    assert_eq!(
        run_collect(&with(&["-dryrun"]), inputs, &["6"]),
        Err(report.clone())
    );
    let mut diagnostics = Vec::<u8>::new();
    assert_eq!(
        run(
            &with(&[]),
            inputs,
            &mut [("6", &mut Vec::<u8>::new())],
            &mut diagnostics
        ),
        Err(report.clone())
    );
    assert!(run_batches(&with(&[]), inputs).is_err());
    let out = dir.join("out.txt");
    assert_eq!(
        run_to_file(&with(&[]), inputs, "6", &out, &mut diagnostics),
        Err(report)
    );
    assert!(!out.exists());
    assert!(diagnostics.is_empty());
    // A name that is not a file.
    let missing = dir.join("missing.ini");
    let critical = format!(
        "Critical: (106.23) Registry file \"{}\" not found, aborting.",
        missing.display()
    );
    let search = validate(Program::Blastp, &["-conffile", path_text(&missing)]).unwrap();
    assert_eq!(
        run_collect(&search, inputs, &["6"]),
        Err(Error::Search(critical.clone()))
    );
    // With `-logfile`: the file gets the post, the error is the same text.
    let log = dir.join("log.txt");
    let search = validate(
        Program::Blastp,
        &[
            "-conffile",
            path_text(&missing),
            "-logfile",
            path_text(&log),
        ],
    )
    .unwrap();
    assert_eq!(
        run_collect(&search, inputs, &["6"]),
        Err(Error::Search(critical.clone()))
    );
    assert_eq!(
        std::fs::read_to_string(&log).unwrap(),
        format!("{critical}\n")
    );
    let _ = std::fs::remove_dir_all(&dir);
}

// No NCBI counterpart: LOSAT refuses the settings that change NCBI's output
// (`ApplicationStart::settings`) after the start, before the program runs (`cli::main_entry`'s
// `start.settings()?`), dry run or not. A library call that runs the start at its beginning
// refuses them the same way (`CallStart::new`); the command line's text is this error's after
// `Error: ` (`tests/ncbi_toolkit_file_options.rs`).
#[test]
fn a_started_call_refuses_the_settings_that_the_command_line_refuses() {
    let dir = scratch("refused");
    let (query, subject) = proteins();
    let inputs = Inputs::new(&query, &subject);
    let conffile = dir.join("c.ini");
    std::fs::write(&conffile, "[BLAST]\nLONG_SEQID=1\n").unwrap();
    let with = |extra: &[&str]| {
        let mut argv = vec!["-conffile", path_text(&conffile)];
        argv.extend_from_slice(extra);
        validate(Program::Blastp, &argv).unwrap()
    };
    let refused = run_collect(&with(&["-dryrun"]), inputs, &["6"]).unwrap_err();
    assert!(
        matches!(&refused, Error::Search(reason) if reason.contains("[BLAST] LONG_SEQID")),
        "{refused:?}"
    );
    let search = with(&[]);
    assert_eq!(run_collect(&search, inputs, &["6"]), Err(refused.clone()));
    let mut diagnostics = Vec::<u8>::new();
    assert_eq!(
        run(
            &search,
            inputs,
            &mut [("6", &mut Vec::<u8>::new())],
            &mut diagnostics
        ),
        Err(refused.clone())
    );
    let out = dir.join("out.txt");
    assert_eq!(
        run_to_file(&search, inputs, "6", &out, &mut diagnostics),
        Err(refused.clone())
    );
    assert!(!out.exists());
    assert_eq!(run_batches(&search, inputs).err(), Some(refused));
    assert!(diagnostics.is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2377-2381 (`-subject_loc`
// excludes `-db`, quoted at `blastp::args`). README: a dry run reports an error in its arguments as
// without `-dryrun`; `validate` parses with one kind of subject, each call with its own
// (`Search::command`).
#[test]
fn a_dry_run_parses_with_the_kind_of_its_subject() {
    let dir = scratch("drykind");
    let (query, subject) = proteins();
    let database = Inputs::database(&query, "nr");
    let refused = run_collect(
        &validate(Program::Blastp, &["-subject_loc", "1-10"]).unwrap(),
        database,
        &["6"],
    )
    .unwrap_err();
    assert!(matches!(refused, Error::Argv(_)), "{refused:?}");
    let search = validate(Program::Blastp, &["-subject_loc", "1-10", "-dryrun"]).unwrap();
    assert_eq!(
        run_collect(&search, database, &["6"]).err(),
        Some(refused.clone())
    );
    let mut diagnostics = Vec::<u8>::new();
    assert_eq!(
        run(
            &search,
            database,
            &mut [("6", &mut Vec::<u8>::new())],
            &mut diagnostics
        ),
        Err(refused.clone())
    );
    // `run_to_file`'s file is the parser's `-out`: its usage line names it, as the command line's
    // does when `-out` is given.
    let out = dir.join("out.txt");
    let Error::Argv(message) = &refused else {
        unreachable!()
    };
    let refused_out = Error::Argv(message.replace(
        "Usage: losat blastp -query <PATH> -db <NAME>\n",
        "Usage: losat blastp -query <PATH> -db <NAME> -out <PATH>\n",
    ));
    assert_ne!(refused_out, refused);
    assert_eq!(
        run_to_file(&search, database, "6", &out, &mut diagnostics),
        Err(refused_out)
    );
    assert_eq!(run_batches(&search, database).err(), Some(refused.clone()));
    assert!(diagnostics.is_empty() && !out.exists());
    // The parse comes before the settings that LOSAT refuses, as on the command line
    // (`cli::main_entry`: the parser, then `start.settings()?`).
    let conffile = dir.join("c.ini");
    std::fs::write(&conffile, "[BLAST]\nLONG_SEQID=1\n").unwrap();
    for dry in [&[][..], &["-dryrun"]] {
        let mut argv = vec!["-subject_loc", "1-10", "-conffile", path_text(&conffile)];
        argv.extend_from_slice(dry);
        let search = validate(Program::Blastp, &argv).unwrap();
        assert_eq!(
            run_collect(&search, database, &["6"]).err(),
            Some(refused.clone()),
            "{dry:?}"
        );
    }
    // The words hold with records.
    let records = Inputs::new(&query, &subject);
    let outcome = run_collect(&search, records, &["6"]).unwrap();
    assert!(outcome.reports[0].bytes.is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (quoted at
// `ncbi_environment::run_messages_follow_log`): tblastx's errors go where its diagnostics go, an
// `-out` that cannot be made too (`losat tblastx ... -out nodir/o.txt -logfile log.txt` writes
// its argument error to `log.txt`); blastn's go to its standard error, the error here.
#[test]
fn an_out_file_that_cannot_be_made_is_an_error_of_the_run() {
    let dir = scratch("outlog");
    let (query, subject) = nucleotides();
    let inputs = Inputs::new(&query, &subject);
    let log = dir.join("log.txt");
    let out = dir.join("nodir").join("out.txt");
    let native = crate::cli::inaccessible("out", &out);
    let message = native
        .downcast_ref::<crate::cli::NativeError>()
        .unwrap()
        .message
        .clone();
    for (program, logged) in [(Program::Tblastx, message.as_str()), (Program::Blastn, "")] {
        std::fs::write(&log, "").unwrap();
        let search = validate(program, &["-logfile", path_text(&log)]).unwrap();
        let mut diagnostics = Vec::<u8>::new();
        assert_eq!(
            run_to_file(&search, inputs, "6", &out, &mut diagnostics),
            Err(Error::search(&native)),
            "{program}"
        );
        assert_eq!(std::fs::read_to_string(&log).unwrap(), logged, "{program}");
        assert!(diagnostics.is_empty(), "{program}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (quoted at
// `ncbi_environment::run_messages_follow_log`): the command line's one run writes its own
// messages (here `-max_target_seqs 1`'s warning) once, to tblastx's `-logfile` file; the library's
// batches, each a run of its own, write them once too, in the file or with the first batch.
#[test]
fn the_batches_write_the_run_s_messages_once() {
    let dir = scratch("batchlog");
    // Two queries of 10,200 letters: a TBLASTX batch ends with the record that brings it to
    // 10,002 letters, so each query is a batch.
    let mut seed = 12345u32;
    let long: Vec<u8> = (0..10_200)
        .map(|_| {
            seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[(seed >> 16) as usize % 4]
        })
        .collect();
    let query = vec![
        FastaRecord::new("Query_1", b"q1", &long[..]),
        FastaRecord::new("Query_2", b"q2", &long[..]),
    ];
    let subject = vec![FastaRecord::new("Subject_1", b"s", NUCLEOTIDES)];
    let inputs = Inputs::new(&query, &subject);
    let warning = few_matches("tblastx");
    let log = dir.join("log.txt");
    let logged = validate(
        Program::Tblastx,
        &["-max_target_seqs", "1", "-logfile", path_text(&log)],
    )
    .unwrap();
    let batches: Vec<_> = run_batches(&logged, inputs)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(batches.len(), 2);
    assert!(batches.iter().all(|batch| batch.diagnostics.is_empty()));
    assert_eq!(std::fs::read_to_string(&log).unwrap(), warning);
    // `run_collect`'s one search writes the same.
    std::fs::remove_file(&log).unwrap();
    let outcome = run_collect(&logged, inputs, &["6"]).unwrap();
    assert!(outcome.diagnostics.is_empty());
    assert_eq!(std::fs::read_to_string(&log).unwrap(), warning);
    // Without `-logfile`: with the first batch.
    let plain = validate(Program::Tblastx, &["-max_target_seqs", "1"]).unwrap();
    let batches: Vec<_> = run_batches(&plain, inputs)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(
        batches
            .iter()
            .map(|batch| String::from_utf8_lossy(&batch.diagnostics).into_owned())
            .collect::<Vec<_>>(),
        [warning, String::new()]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (quoted at
// `ncbi_environment::run_messages_follow_log`): tblastx reads `-outfmt` in its run, after the
// standard error went to the `-logfile` target, so the command line's error of a value goes to
// the file (`losat tblastx ... -outfmt 99 -logfile log.txt`: `Error: Formatting choice is out of
// range`); blastn's goes to its standard error. The library's error is `Error::Argv`, as without
// `-logfile`.
#[test]
fn an_outfmt_error_goes_where_the_run_s_errors_go() {
    let dir = scratch("outfmtlog");
    let (query, subject) = nucleotides();
    let inputs = Inputs::new(&query, &subject);
    let log = dir.join("log.txt");
    for outfmt in ["99", "abc"] {
        let native = crate::blastinput::app::parse_formatting_string(outfmt).unwrap_err();
        let message = native
            .downcast_ref::<crate::cli::NativeError>()
            .unwrap()
            .message
            .clone();
        let refused = Error::Argv(format!("{native:#}"));
        for (program, logged) in [(Program::Tblastx, message.as_str()), (Program::Blastn, "")] {
            let search = validate(program, &["-logfile", path_text(&log)]).unwrap();
            std::fs::write(&log, "").unwrap();
            assert_eq!(
                run_collect(&search, inputs, &[outfmt]).err(),
                Some(refused.clone()),
                "{program} {outfmt}"
            );
            assert_eq!(std::fs::read_to_string(&log).unwrap(), logged, "{program}");
            std::fs::write(&log, "").unwrap();
            let mut diagnostics = Vec::<u8>::new();
            assert_eq!(
                run(
                    &search,
                    inputs,
                    &mut [(outfmt, &mut Vec::<u8>::new())],
                    &mut diagnostics
                ),
                Err(refused.clone()),
                "{program} {outfmt}"
            );
            let out = dir.join("out.txt");
            assert_eq!(
                run_to_file(&search, inputs, outfmt, &out, &mut diagnostics),
                Err(refused.clone()),
                "{program} {outfmt}"
            );
            assert_eq!(
                std::fs::read_to_string(&log).unwrap(),
                logged.repeat(2),
                "{program}"
            );
            assert!(diagnostics.is_empty() && !out.exists(), "{program}");
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}
