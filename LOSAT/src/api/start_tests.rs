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

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

const NUCLEOTIDES: &[u8] =
    b"ATGAAAACCGCGTATATTGCGAAACAGCGCCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGCCAGCTGGAAGAACGCCTGGGC";

/// NCBI's report of a registry file with the line `[x`, which ends the program (exit status 2),
/// and its replay (`tests/ncbi_toolkit_file_options.rs`); LOSAT writes it in its own words
/// (`crate::wording::losat_text`).
const NCBI_SECTION_REPORT: &str = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";

/// [`NCBI_SECTION_REPORT`] in LOSAT's words.
fn section_report() -> String {
    crate::wording::losat_text(NCBI_SECTION_REPORT).into_owned()
}

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
fn scratch(name: &str) -> crate::utils::test_dir::TempDir {
    crate::utils::test_dir::TempDir::unique(&format!("api-start-{name}"))
}

fn path_text(path: &Path) -> &str {
    path.to_str().unwrap()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (cited at
// `ncbi_environment::preparse`): `-dryrun` anywhere before `--`, also as an option's value, but
// not as the value of `-logfile` or `-conffile`, which the start takes with them.
#[test]
fn the_dry_run_is_the_argv_s_word() {
    let dry_run = |argv: &[&str]| validate(Program::Blastp, argv).unwrap().dry_run();
    assert!(!dry_run(&[]));
    assert!(!dry_run(&["-evalue", "1"]));
    assert!(dry_run(&["-dryrun"]));
    assert!(dry_run(&["-evalue", "1", "-dryrun", "-dryrun"]));
    assert!(!dry_run(&["-logfile", "log.txt"]));
    assert!(!dry_run(&["-logfile", "-dryrun"]));
    assert!(!dry_run(&["-conffile", "c.ini"]));
    assert!(!dry_run(&["-conffile=c.ini"]));
    assert!(!dry_run(&["-conffile", "-dryrun"]));
    assert!(!dry_run(&["-evalue", "1", "-logfile"]));
    assert!(!dry_run(&["-evalue", "1", "-conffile"]));
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

// Line RG D3-1 review 2 N1 and review 4 B1: every library call is a run of its own, whose posts no
// other run's failed diagnostic stream mutes (`ncbi_exception::main_diag_failbit`), and whose
// failed stream mutes no other run's: two calls entered in turn on one thread each read their own.
#[test]
fn a_call_begins_a_diagnostic_stream_of_its_own() {
    let mut diagnostics = Vec::<u8>::new();
    let search = validate(Program::Blastp, &["-evalue", "1"]).unwrap();
    let first = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    let second = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    {
        let _first = first.enter();
        crate::ncbi_exception::set_main_diag_failbit();
        assert!(crate::ncbi_exception::main_diag_failbit());
    }
    assert!(!crate::ncbi_exception::main_diag_failbit());
    let _second = second.enter();
    assert!(!crate::ncbi_exception::main_diag_failbit());
    drop(_second);
    let _first = first.enter();
    assert!(crate::ncbi_exception::main_diag_failbit());
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (cited at
// `ncbi_environment::ApplicationStart`: the application loads its configuration before its
// arguments, whatever they are)
// Every search call runs the start at its beginning (line LIB, errpar fix 3): a search of records
// without `-logfile`, `-conffile` and `-dryrun` too, which writes nothing where the command line
// writes nothing and keeps the settings for the call.
#[test]
fn every_search_call_starts_at_its_beginning() {
    let mut diagnostics = Vec::<u8>::new();
    let search = validate(Program::Blastp, &["-evalue", "1"]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    assert!(call.log.is_none());
    assert!(!call.dry_run());
    assert!(diagnostics.is_empty());
    // tblastx's run messages follow the start's target (`run_messages_follow_log`): `cerr`
    // without `-logfile`, which is the call's diagnostics, as for the other programs.
    let search = validate(Program::Tblastx, &["-evalue", "1"]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    assert!(call.log.is_none());
    assert!(diagnostics.is_empty());
    // With `-conffile` the start loads the file, and its settings are kept for the call.
    let dir = scratch("everycall");
    let conffile = dir.join("c.ini");
    std::fs::write(&conffile, "[BLAST]\n").unwrap();
    let search = validate(Program::Blastp, &["-conffile", path_text(&conffile)]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    // Only tblastx's run messages follow `-logfile` (`run_messages_follow_log`).
    assert!(call.log.is_none());
    assert!(diagnostics.is_empty());
    let log = dir.join("log.txt");
    let search = validate(Program::Tblastx, &["-logfile", path_text(&log)]).unwrap();
    let call = CallStart::new(&search, Inputs::new(&[], &[]), &mut diagnostics).unwrap();
    assert!(matches!(call.log, Some(DiagTarget::File(_))));
    assert_eq!(std::fs::read(&log).unwrap(), b"");
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (cited at `cli::main_entry`):
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
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (cited at
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
    let report = Error::Search(section_report().repeat(2).trim_end().into());
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
        Err(Error::Search(critical.clone().into()))
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
        Err(Error::Search(critical.clone().into()))
    );
    assert_eq!(
        std::fs::read_to_string(&log).unwrap(),
        format!("{critical}\n")
    );
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
        matches!(&refused, Error::Search(_)) && refused.to_string().contains("[BLAST] LONG_SEQID"),
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
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2377-2381 (`-subject_loc`
// excludes `-db`, cited at `blastp::args`). README: a dry run reports an error in its arguments as
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
    // `run_to_file`'s file is the reader's `-out`: NCBI's usage is the same with it.
    let out = dir.join("out.txt");
    assert!(matches!(refused, Error::Argv(_)), "{refused:?}");
    assert_eq!(
        run_to_file(&search, database, "6", &out, &mut diagnostics),
        Err(refused.clone())
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
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (cited at
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
    let message = String::from_utf8(message).unwrap();
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
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (cited at
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
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:93-95 (cited at
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
        let message = String::from_utf8(message).unwrap();
        let refused = Error::Argv(format!("{native:#}").into());
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
}

/// Registry files in memory, by path (the start of a caller that gives its own files:
/// `ApplicationStart::with_files`).
struct Files(Vec<(&'static str, &'static [u8])>);

impl crate::blastinput::ncbi_environment::RegistryFiles for Files {
    fn is_file(&self, path: &Path) -> bool {
        self.read(path).is_some()
    }

    fn read(&self, path: &Path) -> Option<std::io::Result<Vec<u8>>> {
        // A file of the working directory: `./NAME` and `NAME` are one file.
        let in_working = |path: &Path| -> PathBuf {
            path.components()
                .skip_while(|part| matches!(part, std::path::Component::CurDir))
                .collect()
        };
        self.0
            .iter()
            .find(|(name, _)| in_working(Path::new(name)) == in_working(path))
            .map(|(_, bytes)| Ok(bytes.to_vec()))
    }
}

/// A receiver that keeps the files handed to it.
#[derive(Default)]
struct Kept(std::sync::Mutex<Vec<(Vec<u8>, Vec<u8>)>>);

impl crate::report::run_files::FileReceiver for Kept {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
        self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
        Ok(())
    }
}

// The start that a caller makes (`begin_with`, LOSAT Web's): the file of `-conffile` from the
// caller's files (found by NCBI's search path, `.` first), a load that ends the run with NCBI's report
// and the `-logfile` file holding its part (TBLASTX: the report and its replay, as on the command
// line), the run's messages and NCBI's error in that file in TBLASTX, and the dry run.
#[test]
fn a_callers_start_reads_its_files_and_hands_over_the_log() {
    use crate::blastinput::ncbi_environment::ApplicationStart;
    use crate::report::run_files::RunFiles;
    use std::io::Write as _;
    let report = section_report();
    let files = Files(vec![("./good.ini", b"[BLAST]\n"), ("./bad.ini", b"[x\n")]);
    let start = |program: &str, words: &[&str], kept: &std::sync::Arc<Kept>| {
        let run_files = RunFiles::receiver(kept.clone());
        let mut stderr = Vec::new();
        let started = begin_with(
            ApplicationStart::with_files(program, words, &[], &files),
            words,
            &run_files,
            &mut stderr,
        )
        .unwrap();
        (started, stderr, run_files)
    };
    let kept = std::sync::Arc::new(Kept::default());
    let (started, stderr, run_files) = start(
        "blastn",
        &["-conffile", "good.ini", "-logfile", "a.log", "-dryrun"],
        &kept,
    );
    assert_eq!((started.ended(), stderr.as_slice()), (None, &b""[..]));
    assert!(started.settings().is_ok() && started.dry_run());
    started.finish();
    run_files.end().unwrap();
    let (started, stderr, run_files) = start(
        "tblastx",
        &["-conffile", "bad.ini", "-logfile", "b.log"],
        &kept,
    );
    assert_eq!(started.ended(), Some(report.repeat(2).trim_end()));
    assert!(stderr.is_empty() && !started.dry_run());
    started.finish();
    run_files.end().unwrap();
    // TBLASTX's run messages and its NCBI error go to the file, the others' to `cerr`.
    for (program, name) in [("tblastx", "c.log"), ("blastn", "d.log")] {
        let (mut started, _, run_files) = start(program, &["-logfile", name], &kept);
        let mut cerr = Vec::new();
        started.sink(&mut cerr).write_all(b"Warning: w\n").unwrap();
        started.post(b"BLAST engine error: e\n");
        started.finish();
        run_files.end().unwrap();
        let expected: &[u8] = if program == "tblastx" {
            b""
        } else {
            b"Warning: w\n"
        };
        assert_eq!(cerr, expected, "{program}");
    }
    let (started, _, _) = start("blastn", &["-conffile", "missing.ini"], &kept);
    assert_eq!(
        started.ended(),
        Some("Critical: (106.23) Registry file \"missing.ini\" not found, aborting.")
    );
    assert_eq!(
        *kept.0.lock().unwrap(),
        vec![
            (b"a.log".to_vec(), Vec::new()),
            (b"b.log".to_vec(), report.repeat(2).into_bytes()),
            (
                b"c.log".to_vec(),
                b"Warning: w\nBLAST engine error: e\n".to_vec()
            ),
            (b"d.log".to_vec(), Vec::new()),
        ]
    );
}

// What a caller's run adds to the command line's messages for its arguments (`CallerSink::add`,
// LOSAT Web's other formats) goes to its diagnostics, but where the run's messages go nowhere
// (TBLASTX's `-logfile ""` and `-logfile /dev/null`), and never to the file of `-logfile` (line
// WQA's sweep, D3 review r3, 2026-11-19); the sink says where it went.
#[test]
fn what_a_callers_run_adds_goes_where_the_run_s_messages_can_go() {
    use crate::blastinput::ncbi_environment::{ApplicationStart, TextPlace};
    use crate::report::run_files::RunFiles;
    let kept = std::sync::Arc::new(Kept::default());
    let files = Files(Vec::new());
    for (program, logfile, shown) in [
        ("tblastx", None, true),
        ("tblastx", Some("run.log"), true),
        ("tblastx", Some("-"), true),
        ("tblastx", Some(""), false),
        ("tblastx", Some("/dev/null"), false),
        ("blastn", Some(""), true),
        ("blastn", Some("/dev/null"), true),
        ("blastp", Some("run.log"), true),
    ] {
        let words: Vec<&str> = logfile.map_or(Vec::new(), |name| vec!["-logfile", name]);
        let mut stderr = Vec::new();
        let run_files = RunFiles::receiver(kept.clone());
        let mut started = begin_with(
            ApplicationStart::with_files(program, &words, &[], &files),
            &words,
            &run_files,
            &mut stderr,
        )
        .unwrap();
        let mut diagnostics = Vec::new();
        let place = started.sink(&mut diagnostics).add(b"Warning: w\n");
        let expected: &[u8] = if shown { b"Warning: w\n" } else { b"" };
        assert_eq!(diagnostics, expected, "{program} {logfile:?}");
        let expected = if shown {
            TextPlace::StandardError
        } else {
            TextPlace::Nowhere
        };
        assert_eq!(place, expected, "{program} {logfile:?}");
        started.finish();
        run_files.end().unwrap();
    }
    // The files of `-logfile` hold nothing of it.
    assert!(kept
        .0
        .lock()
        .unwrap()
        .iter()
        .all(|(_, bytes)| bytes.is_empty()));
}

// Line WQA's sweep, D7 (2026-11-20): the place that a caller's start gives for its run's messages
// (`CallerStart::messages_place`) is where the messages written to its sink end up, for each kind
// of `-logfile` target (none, `-`, a file, `""` and `/dev/null`) in a program whose run's messages
// follow `-logfile` (TBLASTX) and in one whose do not (BLASTN).
#[test]
fn a_callers_start_says_where_its_run_s_messages_go() {
    use crate::blastinput::ncbi_environment::{ApplicationStart, TextPlace};
    use crate::report::run_files::RunFiles;
    use std::io::Write as _;
    let files = Files(Vec::new());
    for program in ["tblastx", "blastn"] {
        for logfile in [
            None,
            Some("-"),
            Some("run.log"),
            Some(""),
            Some("/dev/null"),
        ] {
            let words: Vec<&str> = logfile.map_or(Vec::new(), |name| vec!["-logfile", name]);
            let kept = std::sync::Arc::new(Kept::default());
            let mut stderr = Vec::new();
            let run_files = RunFiles::receiver(kept.clone());
            let mut started = begin_with(
                ApplicationStart::with_files(program, &words, &[], &files),
                &words,
                &run_files,
                &mut stderr,
            )
            .unwrap();
            let place = started.messages_place();
            let mut cerr = Vec::new();
            started.sink(&mut cerr).write_all(b"Warning: w\n").unwrap();
            started.finish();
            run_files.end().unwrap();
            let in_file = kept
                .0
                .lock()
                .unwrap()
                .iter()
                .any(|(_, bytes)| bytes.as_slice() == b"Warning: w\n");
            let label = format!("{program} {logfile:?}");
            assert_eq!(
                (cerr == b"Warning: w\n", in_file),
                match place {
                    TextPlace::StandardError => (true, false),
                    TextPlace::LogFile => (false, true),
                    TextPlace::Nowhere => (false, false),
                },
                "{label}: {place:?}"
            );
            let expected = match (program, logfile) {
                ("tblastx", Some("run.log")) => TextPlace::LogFile,
                ("tblastx", Some("" | "/dev/null")) => TextPlace::Nowhere,
                _ => TextPlace::StandardError,
            };
            assert_eq!(place, expected, "{label}");
        }
    }
}

// Line WQA's sweep, D7 (2026-11-20): the error of an output that a caller's run could not write
// goes where the command line's run writes it, and the place that the sink gives back is where the
// text is. NCBI's message goes where what the run adds goes (`CallerSink::add`): to the standard
// error, also for TBLASTX's `-logfile run.log`, whose file holds what the CLI's file holds for the
// argv, which is not that output's (the CLI's `tblastx -num_alignments 0 -logfile run.log` leaves it
// empty; line WQA's sweep, D3 review r2, 2026-11-19), and nowhere for `""` and `/dev/null`, where
// the CLI's run of that output writes it nowhere (NCBI BLAST+ 2.17.0 and LOSAT's command line,
// `tblastx -outfmt 6 -num_alignments 0`). LOSAT's own refusal goes to the standard error whatever
// `-logfile` is (`cli::main_entry`; LOSAT's `tblastx -outfmt 6` of a subject title that NCBI cannot
// decode writes it there and none to the file, for each of these `-logfile`).
#[test]
fn a_callers_sink_writes_an_output_s_error_where_the_command_line_does() {
    use crate::blastinput::ncbi_environment::{ApplicationStart, TextPlace};
    use crate::report::run_files::RunFiles;
    let files = Files(Vec::new());
    for program in ["tblastx", "blastn"] {
        for logfile in [
            None,
            Some("-"),
            Some("run.log"),
            Some(""),
            Some("/dev/null"),
        ] {
            for (source, text) in [
                (ErrorSource::Ncbi, &b"BLAST query/options error: e\n"[..]),
                (ErrorSource::Losat, b"Error: a refusal\n"),
            ] {
                let words: Vec<&str> = logfile.map_or(Vec::new(), |name| vec!["-logfile", name]);
                let kept = std::sync::Arc::new(Kept::default());
                let mut stderr = Vec::new();
                let run_files = RunFiles::receiver(kept.clone());
                let mut started = begin_with(
                    ApplicationStart::with_files(program, &words, &[], &files),
                    &words,
                    &run_files,
                    &mut stderr,
                )
                .unwrap();
                let messages = started.messages_place();
                let mut cerr = Vec::new();
                let place = started.sink(&mut cerr).error(text, source);
                started.finish();
                run_files.end().unwrap();
                let in_file = kept
                    .0
                    .lock()
                    .unwrap()
                    .iter()
                    .any(|(_, bytes)| bytes.as_slice() == text);
                let label = format!("{program} {logfile:?} {source:?}");
                assert_eq!(
                    (cerr == text, in_file),
                    match place {
                        TextPlace::StandardError => (true, false),
                        TextPlace::LogFile => (false, true),
                        TextPlace::Nowhere => (false, false),
                    },
                    "{label}: {place:?}"
                );
                let expected = match (source, messages) {
                    (ErrorSource::Ncbi, TextPlace::Nowhere) => TextPlace::Nowhere,
                    _ => TextPlace::StandardError,
                };
                assert_eq!(place, expected, "{label}");
            }
        }
    }
}

// Line WQA's sweep, D11 (2026-11-20) and D19: a caller's `-logfile` name is read as the command
// line on Linux reads it in a working directory holding the run's real files
// (`RunFiles::log_target`). NCBI BLAST+ 2.17.0 and LOSAT's command line, with a `.ncbirc` of the
// one byte 0xEF, which NCBI cannot read (the usage report's `Error: (110.4) Error reading the
// registry after line 1: `, written to the standard error at once), and the run's warning
// (`-max_target_seqs 1`; TBLASTX writes it where the diagnostics go): `-logfile -` writes the
// error once to the standard error; `/dev/stderr` is a file there, which gets the error too, so
// the standard error has it twice; `logs/run.log` without `logs` stays on the standard error,
// whereas with `logs` present (`src/run.log` in the crate's directory) it is a file; `/dev/stdout`
// stays named `/dev/stdout`; `logs/` is the file `logs` and ` run.log` is `run.log`; `/dev/fd/2`,
// `/proc/self/fd/2` and `/dev/fd/1` are files on procfs, which reports no free space, so they get
// nothing; `src` names a directory and `/dev/fd/3` a file descriptor that the run has not opened,
// which NCBI cannot open: as `-`. A `-logfile` naming a file of the run's directory is refused.
#[test]
fn a_callers_log_name_may_have_a_directory_or_be_the_standard_error() {
    use crate::blastinput::ncbi_environment::{ApplicationStart, TextPlace};
    use crate::report::run_files::{FileReceiver, RunFiles};
    use std::io::Write as _;
    // The registry file of the working directory: `.ncbirc`, `ncbi.ini` on Windows.
    let registry = if cfg!(windows) {
        "./ncbi.ini"
    } else {
        "./.ncbirc"
    };
    let files = Files(vec![(registry, b"\xef")]);
    let error = "Error: (110.4) Error reading the registry after line 1: \n";
    let warning = "Warning: w\n";
    // Where the name sends the start's messages: `cerr`, a file that is the standard error, a file
    // of the run (under its normalized host name), or nowhere (a file system without free space).
    #[derive(Clone, Copy, PartialEq)]
    enum Kind {
        Stderr,
        StderrFile,
        File(&'static str),
        Nowhere,
    }
    for program in ["blastn", "tblastx"] {
        let follows = program == "tblastx";
        for (logfile, kind) in [
            ("-", Kind::Stderr),
            ("/dev/stderr", Kind::StderrFile),
            ("/proc/self/fd/2", Kind::Nowhere),
            ("/dev/fd/2", Kind::Nowhere),
            ("/dev/fd/1", Kind::Nowhere),
            ("no_such_logs/run.log", Kind::Stderr),
            ("src/run.log", Kind::File("src/run.log")),
            ("src", Kind::Stderr),
            ("src/", Kind::Stderr),
            ("/tmp/run.log", Kind::File("/tmp/run.log")),
            ("/dev/stdout", Kind::File("/dev/stdout")),
            ("logs/", Kind::File("logs")),
            (" run.log", Kind::File("run.log")),
            ("./run.log", Kind::File("run.log")),
            ("logs/.", Kind::Stderr),
            ("/dev/fd/3", Kind::Stderr),
        ] {
            let words = ["-logfile", logfile];
            let kept = std::sync::Arc::new(Kept::default());
            let mut stderr = Vec::new();
            let run_files = RunFiles::receiver(kept.clone());
            let mut started = begin_with(
                ApplicationStart::with_files(program, &words, &[], &files),
                &words,
                &run_files,
                &mut stderr,
            )
            .unwrap();
            let place = started.messages_place();
            started
                .sink(&mut stderr)
                .write_all(warning.as_bytes())
                .unwrap();
            started.finish();
            run_files.end().unwrap();
            let label = format!("{program} -logfile {logfile}");
            let warned = !follows || matches!(kind, Kind::Stderr | Kind::StderrFile);
            let expected_stderr = format!(
                "{}{}",
                error.repeat(if kind == Kind::StderrFile { 2 } else { 1 }),
                if warned { warning } else { "" }
            );
            assert_eq!(
                String::from_utf8(stderr).unwrap(),
                expected_stderr,
                "{label}"
            );
            let expected_files = match kind {
                Kind::File(host_name) => vec![(
                    host_name.as_bytes().to_vec(),
                    format!("{error}{}", if follows { warning } else { "" }).into_bytes(),
                )],
                _ => Vec::new(),
            };
            assert_eq!(*kept.0.lock().unwrap(), expected_files, "{label}");
            let expected_place = match (follows, kind) {
                (false, _) | (true, Kind::Stderr | Kind::StderrFile) => TextPlace::StandardError,
                (true, Kind::File(_)) => TextPlace::LogFile,
                (true, Kind::Nowhere) => TextPlace::Nowhere,
            };
            assert_eq!(place, expected_place, "{label}");
        }
        // A `-logfile` naming a file of the run's directory (registered or on disk) is refused.
        struct WithRegistered(Vec<Vec<u8>>);
        impl FileReceiver for WithRegistered {
            fn receive(&self, _name: &[u8], _bytes: &[u8]) -> std::io::Result<()> {
                Ok(())
            }
            fn registered_files(&self) -> &[Vec<u8>] {
                &self.0
            }
        }
        let receiver = std::sync::Arc::new(WithRegistered(vec![
            b"q.fa".to_vec(),
            b"in/s.fa".to_vec(),
            b".ncbirc".to_vec(),
        ]));
        for existing in ["q.fa", "./q.fa", "in/s.fa", ".ncbirc", "Cargo.toml"] {
            let words = ["-logfile", existing];
            let mut stderr = Vec::new();
            let err = begin_with(
                ApplicationStart::with_files(program, &words, &[], &files),
                &words,
                &RunFiles::receiver(receiver.clone()),
                &mut stderr,
            )
            .err()
            .expect("existing file must be refused");
            assert_eq!(
                err,
                format!(
                    "a -logfile file '{existing}' that is a file of the run's directory (LOSAT \
                     Web reads the run's files and writes none of them) is not supported by \
                     LOSAT's {}",
                    program.to_ascii_uppercase()
                ),
                "{program} {existing}"
            );
            assert_eq!(String::from_utf8(stderr).unwrap(), error);
        }
    }
}

/// NCBI's frozen file `<program>.<extension>` of `tests/fixtures/syntax_logfile` (the standard
/// error, `err`, and the file of `-logfile log`, `log`, of NCBI BLAST+ 2.17.0's
/// `<program> <inputs> -logfile log -bogus 1`; `tests/syntax_cases_cli.rs`), in LOSAT's words.
fn frozen_logfile_case(program: &str, extension: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/syntax_logfile")
        .join(format!("{program}.{extension}"));
    crate::wording::losat_text(&std::fs::read_to_string(path).unwrap()).into_owned()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1086-1112 (cited at
// `ncbi_app::Diag::report_arg_exception`)
// An argument error of a call's argv (line RG item K4i-2b round 2, R-H1): the call's start opens
// the file of `-logfile`, and the reader posts there, as the command line's does (NCBI's case
// `<program>.log.unk`): the error is what the command line writes to its standard error then, the
// usage with NCBI's line that points to the error logs, without its last line feeds
// (`Error::written_argv`), and the file has the posts. `validate`, which opens no file, gives the
// report as the command line writes it without `-logfile`, its posts included, and makes no file.
#[test]
fn an_argument_error_posts_to_the_calls_logfile() {
    let dir = scratch("argv-logfile");
    let log = dir.join("log");
    for program in Program::ALL {
        let name = program.name();
        let words = ["-logfile", path_text(&log), "-bogus", "1"];
        let _ = std::fs::remove_file(&log);
        let mut diagnostics = Vec::new();
        let error = match CallStart::new(
            &Search::unchecked(program, &words),
            Inputs::new(&[], &[]),
            &mut diagnostics,
        ) {
            Err(Error::Argv(bytes)) => String::from_utf8(bytes).unwrap(),
            Err(other) => panic!("{name}: {other:?}"),
            Ok(_) => panic!("{name}: the call started"),
        };
        assert_eq!(
            error,
            frozen_logfile_case(name, "err").trim_end_matches('\n'),
            "{name}"
        );
        assert_eq!(
            std::fs::read_to_string(&log).unwrap(),
            frozen_logfile_case(name, "log"),
            "{name}"
        );
        assert!(diagnostics.is_empty(), "{name}");
        std::fs::remove_file(&log).unwrap();
        let Err(Error::Argv(checked)) = validate(program, &words) else {
            panic!("{name}: validate accepted the argv");
        };
        let rule = format!("{}\n\n", "=".repeat(72));
        let usage = frozen_logfile_case(name, "err")
            .replace(
                "Error in command-line arguments. See error logs for more details.\n",
                "",
            )
            .replace(&rule, "");
        assert_eq!(
            String::from_utf8(checked).unwrap(),
            format!("{usage}{rule}{}", frozen_logfile_case(name, "log")).trim_end_matches('\n'),
            "{name}: validate"
        );
        assert!(!log.exists(), "{name}: validate makes no file");
    }
}
