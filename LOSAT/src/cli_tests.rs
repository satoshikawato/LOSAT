use super::*;
use std::os::unix::ffi::OsStringExt;

// Line DBID pfoff-1 (DECISIONS 62): `CATCH_ALL` reports a `CSeqDBException` (blastp 2.17.0,
// outfmt 13 of `pP pU` with the prefetch off: `BLAST Database error: DB contains no accession
// info.`, exit 2) before any other `CException` (`Error: ` and `what()`, exit 255).
#[test]
fn catch_all_reports_a_seqdb_exception_as_a_database_error() {
    use crate::ncbi_exception::{NcbiException, ThrowSite};
    let site = crate::seqdb::lookup_throw_site(b"DB contains no accession info.").unwrap();
    let error = NativeError::exception(&NcbiException::new(site, "DB contains no accession info."));
    assert_eq!(error.exit, 2);
    assert_eq!(
        error.message,
        b"LOSAT database error: DB contains no accession info.\n"
    );
    let other = NcbiException::new(
        ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line: 852,
            kind: "CStringException::eConvert",
            function: Some("ncbi::NStr::StringToInt8"),
        },
        "x",
    );
    let error = NativeError::exception(&other);
    assert_eq!(error.exit, 255);
    assert_eq!(
        error.message,
        [b"Error: ".as_slice(), &other.what(), b"\n"].concat()
    );
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-229
// `CATCH_ALL`'s clause of a `CBlastException`: the options' code exits 1, a message that names
// a failed allocation exits 4 with the report of running out of memory, any other exits 3.
#[test]
fn catch_all_reports_a_blast_exception_by_its_code_and_message() {
    use crate::ncbi_exception::{NcbiException, ThrowSite};
    let exception = |kind: &'static str, message: &str| {
        NativeError::exception(&NcbiException::new(
            ThrowSite {
                file: "src/algo/blast/api/blast_options_cxx.cpp",
                line: 1,
                kind,
                function: None,
            },
            message,
        ))
    };
    let options = exception("CBlastException::eInvalidOptions", "Out of memory");
    assert_eq!(
        (options.exit, options.message),
        (1, b"LOSAT options error: Out of memory\n".to_vec())
    );
    for message in ["Out of memory", "x: Failed to allocate 8 bytes"] {
        let error = exception("CBlastException::eCoreBlastError", message);
        assert_eq!(error.exit, 4, "{message}");
        assert_eq!(
            error.message,
            format!("BLAST ran out of memory: {message}\n").into_bytes()
        );
    }
    let engine = exception(
        "CRemoteBlastException::eServiceNotAvailable",
        "out of memory",
    );
    assert_eq!(
        (engine.exit, engine.message),
        (3, b"LOSAT engine error: out of memory\n".to_vec())
    );
}
fn blastn_argv(words: Vec<OsString>) -> Vec<OsString> {
    let mut argv = vec![
        OsString::from("losat"),
        OsString::from("blastn"),
        OsString::from("-query"),
        OsString::from("q.fa"),
        OsString::from("-subject"),
        OsString::from("s.fa"),
    ];
    argv.extend(words);
    argv
}

// Without `main`'s `post_errors_in_run` (a library caller of the command line entries such as
// `api::run_blastn`), the run writes no error text of its own: the text stays in the error.
// No unit test calls `post_errors_in_run`.
#[test]
fn a_run_outside_main_posts_no_error_text() {
    let text = "LOSAT query error: Could not construct seq-id from 'gi|abc'\n";
    let error: anyhow::Error = NativeError {
        exit: 1,
        message: text.to_string().into(),
    }
    .into();
    post_native_error(&error);
    assert!(!posted(text.as_bytes()));
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:266-270 keeps the bytes of a string value; LOSAT
// rejects a value that is not UTF-8 explicitly, in both spellings of an option.
#[test]
fn option_values_that_are_not_utf8_are_rejected_explicitly() {
    for (option, value) in [
        ("-query_loc", vec![0xff]),
        ("-subject_loc", b"1-\xff".to_vec()),
        ("-dust", vec![0xfe]),
    ] {
        let separate = blastn_argv(vec![
            OsString::from(option),
            OsString::from_vec(value.clone()),
        ]);
        let mut inline = format!("{option}=").into_bytes();
        inline.extend(&value);
        let joined = blastn_argv(vec![OsString::from_vec(inline)]);
        for argv in [separate, joined] {
            let error = try_parse_from(argv).unwrap_err();
            let message = error.to_string();
            assert!(matches!(error, ArgvError::Refused(_)), "{message}");
            assert!(
                message.contains(&format!("the value of {option} is not UTF-8"))
                    && message.contains("not supported by LOSAT's BLASTN"),
                "{message}"
            );
        }
    }
}

// makeblastdb keeps its values as UTF-8 strings: a value that is not UTF-8 (a file name of
// `-in`, `-out`, `-mask_data`, `-taxid_map`, `-logfile`, a `-title`) is LOSAT's explicit
// rejection after NCBI's reader, in the words of makeblastdb's other refusals (line LIB-B,
// probe 2026-11-21; line RG DECISIONS K4i-2b-2 (1)).
#[test]
fn makeblastdb_rejects_values_that_are_not_utf8_explicitly() {
    for option in [
        "-in",
        "-out",
        "-mask_data",
        "-taxid_map",
        "-logfile",
        "-title",
    ] {
        let argv = vec![
            OsString::from("losat"),
            OsString::from("makeblastdb"),
            OsString::from("-dbtype"),
            OsString::from("nucl"),
            // `-taxid_map` requires it: the parser's errors come before the rejection.
            OsString::from("-parse_seqids"),
            OsString::from(option),
            OsString::from_vec(b"n\xe9.x".to_vec()),
        ];
        let error = try_parse_from(argv).unwrap_err();
        let message = error.to_string();
        assert!(
            matches!(error, ArgvError::Refused(_)),
            "{option}: {message}"
        );
        assert!(
            message.contains(&format!("the value of {option} is not UTF-8"))
                && message.contains("not supported by LOSAT's makeblastdb"),
            "{message}"
        );
    }
}

// BLASTX reads both ranges (line SX, SXd1), excludes `-subject_loc` with `-db` as the other
// programs do, and rejects an option's value that is not UTF-8 with their text (DECISIONS 10
// (c) of line SX).
#[test]
fn blastx_reads_the_ranges_and_rejects_values_that_are_not_utf8() {
    let parsed = try_parse_from([
        "losat",
        "blastx",
        "-query",
        "q.fa",
        "-subject",
        "s.faa",
        "-query_loc",
        "1-10",
        "-subject_loc",
        "2-5",
    ])
    .unwrap();
    let Commands::Blastx(args) = parsed.command else {
        panic!("a BLASTX command line");
    };
    assert_eq!(args.query_loc.as_deref(), Some("1-10"));
    assert_eq!(args.subject_loc.as_deref(), Some("2-5"));
    let error = try_parse_from([
        "losat",
        "blastx",
        "-query",
        "q.fa",
        "-db",
        "d",
        "-subject_loc",
        "2-5",
    ])
    .unwrap_err();
    // NCBI's reader: `-subject_loc` excludes `-db` (`Incompatible with argument`).
    assert!(
        matches!(&error, ArgvError::Report(text) if String::from_utf8_lossy(text).contains("Incompatible with argument")),
        "{error}"
    );
    for (option, value) in [
        ("-query_loc", vec![0xff]),
        ("-subject_loc", b"1-\xff".to_vec()),
        ("-outfmt", vec![0xfe]),
    ] {
        let mut argv = blastn_argv(vec![OsString::from(option), OsString::from_vec(value)]);
        argv[1] = OsString::from("blastx");
        let message = try_parse_from(argv).unwrap_err().to_string();
        assert!(
            message.contains(&format!("the value of {option} is not UTF-8"))
                && message.contains("not supported by LOSAT's BLASTX"),
            "{message}"
        );
    }
}

// The rejection comes after NCBI's reader, which reads the value's bytes: `-help` and the
// reader's errors come first.
#[test]
fn help_and_parser_errors_come_before_a_value_that_is_not_utf8() {
    let mut argv = vec![OsString::from("losat"), OsString::from("blastn")];
    argv.push(OsString::from("-help"));
    argv.push(OsString::from("-query_loc"));
    argv.push(OsString::from_vec(vec![0xff]));
    let error = try_parse_from(argv).unwrap_err();
    assert!(matches!(error, ArgvError::Printed(_)), "{error}");
    let argv = blastn_argv(vec![
        OsString::from("-query_loc"),
        OsString::from_vec(vec![0xff]),
        OsString::from("-no_such_option"),
    ]);
    let error = try_parse_from(argv).unwrap_err();
    assert!(matches!(error, ArgvError::Report(_)), "{error}");
    assert!(!error.to_string().contains("not UTF-8"), "{error}");
}

// A file name keeps its bytes in both spellings (`check_utf8_file_name` reads it).
#[test]
fn file_names_keep_their_bytes_in_both_spellings() {
    let name = b"q\xff.fa".to_vec();
    let mut inline = b"-query=".to_vec();
    inline.extend(&name);
    let argv = vec![
        OsString::from("losat"),
        OsString::from("blastn"),
        OsString::from_vec(inline),
        OsString::from("-subject"),
        OsString::from("s.fa"),
    ];
    let Commands::Blastn(args) = try_parse_from(argv).unwrap().command else {
        panic!("not blastn");
    };
    assert_eq!(args.query.into_os_string(), OsString::from_vec(name));
}

// The words of the top level (`losat` without a program: `top_level`): the reader
// gives what the command line prints there, `losat --version`'s version on the standard output
// (exit status 0) and the error of a word that names no program on its standard error (exit status
// 2), `help` too (the top level has no `help` subcommand). One case each of line RG item K4i-2b's
// review (T-L4), with the command line's bytes of the same commit.
#[test]
fn the_reader_reads_the_top_levels_words_as_the_command_line() {
    let read = |words: &[&str]| {
        let argv: Vec<&str> = std::iter::once("losat")
            .chain(words.iter().copied())
            .collect();
        let error = try_parse_from(argv).map(|_| ()).unwrap_err();
        let code = top_level(
            &std::iter::once("losat")
                .chain(words.iter().copied())
                .map(OsString::from)
                .collect::<Vec<_>>(),
        )
        .status;
        (error, code)
    };
    assert_eq!(
        read(&["--version"]),
        (
            ArgvError::Printed(format!("losat {}\n", env!("CARGO_PKG_VERSION"))),
            0
        )
    );
    for word in ["help", "nosuch"] {
        assert_eq!(
            read(&[word, "blastp"]),
            (
                ArgvError::Report(
                    format!(
                        "error: unknown option or argument '{word}'; use -help to list the options\n"
                    )
                    .into_bytes()
                ),
                2
            ),
            "{word}"
        );
    }
}
