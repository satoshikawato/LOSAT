use super::*;

// makedb-v5 DECISIONS 12: `[BLAST] MAX_SEQID_LENGTH` is refused with `-parse_seqids`
// (build_db.cpp:739-747, cited at `check_registry_settings`), the library's argument
// rules are makeblastdb's (makeblastdb.cpp:322-329, cited at `check_arguments`).
#[test]
fn registry_and_argument_rules_of_parse_seqids() {
    use clap::Parser;
    #[derive(Parser)]
    struct Wrap {
        #[command(flatten)]
        args: MakedbArgs,
    }
    let parse = |words: &[&str]| {
        let mut argv = vec!["makeblastdb"];
        argv.extend(words);
        Wrap::try_parse_from(argv).unwrap().args
    };
    let plain = parse(&["--dbtype", "prot"]);
    let parsed = parse(&["--dbtype", "prot", "--parse_seqids"]);
    assert_eq!(check_registry_settings(&plain, Ok(Some("100"))), Ok(()));
    assert_eq!(check_registry_settings(&plain, Err("unread")), Ok(()));
    assert_eq!(check_registry_settings(&parsed, Ok(None)), Ok(()));
    let error = check_registry_settings(&parsed, Ok(Some("100"))).unwrap_err();
    assert!(error.contains("[BLAST] MAX_SEQID_LENGTH (100)"), "{error}");
    // The values makeblastdb uses (StringToUInt converts them, not 0) are refused; the
    // ones it ignores are accepted (makeblastdb 2.17.0: pr2/fix/n4probe).
    for used in ["5", "+5", "05", "4294967295"] {
        assert!(
            check_registry_settings(&parsed, Ok(Some(used))).is_err(),
            "{used}"
        );
    }
    for ignored in [
        "",
        "0",
        "+0",
        "abc",
        " 5",
        "5 ",
        "-5",
        "+",
        "4294967296",
        "5x",
    ] {
        assert_eq!(
            check_registry_settings(&parsed, Ok(Some(ignored))),
            Ok(()),
            "{ignored:?}"
        );
    }
    let error = check_registry_settings(&parsed, Err("the registry file x sets ...")).unwrap_err();
    assert!(error.starts_with("the registry file x sets ..."), "{error}");
    let mut map = parse(&["--dbtype", "prot"]);
    map.taxid_map = Some("m.txt".to_string());
    let error = check_arguments(&map).unwrap_err();
    assert_eq!(error.exit_code(), 1);
    assert_eq!(
        error.to_string(),
        "Command line argument error: Argument \"parse_seqids\". Must be specified, as it \
             is required by argument:  `taxid_map'"
    );
    map.parse_seqids = true;
    assert!(check_arguments(&map).is_ok());
    map.taxid = Some(1);
    assert_eq!(
        check_arguments(&map).unwrap_err().to_string(),
        "Command line argument error: Argument \"taxid\". Incompatible with argument:  \
             `taxid_map'"
    );
}

#[test]
fn input_names_split_at_spaces_and_quotes() {
    assert_eq!(split_quoted("a.fa"), ["a.fa"]);
    assert_eq!(split_quoted("a.fa b.fa"), ["a.fa", "b.fa"]);
    assert_eq!(split_quoted("\"a b.fa\""), ["a b.fa"]);
    assert_eq!(split_quoted("\"a b.fa\" c.fa"), ["a b.fa", "c.fa"]);
    assert_eq!(split_quoted("  a.fa  "), ["a.fa"]);
    assert!(split_quoted("").is_empty());
}

#[test]
fn the_variable_that_changes_stored_titles_is_refused() {
    const DONT_SCAN: &str = "NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS";
    let with = |name: &'static str, value: &'static str| {
        move |asked: &str| (asked == name).then(|| std::ffi::OsString::from(value))
    };
    assert!(check_environment(&|_: &str| None).is_ok());
    let error = check_environment(&with(DONT_SCAN, "1")).unwrap_err();
    assert!(error.contains(DONT_SCAN), "{error}");
    // An empty value is no value for makeblastdb.
    assert!(check_environment(&with(DONT_SCAN, "")).is_ok());
    // Read by the source of ASN.1 input only.
    assert!(check_environment(&with("NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS", "1")).is_ok());
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1199-1201 (cited at
// `log_of`): the log appends to the `-logfile` file, where the diagnostics append too.
/// The log of `-logfile F` is the diagnostics' open file itself, not a second handle of it
/// (which the WASI targets cannot make): the two streams, both open at once, land in F in the
/// order they are written, after what F held, and nothing goes to the standard error. The
/// empty name and a file in a missing directory have no log.
#[test]
fn the_log_and_the_diagnostics_write_through_one_open_file() {
    let tmp = crate::writedb::tests::TempDir::new("log-one-file");
    let name = tmp.base("mk.log");
    std::fs::write(&name, b"before\n").unwrap();
    let target = diag_target(Some(std::ffi::OsStr::new(&name)));
    let mut opened = None;
    let mut log = log_of(&target, &name, &mut opened).expect("the log of a file that opens");
    let DiagTarget::File(file) = &target else {
        panic!("the target of a file that opens is the file");
    };
    assert!(
        std::ptr::eq(log, file),
        "the log is the target's own handle"
    );
    let mut stderr = Vec::new();
    {
        let mut diag = DiagStream {
            target: &target,
            stderr: &mut stderr,
        };
        log.write_all(b"log 1\n").unwrap();
        diag.write_all(b"diagnostic\n").unwrap();
        log.write_all(b"log 2\n").unwrap();
        diag.write_all(b"diagnostic 2\n").unwrap();
    }
    assert!(opened.is_none(), "no other file is opened");
    assert!(stderr.is_empty());
    assert_eq!(
        std::fs::read(&name).unwrap(),
        b"before\nlog 1\ndiagnostic\nlog 2\ndiagnostic 2\n"
    );
    for missing in [String::new(), tmp.base("nodir/mk.log")] {
        let target = diag_target(Some(std::ffi::OsStr::new(&missing)));
        let mut opened = None;
        assert!(
            log_of(&target, &missing, &mut opened).is_none(),
            "{missing:?}"
        );
    }
}

#[test]
fn errors_end_with_makeblastdb_s_exit_codes_and_prefixes() {
    let io = || std::io::Error::other("no");
    for (error, code, text) in [
        (
            MakedbError::Build(WriteDbError::Options("o".into())),
            1,
            "LOSAT options error: o",
        ),
        (
            MakedbError::Build(WriteDbError::Creation("c".into())),
            1,
            "LOSAT database creation error: c",
        ),
        (
            MakedbError::Database(SeqDbError::NotFound("d".into())),
            2,
            "LOSAT database error: d",
        ),
        (
            MakedbError::Build(WriteDbError::OutOfMemory),
            4,
            "BLAST ran out of memory",
        ),
        (
            MakedbError::Build(WriteDbError::Exception(
                crate::ncbi_exception::NcbiException::new(
                    crate::ncbi_exception::ThrowSite {
                        file: "src/f.cpp",
                        line: 7,
                        kind: "CException::eUnknown",
                        function: None,
                    },
                    "m",
                ),
            )),
            255,
            // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
            &*crate::wording::losat_text(&format!(
                "Error: NCBI C++ Exception:\n    T0 \"{}/src/f.cpp\", line 7: Error: \
                     (CException::eUnknown) m\n",
                crate::ncbi_exception::NCBI_BUILD_ROOT
            )),
        ),
        (
            MakedbError::Build(WriteDbError::Other("x".into())),
            255,
            "Error: x",
        ),
        (
            MakedbError::Build(WriteDbError::Io {
                path: "p".into(),
                source: io(),
            }),
            255,
            "Error: p: no",
        ),
    ] {
        assert_eq!(error.exit_code(), code, "{error}");
        assert_eq!(error.to_string(), text);
    }
}
