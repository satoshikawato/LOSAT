use super::*;
use crate::blastinput::ncbi_environment::DiagTarget;
use crate::utils::test_dir::TempDir;
use crate::writedb::tests::base_in;

/// The refusals of an index file that NCBI cannot read as a database are LOSAT's own, printed
/// with makeblastdb's prefix and ended with exit status 1, and their text is the one of
/// `SeqDbError` (one producer: `NCBI_CRASHES_ON_IT` ends the one where NCBI crashes).
#[test]
fn a_refused_index_file_is_printed_as_a_refusal() {
    let path = std::path::PathBuf::from("o/db.pin");
    let crashes = MakedbError::Database(SeqDbError::NcbiCrashes {
        path: path.clone(),
        message: "the index file has no bytes".to_string(),
    });
    assert_eq!(crashes.exit_code(), 1);
    let text = crashes.to_string();
    assert!(text.starts_with("losat makeblastdb: o/db.pin: "), "{text}");
    assert!(text.ends_with(crate::seqdb::NCBI_CRASHES_ON_IT), "{text}");
    let depends = MakedbError::Database(SeqDbError::Unsupported {
        path,
        message: "a field past the mapping".to_string(),
    });
    assert_eq!(depends.exit_code(), 1);
    assert!(depends
        .to_string()
        .starts_with("losat makeblastdb: o/db.pin: "));
    // NCBI's own error stays NCBI's: exit status 2, "BLAST Database error: ".
    let ncbi = MakedbError::Database(SeqDbError::Format {
        path: std::path::PathBuf::from("o/db.pin"),
        message: "Error: Not a valid version 4 or 5 database.".to_string(),
    });
    assert_eq!(ncbi.exit_code(), 2);
}

// makedb-v5 DECISIONS 12: `[BLAST] MAX_SEQID_LENGTH` is refused with `-parse_seqids`
// (build_db.cpp:739-747, cited at `check_registry_settings`), the library's argument
// rules are makeblastdb's (makeblastdb.cpp:322-329, cited at `MakedbArgs`).
#[test]
fn registry_and_argument_rules_of_parse_seqids() {
    let parse = |words: &[&str]| -> MakedbArgs {
        let mut argv = vec!["losat", "makeblastdb"];
        argv.extend(words);
        match crate::cli::try_parse_from(argv).unwrap().command {
            crate::cli::Commands::Makedb(args) => args,
            command => panic!("not makeblastdb: {command:?}"),
        }
    };
    let plain = parse(&["-dbtype", "prot"]);
    let parsed = parse(&["-dbtype", "prot", "-parse_seqids"]);
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
    // The FASTA reader is the only reader of the entry.
    let mut copied = parsed.clone();
    copied.input_type = "blastdb".to_string();
    assert_eq!(check_registry_settings(&copied, Ok(Some("100"))), Ok(()));
}

/// The dependencies of `-gi_mask` and `-gi_mask_name` among the others, in the order and
/// with the words of NCBI's argument reader (makeblastdb 2.17.0, survey M4c probes b5/o01-o08:
/// `Error:  (CArgException::eConstraint) <message>` after the usage, status 1), which reads
/// makeblastdb's command line and its library calls (`cli::reader`).
#[test]
fn argument_dependencies_follow_ncbis_order() {
    let required = |name: &str, by: &str| {
        format!("Argument \"{name}\". Must be specified, as it is required by argument:  `{by}'")
    };
    let excluded =
        |name: &str, by: &str| format!("Argument \"{name}\". Incompatible with argument:  `{by}'");
    for (words, message) in [
        (
            &["-mask_data", "m", "-gi_mask", "-mask_id", "foo"][..],
            excluded("mask_id", "gi_mask"),
        ),
        (
            &["-gi_mask", "-gi_mask_name", "x"],
            required("parse_seqids", "gi_mask"),
        ),
        (
            &["-gi_mask_name", "x", "-taxid_map", "tm.txt"],
            required("parse_seqids", "taxid_map"),
        ),
        (
            &["-gi_mask", "-taxid_map", "tm.txt", "-mask_data", "m"],
            required("parse_seqids", "gi_mask"),
        ),
        (
            &[
                "-parse_seqids",
                "-mask_desc",
                "d",
                "-gi_mask_name",
                "x",
                "-gi_mask",
            ],
            required("mask_data", "gi_mask_name"),
        ),
        (
            &["-gi_mask_name", "x"],
            required("mask_data", "gi_mask_name"),
        ),
        (
            &["-parse_seqids", "-gi_mask", "-mask_id", "foo"],
            excluded("mask_id", "gi_mask"),
        ),
        (
            &[
                "-gi_mask",
                "-mask_id",
                "foo",
                "-taxid",
                "9606",
                "-taxid_map",
                "tm.txt",
                "-mask_data",
                "m",
            ],
            excluded("mask_id", "gi_mask"),
        ),
        (
            &["-taxid_map", "m.txt"],
            required("parse_seqids", "taxid_map"),
        ),
        (
            &["-parse_seqids", "-taxid_map", "m.txt", "-taxid", "1"],
            excluded("taxid_map", "taxid"),
        ),
    ] {
        let mut argv = vec!["losat", "makeblastdb", "-dbtype", "nucl"];
        argv.extend(words);
        let error = crate::cli::try_parse_from(argv).unwrap_err();
        let crate::cli::ArgvError::Report(stderr) = &error else {
            panic!("{words:?}: {error}")
        };
        let stderr = String::from_utf8_lossy(stderr);
        assert!(
            stderr.ends_with(&format!("\nError: {message}\nError:  {message}\n")),
            "{words:?}: {stderr}"
        );
    }
    let words = [
        "-parse_seqids",
        "-gi_mask",
        "-mask_data",
        "m",
        "-gi_mask_name",
        "x",
    ];
    let mut argv = vec!["losat", "makeblastdb", "-dbtype", "nucl"];
    argv.extend(words);
    assert!(crate::cli::try_parse_from(argv).is_ok());
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
    let fasta = crate::writedb::tests::parse_makedb_words(&["-dbtype".into(), "nucl".into()]);
    assert!(check_environment(&fasta, &|_: &str| None).is_ok());
    let error = check_environment(&fasta, &with(DONT_SCAN, "1")).unwrap_err();
    assert!(error.contains(DONT_SCAN), "{error}");
    // An empty value is no value for makeblastdb.
    assert!(check_environment(&fasta, &with(DONT_SCAN, "")).is_ok());
    // Read by the source of ASN.1 input only.
    assert!(check_environment(
        &fasta,
        &with("NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS", "1")
    )
    .is_ok());
    // The sequences of a database are copied as stored: the variable changes nothing there.
    let mut blastdb = fasta.clone();
    blastdb.input_type = "blastdb".to_string();
    assert!(check_environment(&blastdb, &with(DONT_SCAN, "1")).is_ok());
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1199-1201 (cited at
// `log_of`): the log appends to the `-logfile` file, where the diagnostics append too.
/// The log of `-logfile F` is F opened by the system, beside the diagnostics' file: where the
/// two are one file, the two streams land in F in the order they are written, after what F
/// held, and nothing goes to the standard error. The empty name and a file in a missing
/// directory have no log.
#[test]
fn the_log_and_the_diagnostics_append_to_one_file_in_order() {
    let tmp = TempDir::unique("writedb-log-one-file");
    let name = base_in(&tmp, "mk.log");
    std::fs::write(&name, b"before\n").unwrap();
    let target = diag_target(Some(std::ffi::OsStr::new(&name)));
    assert!(matches!(target, DiagTarget::File(_)));
    let mut log = log_of(&name).expect("the log of a file that opens");
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
    assert!(stderr.is_empty());
    assert_eq!(
        std::fs::read(&name).unwrap(),
        b"before\nlog 1\ndiagnostic\nlog 2\ndiagnostic 2\n"
    );
    for missing in [String::new(), base_in(&tmp, "nodir/mk.log")] {
        assert!(log_of(&missing).is_none(), "{missing:?}");
    }
}

/// Where NCBI's diagnostics read a `-logfile` name otherwise than the system
/// (`ncbi_environment::log_file_name`), the log and the diagnostics are two files: the log
/// is the name as the system opens it (line WQA's D11 round 2, 2026-11-20; NCBI BLAST+
/// 2.17.0 makeblastdb `-logfile a/../mk.log` with a link `a` to `sub/deep` writes its log to
/// `sub/mk.log` and its diagnostics to `mk.log`, and `-logfile logs/` makes the file `logs`
/// for its diagnostics and stops with "File is not accessible" for the log).
#[cfg(unix)]
#[test]
fn the_log_is_the_name_as_the_system_opens_it() {
    let tmp = TempDir::unique("writedb-log-two-files");
    std::fs::create_dir_all(tmp.join("sub/deep")).unwrap();
    std::os::unix::fs::symlink(tmp.join("sub/deep"), tmp.join("a")).unwrap();
    let name = base_in(&tmp, "a/../mk.log");
    let target = diag_target(Some(std::ffi::OsStr::new(&name)));
    let mut log = log_of(&name).expect("the log of a file that opens");
    target.write_text(b"diagnostic\n", &mut std::io::sink());
    log.write_all(b"log\n").unwrap();
    assert_eq!(std::fs::read(tmp.join("mk.log")).unwrap(), b"diagnostic\n");
    assert_eq!(std::fs::read(tmp.join("sub/mk.log")).unwrap(), b"log\n");
    let name = base_in(&tmp, "logs/");
    assert!(matches!(
        diag_target(Some(std::ffi::OsStr::new(&name))),
        DiagTarget::File(_)
    ));
    assert!(tmp.join("logs").is_file());
    assert!(log_of(&name).is_none());
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

/// A device that records each write and fails the writes whose numbers (from 1) are listed,
/// once each, with `EAGAIN` (a failure that passes).
struct Writes {
    writes: Vec<Vec<u8>>,
    fail: Vec<usize>,
    count: usize,
}

impl Write for Writes {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.count += 1;
        if self.fail.contains(&self.count) {
            return Err(std::io::ErrorKind::WouldBlock.into());
        }
        self.writes.push(buf.to_vec());
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Writes the first lines of a log as `BuildDatabase::new` inserts them, with a title of
/// `title`.
fn write_log_head(log: &mut dyn Write, title: &[u8]) {
    let (time, name) = ("T", "/d/x");
    let _ = write!(log, "\n\nBuilding a new DB, current time: {time}");
    let _ = super::super::build_db::endl(log);
    let _ = write!(log, "New DB name:   {name}");
    let _ = super::super::build_db::endl(log);
    let _ = log
        .write_all(b"New DB title:  ")
        .and_then(|()| log.write_all(title));
    let _ = super::super::build_db::endl(log);
    let _ = write!(log, "Keep MBits: T");
    let _ = super::super::build_db::endl(log);
    let _ = super::super::build_db::endl(log);
}

/// The log on standard output (a pipe or a file) as NCBI's `cout`, over `device`.
fn cout_log(device: &mut dyn Write) -> CoutLog<'_> {
    let stream = crate::cli::NcbiStdio::Cout {
        block: 4096,
        lines: false,
    };
    CoutLog(super::super::cout::cout(device, stream))
}

/// The log on standard output is NCBI's `cout` with a flush at each line's end: a pipe or a
/// file gets each line in one write, a title with a new line in it included (line RG item
/// KC-11b-2).
#[test]
fn each_line_of_the_log_reaches_standard_output_in_one_write() {
    for (title, line) in [
        (b"x.fa".as_slice(), b"New DB title:  x.fa\n".as_slice()),
        (b"t\n", b"New DB title:  t\n\n"),
    ] {
        let mut device = Writes {
            writes: Vec::new(),
            fail: Vec::new(),
            count: 0,
        };
        let mut log = cout_log(&mut device);
        write_log_head(&mut log, title);
        log.flush().unwrap();
        drop(log);
        let lines: Vec<&[u8]> = vec![
            b"\n\nBuilding a new DB, current time: T\n",
            b"New DB name:   /d/x\n",
            line,
            b"Keep MBits: T\n",
            b"\n",
        ];
        assert_eq!(device.writes, lines);
    }
}

/// A write of the log that finds the pipe closed ends the build at the end of its line
/// (`build_db::endl`), and so does every later line, as SIGPIPE ends NCBI's process at that
/// line's flush; a write that fails otherwise (`EAGAIN`) loses its line and ends nothing
/// (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10).
#[cfg(not(windows))]
#[test]
fn a_closed_pipe_ends_the_log_at_its_line() {
    /// A device whose pipe is closed from the write `closed_at` (from 1) on.
    struct Closed {
        count: usize,
        closed_at: usize,
    }
    impl Write for Closed {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.count += 1;
            if self.count >= self.closed_at {
                return Err(std::io::ErrorKind::BrokenPipe.into());
            }
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let lines = ["a", "b", "c"];
    for closed_at in 1..=lines.len() {
        let mut device = Closed {
            count: 0,
            closed_at,
        };
        let mut log = cout_log(&mut device);
        let ends: Vec<bool> = lines
            .iter()
            .map(|line| {
                let _ = log.write_all(line.as_bytes());
                super::super::build_db::endl(&mut log).is_err_and(|error| error.is_closed_pipe())
            })
            .collect();
        // Each line reaches the device in one write, at its flush.
        let expected: Vec<bool> = (1..=lines.len()).map(|line| line >= closed_at).collect();
        assert_eq!(ends, expected, "closed at write {closed_at}");
    }
    let mut device = Writes {
        writes: Vec::new(),
        fail: vec![2],
        count: 0,
    };
    let mut log = cout_log(&mut device);
    for line in lines {
        let _ = log.write_all(line.as_bytes());
        assert!(super::super::build_db::endl(&mut log).is_ok(), "{line}");
    }
}

/// A write of the log that fails loses its line, and the stream writes nothing after it
/// (`badbit`), even where the device takes the later writes (line RG item KC-11b-2).
#[test]
fn after_a_failed_write_the_log_writes_nothing() {
    for failed in 1..=5 {
        let mut device = Writes {
            writes: Vec::new(),
            fail: vec![failed],
            count: 0,
        };
        let mut log = cout_log(&mut device);
        write_log_head(&mut log, b"x.fa");
        log.flush().unwrap();
        drop(log);
        assert_eq!(device.count, failed, "write {failed}");
        assert_eq!(device.writes.len(), failed - 1, "write {failed}");
    }
}
