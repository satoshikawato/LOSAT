//! `losat makeprofiledb` whose standard output or `-logfile` is a pipe that its reader has closed,
//! against NCBI BLAST+ 2.17.0 makeprofiledb (probes of 2026-11-26 for line RG item KD).
//!
//! A fresh `-dbtype rps` build writes nothing to standard output, so a closed standard-output pipe
//! builds the full database and exits with status 0 in both NCBI and LOSAT. For `-dbtype delta`,
//! exclusion lines are written without `endl` to `m_LogFile` (`cout` by default or with
//! `-logfile -`, 4096-byte pipe buffer; `OfStream` with `-logfile /dev/stdout`, 8191-byte put
//! area): exclusion text within the stream buffer stays buffered so all volumes are built (14
//! files, exit 0 in LOSAT), whereas exclusion text that overflows the buffer hits the closed pipe
//! before any volume is created (0 files; NCBI SIGPIPE 141, LOSAT exit 0). Overwriting an
//! existing database of the same name deletes its files (`s_DeleteMakeprofileDb`) and then writes
//! `Deleted existing BLAST database with identical name.` with `endl` to `m_LogFile`: when that
//! stream is a closed pipe outside Windows, NCBI 2.17.0 is ended by SIGPIPE (exit status 141, no
//! message) before creating any new volume, and LOSAT stops there with exit status 0 and no
//! message (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 10). When `CWriteDB` writes an
//! `ib=...,n4sz=...` line to `cout` (`writedb_convert.cpp:464-466`) on a closed standard-output
//! pipe (even with `-logfile`), NCBI is ended by SIGPIPE inside `AddSequence` without closing the
//! database writer, flushing buffered `-logfile` bytes, or running `~CMakeProfileDBApp`, and
//! LOSAT stops there with exit status 0, empty stderr, and unflushed `-logfile` bytes discarded.
//! Non-pipe write failures (`/dev/full`) and `-logfile` naming a regular file finish with exit
//! status 0 in both.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:490-493,1330-1333,1804-1816
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:464-466

#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

mod temp_support;
use temp_support::unique_dir;

fn fixtures_in() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/makeprofiledb/in")
}

fn closed_pipe() -> Stdio {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

fn run_makeprofiledb(dir: &Path, words: &[&str], stdout: Stdio) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_losat"));
    cmd.current_dir(dir)
        .arg("makeprofiledb")
        .args(words)
        .env("HOME", dir)
        .env("BLAST_USAGE_REPORT", "false")
        .env("LOSAT_BLASTX_WORKER_LOG", "/dev/null")
        .env_remove("BLASTDB")
        .env_remove("NCBI_CONFIG_PATH")
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(Stdio::piped());
    cmd.output().expect("run losat makeprofiledb")
}

fn copy_profiles(dir: &Path, names: &[&str]) {
    let in_dir = dir.join("in");
    fs::create_dir_all(&in_dir).unwrap();
    for name in names {
        fs::copy(fixtures_in().join(name), in_dir.join(name)).unwrap();
    }
}

fn db_files(dir: &Path, prefix: &str) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with(prefix))
        .collect();
    names.sort();
    names
}

#[test]
fn fresh_build_completes_with_a_closed_stdout_pipe() {
    let dir = unique_dir("losat_mp_pipe_fresh");
    copy_profiles(dir.path(), &["l_p3.lst", "pc.smp", "pe.smp", "pf.smp"]);
    let out = run_makeprofiledb(
        dir.path(),
        &["-in", "in/l_p3.lst", "-out", "db1", "-dbtype", "rps"],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert!(
        out.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let files = db_files(dir.path(), "db1.");
    assert!(files.contains(&"db1.rps".to_string()), "{files:?}");
    assert!(files.contains(&"db1.pin".to_string()), "{files:?}");
}

#[test]
fn overwrite_on_closed_stdout_pipe_stops_after_deleting_existing_db() {
    let dir = unique_dir("losat_mp_pipe_overwrite");
    copy_profiles(dir.path(), &["l_p3.lst", "pc.smp", "pe.smp", "pf.smp"]);

    // Build initial single-volume and multi-volume databases.
    let init = run_makeprofiledb(
        dir.path(),
        &["-in", "in/l_p3.lst", "-out", "db2", "-dbtype", "rps"],
        Stdio::piped(),
    );
    assert_eq!(init.status.code(), Some(0));
    assert!(!db_files(dir.path(), "db2.").is_empty());

    // Overwriting with a closed stdout pipe deletes db2.* and stops at the log line.
    let out = run_makeprofiledb(
        dir.path(),
        &["-in", "in/l_p3.lst", "-out", "db2", "-dbtype", "rps"],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert!(
        out.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        db_files(dir.path(), "db2.").is_empty(),
        "expected old DB deleted and no new files, found {:?}",
        db_files(dir.path(), "db2.")
    );

    // Same for `-logfile -` and `-logfile /dev/stdout`.
    for flag in ["-", "/dev/stdout"] {
        let init = run_makeprofiledb(
            dir.path(),
            &["-in", "in/l_p3.lst", "-out", "db_lf", "-dbtype", "rps"],
            Stdio::piped(),
        );
        assert_eq!(init.status.code(), Some(0));
        let out = run_makeprofiledb(
            dir.path(),
            &[
                "-in",
                "in/l_p3.lst",
                "-out",
                "db_lf",
                "-dbtype",
                "rps",
                "-logfile",
                flag,
            ],
            closed_pipe(),
        );
        assert_eq!(out.status.code(), Some(0), "logfile={flag}");
        assert!(out.stderr.is_empty(), "logfile={flag}");
        assert!(
            db_files(dir.path(), "db_lf.").is_empty(),
            "logfile={flag}: {:?}",
            db_files(dir.path(), "db_lf.")
        );
    }

    // Multi-volume overwrite with closed stdout pipe also deletes old volumes and stops.
    let mv_init = run_makeprofiledb(
        dir.path(),
        &[
            "-in",
            "in/l_p3.lst",
            "-out",
            "db_mv",
            "-max_smp_vol",
            "2",
            "-dbtype",
            "rps",
        ],
        Stdio::piped(),
    );
    assert_eq!(mv_init.status.code(), Some(0));
    assert!(db_files(dir.path(), "db_mv").contains(&"db_mv.pal".to_string()));
    let mv_out = run_makeprofiledb(
        dir.path(),
        &[
            "-in",
            "in/l_p3.lst",
            "-out",
            "db_mv",
            "-max_smp_vol",
            "2",
            "-dbtype",
            "rps",
        ],
        closed_pipe(),
    );
    assert_eq!(mv_out.status.code(), Some(0));
    assert!(mv_out.stderr.is_empty());
    assert!(
        db_files(dir.path(), "db_mv").is_empty(),
        "{:?}",
        db_files(dir.path(), "db_mv")
    );
}

#[test]
fn overwrite_with_regular_logfile_or_dev_full_completes() {
    let dir = unique_dir("losat_mp_pipe_logfile");
    copy_profiles(dir.path(), &["l_p3.lst", "pc.smp", "pe.smp", "pf.smp"]);

    let init = run_makeprofiledb(
        dir.path(),
        &["-in", "in/l_p3.lst", "-out", "db_log", "-dbtype", "rps"],
        Stdio::piped(),
    );
    assert_eq!(init.status.code(), Some(0));
    let out = run_makeprofiledb(
        dir.path(),
        &[
            "-in",
            "in/l_p3.lst",
            "-out",
            "db_log",
            "-dbtype",
            "rps",
            "-logfile",
            "build.log",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert!(out.stderr.is_empty());
    assert!(db_files(dir.path(), "db_log.").contains(&"db_log.rps".to_string()));
    assert_eq!(
        fs::read_to_string(dir.path().join("build.log")).unwrap(),
        "Deleted existing BLAST database with identical name.\n"
    );

    if let Ok(full) = fs::OpenOptions::new().write(true).open("/dev/full") {
        let out = run_makeprofiledb(
            dir.path(),
            &["-in", "in/l_p3.lst", "-out", "db_log", "-dbtype", "rps"],
            Stdio::from(full),
        );
        assert_eq!(out.status.code(), Some(0));
        assert!(out.stderr.is_empty());
        assert!(db_files(dir.path(), "db_log.").contains(&"db_log.rps".to_string()));
    }
}

#[test]
fn writer_cout_on_closed_stdout_pipe_stops_without_cleanup_even_with_logfile() {
    let dir = unique_dir("losat_mp_pipe_cout");
    copy_profiles(dir.path(), &["d2_n4_more2.smp", "pc.smp"]);
    fs::write(
        dir.path().join("list.pn"),
        "in/d2_n4_more2.smp\nin/pc.smp\n",
    )
    .unwrap();

    for (db, extra) in [("db_pipe", &[][..]), ("db_log", &["-logfile", "f.log"][..])] {
        let mut args = vec!["-in", "list.pn", "-out", db, "-dbtype", "rps"];
        args.extend_from_slice(extra);
        let out = run_makeprofiledb(dir.path(), &args, closed_pipe());
        assert_eq!(out.status.code(), Some(0), "{db}");
        assert!(
            out.stderr.is_empty(),
            "{db}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let files = db_files(dir.path(), &format!("{db}."));
        // SIGPIPE ends NCBI inside AddSequence: RPS files stay, WriteDB was not closed (.pin absent).
        assert!(files.contains(&format!("{db}.rps")), "{db}: {files:?}");
        assert!(!files.contains(&format!("{db}.pin")), "{db}: {files:?}");
    }
}

#[test]
fn help_version_and_dryrun_on_closed_stdout_pipe() {
    let dir = unique_dir("losat_mp_pipe_meta");
    copy_profiles(dir.path(), &["l_p3.lst", "pc.smp", "pe.smp", "pf.smp"]);
    for words in [
        &["-help"][..],
        &["-version"][..],
        &["-in", "in/l_p3.lst", "-out", "db_dry", "-dryrun"][..],
    ] {
        let out = run_makeprofiledb(dir.path(), words, closed_pipe());
        assert_eq!(out.status.code(), Some(0), "{words:?}");
        assert!(out.stderr.is_empty(), "{words:?}");
    }
}

/// D11 (`-dbtype delta -exclude_invalid false -in in/l_D10.lst -out db`) writes one 70-byte
/// exclusion line without `endl`: on a closed stdout pipe (default, `-logfile -`, and
/// `-logfile /dev/stdout`), NCBI 2.17.0 keeps it in its stream buffer and builds all 14 files.
#[test]
fn delta_buffered_exclusion_completes_on_closed_stdout_pipe() {
    let dir = unique_dir("losat_mp_pipe_delta_d11");
    copy_profiles(
        dir.path(),
        &["l_D10.lst", "d_obsr5.smp", "d_obsr49.smp", "c2.smp"],
    );
    let expected_exts = [
        "aux", "freq", "loo", "obsr", "pdb", "phr", "pin", "pos", "pot", "psq", "ptf", "pto",
        "rps", "wcounts",
    ];
    for (db, extra) in [
        ("db_def", &[][..]),
        ("db_dash", &["-logfile", "-"][..]),
        ("db_dev", &["-logfile", "/dev/stdout"][..]),
    ] {
        let mut args = vec![
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            "in/l_D10.lst",
            "-out",
            db,
        ];
        args.extend_from_slice(extra);
        let out = run_makeprofiledb(dir.path(), &args, closed_pipe());
        assert_eq!(out.status.code(), Some(0), "{db}");
        assert!(
            out.stderr.is_empty(),
            "{db}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let files = db_files(dir.path(), &format!("{db}."));
        let want: Vec<String> = expected_exts.iter().map(|e| format!("{db}.{e}")).collect();
        assert_eq!(files, want, "{db}");
    }
}

/// Many excluded delta profiles (70 bytes per line) test NCBI's exact stream buffer thresholds on a
/// closed stdout pipe: `cout` / `-logfile -` builds all 14 files at 58 lines (4060 <= 4096 bytes)
/// and stops with 0 files at 59 lines (4130 > 4096 bytes); `-logfile /dev/stdout` (`OfStream`,
/// 8191-byte put area) builds all 14 files at 117 lines (8190 < 8191) and stops with 0 files at
/// 118 lines (8260 >= 8191).
#[test]
fn delta_exclusion_overflow_matches_ncbi_buffer_boundaries() {
    let dir = unique_dir("losat_mp_pipe_delta_overflow");
    copy_profiles(dir.path(), &["d_obsr49.smp", "c2.smp"]);

    let cases: &[(&str, &[&str], usize, usize)] = &[
        ("def_58", &[], 58, 14),
        ("def_59", &[], 59, 0),
        ("dash_58", &["-logfile", "-"], 58, 14),
        ("dash_59", &["-logfile", "-"], 59, 0),
        ("dev_117", &["-logfile", "/dev/stdout"], 117, 14),
        ("dev_118", &["-logfile", "/dev/stdout"], 118, 0),
    ];
    for &(label, extra, excluded_count, want_files) in cases {
        let list_name = format!("list_{label}.pn");
        let mut content = "in/d_obsr49.smp\n".repeat(excluded_count);
        content.push_str("in/c2.smp\n");
        fs::write(dir.path().join(&list_name), content).unwrap();

        let mut args = vec![
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            &list_name,
            "-out",
            label,
        ];
        args.extend_from_slice(extra);
        let out = run_makeprofiledb(dir.path(), &args, closed_pipe());
        assert_eq!(out.status.code(), Some(0), "{label}");
        assert!(
            out.stderr.is_empty(),
            "{label}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let files = db_files(dir.path(), &format!("{label}."));
        assert_eq!(files.len(), want_files, "{label}: {files:?}");
    }
}

/// Overwriting an existing delta DB with an exclusion list on a closed stdout pipe: the
/// `Deleted existing ...` line with `endl` comes first, so old files are deleted and no new files
/// are built (NCBI 2.17.0: SIGPIPE 141, 0 files).
#[test]
fn delta_overwrite_plus_exclusion_stops_at_delete_line() {
    let dir = unique_dir("losat_mp_pipe_delta_overwrite");
    copy_profiles(
        dir.path(),
        &["l_D10.lst", "d_obsr5.smp", "d_obsr49.smp", "c2.smp"],
    );
    let init = run_makeprofiledb(
        dir.path(),
        &[
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            "in/l_D10.lst",
            "-out",
            "db_ov",
        ],
        Stdio::piped(),
    );
    assert_eq!(init.status.code(), Some(0));
    assert_eq!(db_files(dir.path(), "db_ov.").len(), 14);

    let out = run_makeprofiledb(
        dir.path(),
        &[
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            "in/l_D10.lst",
            "-out",
            "db_ov",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert!(out.stderr.is_empty());
    assert!(
        db_files(dir.path(), "db_ov.").is_empty(),
        "{:?}",
        db_files(dir.path(), "db_ov.")
    );
}

/// A delta build with an excluded profile buffered in a regular `-logfile`, followed by a
/// delta-valid profile with `ncbi4na` length mismatch whose `CWriteDB` `ib=9,n4sz=11` line hits a
/// closed stdout pipe: NCBI 2.17.0 is ended by SIGPIPE inside `AddSequence` without flushing the
/// `-logfile` buffer (leaving the log file empty) and without closing `WriteDB` (`.pin` absent).
#[test]
fn delta_writer_closed_pipe_does_not_flush_regular_logfile() {
    let dir = unique_dir("losat_mp_pipe_delta_logfile_discard");
    copy_profiles(dir.path(), &["d_obsr49.smp", "c2.smp"]);

    // Derive a delta-valid profile from c2.smp (length 18) with ncbi4na having 2 extra bytes (11B)
    // so delta validation accepts it and CWriteDB prints `ib=9,n4sz=11` to cout.
    let c2 = fs::read_to_string(dir.path().join("in/c2.smp")).unwrap();
    let needle = "seq-data ncbieaa \"KTPAIRLLDKAKIPYELH\"";
    assert!(c2.contains(needle));
    let c2_n4 = c2.replace(needle, "seq-data ncbi4na '4141414141414141414141'H");
    fs::write(dir.path().join("in/c2_n4.smp"), c2_n4).unwrap();
    fs::write(
        dir.path().join("list.pn"),
        "in/d_obsr49.smp\nin/c2_n4.smp\n",
    )
    .unwrap();

    // Baseline with normal stdout: stdout gets `ib=9,n4sz=11\n`, norm.log gets the exclusion line.
    let norm = run_makeprofiledb(
        dir.path(),
        &[
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            "list.pn",
            "-out",
            "db_norm",
            "-logfile",
            "norm.log",
        ],
        Stdio::piped(),
    );
    assert_eq!(norm.status.code(), Some(0));
    assert_eq!(norm.stdout, b"ib=9,n4sz=11\n");
    assert_eq!(
        fs::read(dir.path().join("norm.log")).unwrap(),
        b"in/d_obsr49.smp was excluded: due to too few independent observations\n"
    );

    // Closed stdout pipe: stops inside AddSequence with empty stderr, no .pin, and empty pipe.log.
    let pipe = run_makeprofiledb(
        dir.path(),
        &[
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            "list.pn",
            "-out",
            "db_pipe",
            "-logfile",
            "pipe.log",
        ],
        closed_pipe(),
    );
    assert_eq!(pipe.status.code(), Some(0));
    assert!(
        pipe.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&pipe.stderr)
    );
    let files = db_files(dir.path(), "db_pipe.");
    assert!(files.contains(&"db_pipe.rps".to_string()), "{files:?}");
    assert!(files.contains(&"db_pipe.wcounts".to_string()), "{files:?}");
    assert!(files.contains(&"db_pipe.obsr".to_string()), "{files:?}");
    for ext in ["phr", "pin", "psq", "pog"] {
        assert!(!files.contains(&format!("db_pipe.{ext}")), "{files:?}");
    }
    assert_eq!(
        fs::read(dir.path().join("pipe.log")).unwrap(),
        b"",
        "unflushed -logfile buffer must not be written when closed pipe stops the run"
    );
}

/// An excluded delta profile buffered in `cout` followed by a missing profile on a closed stdout
/// pipe: NCBI's `cerr` is tied to `cout`, so flushing `cout` before `ERR_POST` hits SIGPIPE
/// before any error line is written to stderr.
#[test]
fn delta_buffered_exclusion_before_error_on_closed_stdout_pipe() {
    let dir = unique_dir("losat_mp_pipe_delta_tie");
    copy_profiles(dir.path(), &["d_obsr49.smp"]);
    fs::write(
        dir.path().join("list.pn"),
        "in/d_obsr49.smp\nin/missing.smp\n",
    )
    .unwrap();
    let out = run_makeprofiledb(
        dir.path(),
        &[
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-in",
            "list.pn",
            "-out",
            "db_tie",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert!(
        out.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// Runs `losat makeprofiledb` with both standard output and standard error directed to the same
/// file (`> both.txt 2>&1`) and returns `(exit_code, combined_bytes)`.
fn run_makeprofiledb_combined(dir: &Path, words: &[&str]) -> (Option<i32>, Vec<u8>) {
    let both_path = dir.join("both.txt");
    let both = fs::File::create(&both_path).unwrap();
    let both_err = both.try_clone().unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_losat"));
    cmd.current_dir(dir)
        .arg("makeprofiledb")
        .args(words)
        .env("HOME", dir)
        .env("TZ", "UTC")
        .env("BLAST_USAGE_REPORT", "false")
        .env("LOSAT_BLASTX_WORKER_LOG", "/dev/null")
        .env_remove("BLASTDB")
        .env_remove("NCBI_CONFIG_PATH")
        .stdin(Stdio::null())
        .stdout(Stdio::from(both))
        .stderr(Stdio::from(both_err));
    let status = cmd.status().expect("run losat makeprofiledb");
    let bytes = fs::read(&both_path).unwrap();
    (status.code(), bytes)
}

/// D31, D25, D31 with `-logfile -`, and D31 with `-logfile /dev/stdout` under `2>&1` interleave
/// stdout exclusion lines and cerr diagnostics in the exact byte order of NCBI BLAST+ 2.17.0
/// (`cerr` is tied to `cout`, whereas `-logfile /dev/stdout` is not tied to `cout`).
#[test]
fn cerr_tie_matches_ncbi_combined_stream_order() {
    let d31_want = b"in/d_obsr49.smp was excluded: due to too few independent observations\n\
                     Created protein BLAST (alias) database db with 3 sequences\n";
    let d31_dev_want = b"in/d_obsr49.smp was excluded: due to too few independent observations\n";
    let badchar = b"Error: (803.7) PssmWithParameters.pssm.query.seq.descr.E.title\n\
                    Bad char [0x1] in string at line 15\n\
                    \x01cd00001, PTS_IIB_man, PTS_IIB, PTS system, Mannose/sorbose\n";
    let mut d25_want = Vec::new();
    d25_want.extend_from_slice(badchar);
    d25_want.extend_from_slice(
        b"in/d_badchar_ex.smp was excluded: due to too few independent observations\n",
    );
    for _ in 0..3 {
        d25_want.extend_from_slice(badchar);
    }

    for (label, extra, want) in [
        ("d31", &[][..], &d31_want[..]),
        ("d31_dash", &["-logfile", "-"][..], &d31_want[..]),
        (
            "d31_devstdout",
            &["-logfile", "/dev/stdout"][..],
            &d31_dev_want[..],
        ),
    ] {
        let dir = unique_dir(&format!("losat_mp_order_{label}"));
        copy_profiles(
            dir.path(),
            &["l_D29.lst", "d_obsr49.smp", "c1.smp", "c2.smp", "c3.smp"],
        );
        let mut args = vec![
            "-dbtype",
            "delta",
            "-max_smp_vol",
            "2",
            "-in",
            "in/l_D29.lst",
            "-out",
            "db",
        ];
        args.extend_from_slice(extra);
        let (rc, got) = run_makeprofiledb_combined(dir.path(), &args);
        assert_eq!(rc, Some(0), "{label}");
        assert_eq!(
            got,
            want,
            "{label}:\ngot:  {:?}\nwant: {:?}",
            String::from_utf8_lossy(&got),
            String::from_utf8_lossy(want)
        );
    }

    let dir = unique_dir("losat_mp_order_d25");
    copy_profiles(
        dir.path(),
        &["l_D25.lst", "d_badchar_ex.smp", "d_badchar_in.smp"],
    );
    let (rc, got) = run_makeprofiledb_combined(
        dir.path(),
        &["-dbtype", "delta", "-in", "in/l_D25.lst", "-out", "db"],
    );
    assert_eq!(rc, Some(0));
    assert_eq!(
        got,
        d25_want,
        "D25:\ngot:  {:?}\nwant: {:?}",
        String::from_utf8_lossy(&got),
        String::from_utf8_lossy(&d25_want)
    );
}

/// D25 and D31 on a closed stdout pipe stop at the `cerr`-to-`cout` tie flush where NCBI BLAST+
/// 2.17.0 is ended by SIGPIPE (141): D25 leaves 0 database files and only the first Bad char post
/// on stderr (because `cout` was empty at the first post); D31 leaves all 31 files (including both
/// `.pog` files and `db.pal`, no destructor clean-up) and empty stderr.
#[test]
fn delta_cerr_tie_stops_at_ncbi_point_on_closed_stdout_pipe() {
    let badchar = b"Error: (803.7) PssmWithParameters.pssm.query.seq.descr.E.title\n\
                    Bad char [0x1] in string at line 15\n\
                    \x01cd00001, PTS_IIB_man, PTS_IIB, PTS system, Mannose/sorbose\n";

    let d25_dir = unique_dir("losat_mp_pipe_d25");
    copy_profiles(
        d25_dir.path(),
        &["l_D25.lst", "d_badchar_ex.smp", "d_badchar_in.smp"],
    );
    let d25_out = run_makeprofiledb(
        d25_dir.path(),
        &["-dbtype", "delta", "-in", "in/l_D25.lst", "-out", "db"],
        closed_pipe(),
    );
    assert_eq!(d25_out.status.code(), Some(0));
    assert_eq!(d25_out.stderr, badchar);
    assert!(
        db_files(d25_dir.path(), "db.").is_empty(),
        "{:?}",
        db_files(d25_dir.path(), "db.")
    );

    let d31_dir = unique_dir("losat_mp_pipe_d31");
    copy_profiles(
        d31_dir.path(),
        &["l_D29.lst", "d_obsr49.smp", "c1.smp", "c2.smp", "c3.smp"],
    );
    let d31_out = run_makeprofiledb(
        d31_dir.path(),
        &[
            "-dbtype",
            "delta",
            "-max_smp_vol",
            "2",
            "-in",
            "in/l_D29.lst",
            "-out",
            "db",
        ],
        closed_pipe(),
    );
    assert_eq!(d31_out.status.code(), Some(0));
    assert!(
        d31_out.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&d31_out.stderr)
    );
    let d31_files = db_files(d31_dir.path(), "db.");
    assert_eq!(d31_files.len(), 31, "{d31_files:?}");
    assert!(
        d31_files.contains(&"db.00.pog".to_string()),
        "{d31_files:?}"
    );
    assert!(
        d31_files.contains(&"db.01.pog".to_string()),
        "{d31_files:?}"
    );
    assert!(d31_files.contains(&"db.pal".to_string()), "{d31_files:?}");
}

/// D31 and D25 with `-logfile run.log`: diagnostics go to the file (not `cerr`), so there is no
/// tie flush and both runs complete with the exact NCBI 2.17.0 `run.log` bytes and full database
/// files even on a closed stdout pipe.
#[test]
fn delta_regular_logfile_has_no_cerr_tie_and_matches_ncbi() {
    let badchar = b"Error: (803.7) PssmWithParameters.pssm.query.seq.descr.E.title\n\
                    Bad char [0x1] in string at line 15\n\
                    \x01cd00001, PTS_IIB_man, PTS_IIB, PTS system, Mannose/sorbose\n";
    let d25_excl = b"in/d_badchar_ex.smp was excluded: due to too few independent observations\n";
    let mut d25_log_want = Vec::new();
    for _ in 0..4 {
        d25_log_want.extend_from_slice(badchar);
    }
    // `m_LogFile` opens `run.log` at offset 0 and flushes the exclusion line at program end,
    // overwriting the first `d25_excl.len()` bytes of the 4 appended diagnostics (NCBI 2.17.0).
    d25_log_want[..d25_excl.len()].copy_from_slice(d25_excl);

    let d25_dir = unique_dir("losat_mp_logfile_d25");
    copy_profiles(
        d25_dir.path(),
        &["l_D25.lst", "d_badchar_ex.smp", "d_badchar_in.smp"],
    );
    let d25_out = run_makeprofiledb(
        d25_dir.path(),
        &[
            "-dbtype",
            "delta",
            "-in",
            "in/l_D25.lst",
            "-out",
            "db",
            "-logfile",
            "run.log",
        ],
        closed_pipe(),
    );
    assert_eq!(d25_out.status.code(), Some(0));
    assert!(d25_out.stderr.is_empty());
    assert_eq!(db_files(d25_dir.path(), "db.").len(), 14);
    assert_eq!(
        fs::read(d25_dir.path().join("run.log")).unwrap(),
        d25_log_want
    );

    let d31_dir = unique_dir("losat_mp_logfile_d31");
    copy_profiles(
        d31_dir.path(),
        &["l_D29.lst", "d_obsr49.smp", "c1.smp", "c2.smp", "c3.smp"],
    );
    let d31_out = run_makeprofiledb(
        d31_dir.path(),
        &[
            "-dbtype",
            "delta",
            "-max_smp_vol",
            "2",
            "-in",
            "in/l_D29.lst",
            "-out",
            "db",
            "-logfile",
            "run.log",
        ],
        closed_pipe(),
    );
    assert_eq!(d31_out.status.code(), Some(0));
    assert!(d31_out.stderr.is_empty());
    assert_eq!(db_files(d31_dir.path(), "db.").len(), 29);
    assert_eq!(
        fs::read(d31_dir.path().join("run.log")).unwrap(),
        b"in/d_obsr49.smp was excluded: due to too few independent observations\n"
    );
}
