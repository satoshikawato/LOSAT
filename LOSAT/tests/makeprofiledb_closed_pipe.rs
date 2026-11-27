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
//! What the files hold where the run ends is NCBI's too (probes of 2026-11-27 for line RG item
//! T8c, `a_closed_pipe_leaves_ncbis_bytes_where_the_run_ends`): NCBI's ended process loses what
//! its file streams hold, and LOSAT's file buffers drop it
//! (`writedb::ofstream::end_at_closed_pipe`).
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

/// The post of `d_badchar_in.smp`'s title (a control character), which makeprofiledb writes at
/// each read of the profile.
const BADCHAR: &[u8] = b"Error: (803.7) PssmWithParameters.pssm.query.seq.descr.E.title\n\
                         Bad char [0x1] in string at line 15\n\
                         \x01cd00001, PTS_IIB_man, PTS_IIB, PTS system, Mannose/sorbose\n";

/// FNV-1a 64 of `bytes`.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// A run that ends at a closed standard output, and the files NCBI makeprofiledb 2.17.0 left
/// there (line RG item T8c, `remaining-gaps-20261112/T/t8c/probe`, 2026-11-27: a reader gone
/// before the start, each case run 3 times with the same result).
struct Stop {
    /// The case; its list is `in/l_<name>.lst` (the title of a complete volume).
    name: &'static str,
    /// The profiles of the list (`in/<name>.smp`), each the given number of times.
    list: &'static [(&'static str, usize)],
    words: &'static [&'static str],
    /// Whether standard error holds `BADCHAR` (the post before the one that ended NCBI).
    badchar: bool,
    /// Each file left: name, size and FNV-1a 64. LMDB's data files (`.pdb`, `.ptf`, sized by the
    /// host's pages) and a written `.pin` (the build's date) are compared by name only.
    files: &'static [(&'static str, u64, u64)],
}

const EMPTY: u64 = 0xcbf2_9ce4_8422_2325;

/// LMDB's lock file as NCBI's ended process leaves it (8192 bytes, no transaction has used it).
const LOCK: u64 = 0x36a8_ac0a_36ab_9e66;

/// The point where each case ends: `ib=` = the writer's `cout` line for the profile
/// `d2_n4_more2.smp` (`c2_n4.smp` in delta), written when the next profile is added or at the
/// writer's `Close`; the tie = the post of `d_badchar_in.smp` at its read in the volume, which
/// flushes `cout` holding an exclusion line.
const STOPS: &[Stop] = &[
    // `ib=` at the second profile: the first profile's rows on disk (its offset, 0, held).
    Stop {
        name: "r_n4_first",
        list: &[("d2_n4_more2", 1), ("c1", 1)],
        words: &["-dbtype", "rps"],
        badchar: false,
        files: &[
            ("db.aux", 0, EMPTY),
            ("db.freq", 692, 0xd3cd_6aaa_8a64_9bcf),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 692, 0xc17b_32ea_2639_bd84),
        ],
    },
    // `ib=` at the fourth profile: the rows of three, the offsets of two (the third's held), no
    // `.aux` text and no volume bytes (held).
    Stop {
        name: "r_n4_third",
        list: &[("c1", 1), ("c2", 1), ("d2_n4_more2", 1), ("c3", 1)],
        words: &["-dbtype", "rps"],
        badchar: false,
        files: &[
            ("db.00.phr", 0, EMPTY),
            ("db.00.pin", 0, EMPTY),
            ("db.00.psq", 0, EMPTY),
            ("db.aux", 0, EMPTY),
            ("db.freq", 5404, 0x322e_92b2_7989_b3fe),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 5404, 0x5340_64f7_2616_94af),
        ],
    },
    // The same with `-logfile`: the log's file stays empty.
    Stop {
        name: "r_n4_third",
        list: &[("c1", 1), ("c2", 1), ("d2_n4_more2", 1), ("c3", 1)],
        words: &["-dbtype", "rps", "-logfile", "f.log"],
        badchar: false,
        files: &[
            ("db.00.phr", 0, EMPTY),
            ("db.00.pin", 0, EMPTY),
            ("db.00.psq", 0, EMPTY),
            ("db.aux", 0, EMPTY),
            ("db.freq", 5404, 0x322e_92b2_7989_b3fe),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 5404, 0x5340_64f7_2616_94af),
            ("f.log", 0, EMPTY),
        ],
    },
    // `ib=` in the writer's `Close` after the last profile, before the RPS files close.
    Stop {
        name: "r_n4_last",
        list: &[("c1", 1), ("c2", 1), ("c3", 1), ("d2_n4_more2", 1)],
        words: &["-dbtype", "rps", "-blastdb_version", "4"],
        badchar: false,
        files: &[
            ("db.00.phr", 0, EMPTY),
            ("db.00.pin", 0, EMPTY),
            ("db.00.psq", 0, EMPTY),
            ("db.aux", 0, EMPTY),
            ("db.freq", 7756, 0xa76a_44b1_093a_b110),
            ("db.loo", 0, EMPTY),
            ("db.rps", 7756, 0xdcf9_405b_8710_9f09),
        ],
    },
    // 601 profiles (`-index false`): `.aux` and the volume's files past their buffers (8191 B).
    Stop {
        name: "r_n4_aux",
        list: &[("c1", 600), ("d2_n4_more2", 1), ("c2", 1)],
        words: &["-dbtype", "rps", "-index", "false"],
        badchar: false,
        files: &[
            ("db.00.phr", 386_778, 0x0f8b_61a9_2e84_6938),
            ("db.00.pin", 0, EMPTY),
            ("db.00.psq", 8211, 0x814a_04ad_1af6_b8c5),
            ("db.aux", 8198, 0xdf78_abec_8c8e_052f),
            ("db.freq", 1_548_692, 0x22a5_2389_e0ca_ebfc),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 1_548_692, 0x0bc7_7ada_d9dc_7adb),
        ],
    },
    // `ib=` in the third volume: the first two complete (their `.pog` kept, no clean-up).
    Stop {
        name: "r_n4_vol2",
        list: &[
            ("c1", 1),
            ("c2", 1),
            ("c3", 1),
            ("d2_n4_more2", 1),
            ("c2", 1),
        ],
        words: &["-dbtype", "rps", "-max_smp_vol", "2"],
        badchar: false,
        files: &[
            ("db.00.aux", 71, 0xb1e1_1215_4650_6fdf),
            ("db.00.freq", 2592, 0x20e2_e13d_e51f_599f),
            ("db.00.loo", 524_344, 0x921a_8f87_5782_302c),
            ("db.00.pdb", 32768, 0),
            ("db.00.phr", 642, 0xa75f_d20b_fa07_f754),
            ("db.00.pin", 104, 0),
            ("db.00.pog", 36, 0x4501_919c_e010_9275),
            ("db.00.pos", 27, 0xdb9a_1f64_7a43_b26d),
            ("db.00.pot", 20, 0x4192_1acc_144d_4237),
            ("db.00.psq", 24, 0x70df_83e3_0d5c_67bf),
            ("db.00.ptf", 16384, 0),
            ("db.00.pto", 8, 0x89cd_3129_1d2a_efa4),
            ("db.00.rps", 2592, 0x6dd6_ec7a_e800_e46f),
            ("db.01.aux", 87, 0x2046_edd3_80be_ab85),
            ("db.01.freq", 4500, 0x3c74_43bb_76e9_c47c),
            ("db.01.loo", 524_344, 0x0d4c_c4e1_2af8_d719),
            ("db.01.pdb", 32768, 0),
            ("db.01.phr", 1877, 0xa0bc_e0f0_9f67_109c),
            ("db.01.pin", 112, 0),
            ("db.01.pog", 40, 0x756c_4ac6_c699_b128),
            ("db.01.pos", 46, 0x0662_65eb_537e_12df),
            ("db.01.pot", 32, 0x69c4_6984_6f4e_aa5c),
            ("db.01.psq", 41, 0x3c42_26ba_404e_6fdf),
            ("db.01.ptf", 16384, 0),
            ("db.01.pto", 16, 0xb822_24f2_02a3_d9b4),
            ("db.01.rps", 4500, 0x1b9c_bd43_15ff_1446),
            ("db.02.aux", 0, EMPTY),
            ("db.02.freq", 692, 0xd3cd_6aaa_8a64_9bcf),
            ("db.02.loo", 0, EMPTY),
            ("db.02.pdb", 8192, 0),
            ("db.02.pdb-lock", 8192, LOCK),
            ("db.02.ptf", 8192, 0),
            ("db.02.ptf-lock", 8192, LOCK),
            ("db.02.rps", 692, 0xc17b_32ea_2639_bd84),
        ],
    },
    // The tie at the volume's first profile: the header of `.rps` and `.freq` held.
    Stop {
        name: "d_tie_first",
        list: &[("d_badchar_in", 1), ("d_obsr49", 1)],
        words: &["-dbtype", "delta", "-exclude_invalid", "false"],
        badchar: true,
        files: &[
            ("db.aux", 0, EMPTY),
            ("db.freq", 0, EMPTY),
            ("db.loo", 0, EMPTY),
            ("db.rps", 0, EMPTY),
        ],
    },
    // The tie at the volume's first profile with a header of 8,416 bytes (2,101 profiles,
    // `d_obsr49.smp` excluded): the header's first 8,192 bytes, which did not fit the stream's
    // buffer, are in `.rps` and `.freq`.
    Stop {
        name: "d_tie_bighdr",
        list: &[("d_badchar_in", 1), ("c2", 2100), ("d_obsr49", 1)],
        words: &["-dbtype", "delta", "-exclude_invalid", "false"],
        badchar: true,
        files: &[
            ("db.aux", 0, EMPTY),
            ("db.freq", 8192, 0xf311_da4e_821e_93a8),
            ("db.loo", 0, EMPTY),
            ("db.rps", 8192, 0xf311_da4e_821e_93a8),
        ],
    },
    // The tie at the fourth profile: three profiles' rows, two offsets.
    Stop {
        name: "d_tie_third",
        list: &[
            ("c1", 1),
            ("c2", 1),
            ("c3", 1),
            ("d_badchar_in", 1),
            ("d_obsr49", 1),
        ],
        words: &["-dbtype", "delta", "-exclude_invalid", "false"],
        badchar: true,
        files: &[
            ("db.00.phr", 0, EMPTY),
            ("db.00.pin", 0, EMPTY),
            ("db.00.psq", 0, EMPTY),
            ("db.aux", 0, EMPTY),
            ("db.freq", 7084, 0xb869_86b7_d452_345e),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 7084, 0x94b9_bcb8_48f2_9b68),
        ],
    },
    // The tie in the third volume, the first two complete with their delta files.
    Stop {
        name: "d_tie_third",
        list: &[
            ("c1", 1),
            ("c2", 1),
            ("c3", 1),
            ("d_badchar_in", 1),
            ("d_obsr49", 1),
        ],
        words: &[
            "-dbtype",
            "delta",
            "-exclude_invalid",
            "false",
            "-max_smp_vol",
            "2",
        ],
        badchar: true,
        files: &[
            ("db.00.aux", 71, 0xb1e1_1215_4650_6fdf),
            ("db.00.freq", 2592, 0x20e2_e13d_e51f_599f),
            ("db.00.loo", 524_344, 0x921a_8f87_5782_302c),
            ("db.00.obsr", 24, 0x4725_540c_7de7_456a),
            ("db.00.pdb", 32768, 0),
            ("db.00.phr", 642, 0xa75f_d20b_fa07_f754),
            ("db.00.pin", 112, 0),
            ("db.00.pog", 36, 0x4501_919c_e010_9275),
            ("db.00.pos", 27, 0xdb9a_1f64_7a43_b26d),
            ("db.00.pot", 20, 0x4192_1acc_144d_4237),
            ("db.00.psq", 24, 0x70df_83e3_0d5c_67bf),
            ("db.00.ptf", 16384, 0),
            ("db.00.pto", 8, 0x89cd_3129_1d2a_efa4),
            ("db.00.rps", 2592, 0x6dd6_ec7a_e800_e46f),
            ("db.00.wcounts", 2592, 0x9e37_011c_2a3c_c013),
            ("db.01.aux", 87, 0x2046_edd3_80be_ab85),
            ("db.01.freq", 4500, 0x3c74_43bb_76e9_c47c),
            ("db.01.loo", 524_344, 0x0d4c_c4e1_2af8_d719),
            ("db.01.obsr", 36, 0x367f_abc8_bb22_8e96),
            ("db.01.pdb", 32768, 0),
            ("db.01.phr", 1877, 0xa0bc_e0f0_9f67_109c),
            ("db.01.pin", 120, 0),
            ("db.01.pog", 40, 0x756c_4ac6_c699_b128),
            ("db.01.pos", 46, 0x0662_65eb_537e_12df),
            ("db.01.pot", 32, 0x69c4_6984_6f4e_aa5c),
            ("db.01.psq", 41, 0x3c42_26ba_404e_6fdf),
            ("db.01.ptf", 16384, 0),
            ("db.01.pto", 16, 0xb822_24f2_02a3_d9b4),
            ("db.01.rps", 4500, 0x1b9c_bd43_15ff_1446),
            ("db.01.wcounts", 4500, 0xdf2f_cdf1_8f48_2750),
            ("db.02.aux", 0, EMPTY),
            ("db.02.freq", 0, EMPTY),
            ("db.02.loo", 0, EMPTY),
            ("db.02.rps", 0, EMPTY),
        ],
    },
    // Delta, `ib=` at the next profile.
    Stop {
        name: "d_n4_mid",
        list: &[("d_obsr49", 1), ("c2_n4", 1), ("c1", 1), ("c3", 1)],
        words: &["-dbtype", "delta", "-exclude_invalid", "false"],
        badchar: false,
        files: &[
            ("db.aux", 0, EMPTY),
            ("db.freq", 2152, 0x392e_b140_0340_31f0),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 2152, 0xcca0_4cc9_3340_9532),
        ],
    },
    // Delta, `ib=` at the writer's `Close`: `.wcounts` and `.obsr` written before it.
    Stop {
        name: "d_n4_close",
        list: &[("d_obsr49", 1), ("c1", 1), ("c2_n4", 1)],
        words: &["-dbtype", "delta", "-exclude_invalid", "false"],
        badchar: false,
        files: &[
            ("db.00.phr", 0, EMPTY),
            ("db.00.pin", 0, EMPTY),
            ("db.00.psq", 0, EMPTY),
            ("db.aux", 0, EMPTY),
            ("db.freq", 4724, 0x6266_e969_76f6_6acb),
            ("db.loo", 0, EMPTY),
            ("db.obsr", 36, 0x1080_8983_2b60_61a3),
            ("db.pdb", 8192, 0),
            ("db.pdb-lock", 8192, LOCK),
            ("db.ptf", 8192, 0),
            ("db.ptf-lock", 8192, LOCK),
            ("db.rps", 4724, 0x82c2_4804_0c17_a739),
            ("db.wcounts", 4724, 0xc2b2_8ad5_6e3c_a834),
        ],
    },
];

/// Each run that ends at a closed standard output leaves the files NCBI's process ended by
/// SIGPIPE leaves, byte for byte: the RPS files hold the rows of the profiles added and every
/// offset but the last one's, `.aux` and the volume's files what their buffers wrote, the
/// complete volumes before as they are (line RG item T8c).
#[test]
fn a_closed_pipe_leaves_ncbis_bytes_where_the_run_ends() {
    for (index, stop) in STOPS.iter().enumerate() {
        let dir = unique_dir(&format!("losat_mp_pipe_stop_{index}"));
        let names: Vec<&str> = stop
            .list
            .iter()
            .map(|(name, _)| *name)
            .filter(|name| *name != "c2_n4")
            .collect();
        let mut files: Vec<String> = names.iter().map(|name| format!("{name}.smp")).collect();
        files.push("c2.smp".to_string());
        let files: Vec<&str> = files.iter().map(String::as_str).collect();
        copy_profiles(dir.path(), &files);
        let c2 = fs::read_to_string(dir.path().join("in/c2.smp")).unwrap();
        let needle = "seq-data ncbieaa \"KTPAIRLLDKAKIPYELH\"";
        fs::write(
            dir.path().join("in/c2_n4.smp"),
            c2.replace(needle, "seq-data ncbi4na '4141414141414141414141'H"),
        )
        .unwrap();
        let list = format!("in/l_{}.lst", stop.name);
        let mut text = String::new();
        for (name, count) in stop.list {
            text.push_str(&format!("in/{name}.smp\n").repeat(*count));
        }
        fs::write(dir.path().join(&list), text).unwrap();
        let mut words = vec!["-in", &list, "-out", "db"];
        words.extend_from_slice(stop.words);
        let out = run_makeprofiledb(dir.path(), &words, closed_pipe());
        let case = format!("{} {:?}", stop.name, stop.words);
        assert_eq!(out.status.code(), Some(0), "{case}");
        let want_err: &[u8] = if stop.badchar { BADCHAR } else { b"" };
        assert_eq!(
            out.stderr,
            want_err,
            "{case}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let mut left: Vec<String> = fs::read_dir(dir.path())
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .filter(|name| name != "in")
            .collect();
        left.sort();
        let want: Vec<&str> = stop.files.iter().map(|(name, _, _)| *name).collect();
        assert_eq!(left, want, "{case}");
        for &(name, size, hash) in stop.files {
            let bytes = fs::read(dir.path().join(name)).unwrap();
            let by_name = name.ends_with(".pdb")
                || name.ends_with(".ptf")
                || (name.ends_with(".pin") && size > 0);
            if !by_name {
                assert_eq!(
                    (bytes.len() as u64, fnv1a64(&bytes)),
                    (size, hash),
                    "{case}: {name}"
                );
            }
        }
    }
}
