//! `losat makeprofiledb` whose standard output or `-logfile` is a pipe that its reader has closed,
//! against NCBI BLAST+ 2.17.0 makeprofiledb (probes of 2026-11-26 for line RG item KD).
//!
//! A fresh `-dbtype rps` build writes nothing to standard output, so a closed standard-output pipe
//! builds the full database and exits with status 0 in both NCBI and LOSAT. Overwriting an
//! existing database of the same name deletes its files (`s_DeleteMakeprofileDb`) and then writes
//! `Deleted existing BLAST database with identical name.` with `endl` to `m_LogFile` (standard
//! output unless `-logfile` names a file): when that stream is a closed pipe outside Windows, NCBI
//! 2.17.0 is ended by SIGPIPE (exit status 141, no message) before creating any new volume, and
//! LOSAT stops there with exit status 0 and no message (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
//! exception 10). When `CWriteDB` writes an `ib=...,n4sz=...` line to `cout`
//! (`writedb_convert.cpp:464-466`) on a closed standard-output pipe (even with `-logfile`), NCBI
//! is ended by SIGPIPE inside `AddSequence` without closing the database writer or running
//! `~CMakeProfileDBApp`, and LOSAT stops there with exit status 0 and no message. Non-pipe write
//! failures (`/dev/full`) and `-logfile` naming a regular file finish with exit status 0 in both.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:490-493,1330-1333
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
