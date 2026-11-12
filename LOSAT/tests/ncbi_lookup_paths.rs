//! Where NCBI BLAST+ looks for files of its own around a search (env parity S1 and S4): the
//! registry files (`.ncbirc`, `<program>.ini`) on its search path, which depends on when NCBI
//! makes the path (`ncbi_environment::registry_search_path`), and the file of a `-matrix` that
//! its tables do not have (`ApplicationSettings::matrix_file`), which NCBI 2.17.0 reads and then
//! aborts on (rc 134, `free(): invalid size`; tblastx, and blastx `-ungapped -comp_based_stats 0`)
//! when it is a matrix. And what NCBI writes about a registry file it reports as malformed (env
//! parity G, `ncbi_environment::ApplicationStart`): the messages of each read of `.ncbirc`,
//! where they go with `-logfile`, and the exit status 2 of a malformed `<program>.ini`.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const INPUTS: &str = "tests/fasta/outfmt11";

/// A directory of the test's own, made empty (tests of one binary run in parallel).
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat_ncbi_lookup_paths_{name}_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// `program args` from `dir` with `envs` and nothing else of NCBI's or LOSAT's environment,
/// `HOME` a directory without files, and `stdin` on a pipe.
fn run_in(
    program: &Path,
    args: &[&str],
    dir: &Path,
    envs: &[(&str, &str)],
    stdin: &[u8],
) -> Output {
    let mut command = Command::new(program);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(
                text.as_ref(),
                "BLASTMAT" | "BLASTDB" | "BLAST_USAGE_REPORT" | "BATCH_SIZE" | "CHUNK_SIZE"
            )
        {
            command.env_remove(&name);
        }
    }
    let home = dir.join("home");
    fs::create_dir_all(&home).unwrap();
    command
        .env("HOME", &home)
        .env("APPDATA", &home)
        .env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat_ncbi_lookup_paths_worker.log"),
        )
        .envs(envs.iter().copied())
        .current_dir(dir)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().expect("run LOSAT");
    child.stdin.take().unwrap().write_all(stdin).unwrap();
    child.wait_with_output().expect("LOSAT output")
}

fn input(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(INPUTS)
        .join(name)
        .display()
        .to_string()
}

/// The two programs that take a `-matrix` without a table, with their query and subject.
fn matrix_runs() -> [(&'static str, Vec<String>); 2] {
    [
        (
            "TBLASTX",
            ["tblastx", "-subject", &input("arch_tbx_s.fna")]
                .map(String::from)
                .to_vec(),
        ),
        (
            "BLASTX",
            [
                "blastx",
                "-subject",
                &input("arch_bx_s.faa"),
                "-ungapped",
                "-comp_based_stats",
                "0",
            ]
            .map(String::from)
            .to_vec(),
        ),
    ]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1381
// (`BlastFindMatrixPath`, quoted at `ApplicationSettings::matrix_file`). NCBI BLAST+ 2.17.0
// (`-matrix mymat`, env parity S4): with the file `MYMAT` in `BLASTMAT`, in `BLASTMAT/nt`, in a
// directory of `NCBI_DATA_PATH` or `NCBI_CONFIG__NCBI__DATA`, or in `./data`, tblastx writes its
// first line and two empty ones and aborts (rc 134); with `BLASTMAT` naming a file, or no file,
// it stops with "BLAST engine error: Error: Unknown error code -1" (rc 3); an empty query reads
// no matrix (rc 0). blastx `-ungapped -comp_based_stats 0` does the same.
#[test]
fn a_matrix_file_that_ncbi_reads_is_refused() {
    let dir = scratch("matrix");
    for sub in ["bm/nt", "dp", "data", "cwd/data"] {
        fs::create_dir_all(dir.join(sub)).unwrap();
    }
    let bm = dir.join("bm").display().to_string();
    let dp = dir.join("dp").display().to_string();
    let file = |at: &str| {
        for place in ["bm/MYMAT", "bm/nt/MYMAT", "dp/MYMAT", "data/MYMAT"] {
            let _ = fs::remove_file(dir.join(place));
        }
        if !at.is_empty() {
            fs::write(dir.join(at), "# any text\n").unwrap();
        }
    };
    let data_path = format!("{}{dp}", if cfg!(windows) { ";;" } else { "::" });
    for (program, words) in matrix_runs() {
        let mut args: Vec<&str> = words.iter().map(String::as_str).collect();
        let query = input(if program == "TBLASTX" {
            "arch_tbx_q.fna"
        } else {
            "arch_bx_q.fna"
        });
        args.extend(["-matrix", "mymat", "-num_threads", "1"]);
        let mut with_query = args.clone();
        with_query.extend(["-query", &query]);
        for (at, envs) in [
            ("bm/MYMAT", vec![("BLASTMAT", bm.as_str())]),
            ("bm/nt/MYMAT", vec![("BLASTMAT", bm.as_str())]),
            ("dp/MYMAT", vec![("NCBI_DATA_PATH", data_path.as_str())]),
            ("dp/MYMAT", vec![("NCBI_CONFIG__NCBI__DATA", dp.as_str())]),
            ("data/MYMAT", vec![]),
        ] {
            file(at);
            let output = run_in(Path::new(BIN), &with_query, &dir, &envs, b"");
            let stderr = String::from_utf8_lossy(&output.stderr);
            // TBLASTX refuses before it writes anything; BLASTX, which reads its query as it
            // searches, refuses the first batch where NCBI reads the file, after the prolog, in
            // the frame of the batch's error.
            if program == "TBLASTX" {
                assert_eq!(output.status.code(), Some(1), "{program} {at}: {stderr}");
                assert!(output.stdout.is_empty(), "{program} {at}");
            } else {
                assert_eq!(output.status.code(), Some(3), "{program} {at}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST engine error: -matrix mymat"),
                    "{program} {at}: {stderr}"
                );
            }
            assert!(
                stderr.contains("-matrix mymat is not a matrix of NCBI BLAST+'s tables")
                    && stderr.contains("MYMAT")
                    && stderr.contains(&format!("not supported by LOSAT's {program}")),
                "{program} {at}: {stderr}"
            );
            // No query, no batch: NCBI reads no matrix (an empty pipe: the report without
            // queries).
            let mut empty = args.clone();
            empty.extend(["-query", "-"]);
            let output = run_in(Path::new(BIN), &empty, &dir, &envs, b"");
            assert_eq!(output.status.code(), Some(0), "{program} {at} empty query");
        }
        // BLASTMAT that is a file, and no file at all: NCBI's error of the first batch.
        file("bm/MYMAT");
        let blastmat_file = dir.join("bm/MYMAT").display().to_string();
        for envs in [vec![("BLASTMAT", blastmat_file.as_str())], vec![]] {
            let at = dir.join("cwd");
            let output = run_in(Path::new(BIN), &with_query, &at, &envs, b"");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                output.status.code(),
                Some(3),
                "{program} {envs:?}: {stderr}"
            );
            assert!(
                stderr.starts_with("BLAST engine error: Error: Unknown error code -1"),
                "{program} {envs:?}: {stderr}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1650-1652 (a file found that
// `fopen` cannot open: -1, "Unknown error code -1"). NCBI BLAST+ 2.17.0 with `MYMAT` of mode 000 in
// `BLASTMAT`: rc 3 and that error (env parity S4 review).
#[cfg(unix)]
#[test]
fn a_matrix_file_that_ncbi_cannot_open_gives_ncbis_error() {
    use std::os::unix::fs::PermissionsExt;
    let dir = scratch("unreadable");
    let file = dir.join("MYMAT");
    fs::write(&file, "# any text\n").unwrap();
    fs::set_permissions(&file, fs::Permissions::from_mode(0o000)).unwrap();
    if fs::File::open(&file).is_ok() {
        // The superuser opens it all the same, as NCBI would then.
        let _ = fs::remove_dir_all(&dir);
        return;
    }
    let blastmat = dir.display().to_string();
    for (program, words) in matrix_runs() {
        let query = input(if program == "TBLASTX" {
            "arch_tbx_q.fna"
        } else {
            "arch_bx_q.fna"
        });
        let mut args: Vec<&str> = words.iter().map(String::as_str).collect();
        args.extend(["-matrix", "mymat", "-num_threads", "1", "-query", &query]);
        let output = run_in(Path::new(BIN), &args, &dir, &[("BLASTMAT", &blastmat)], b"");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(3), "{program}: {stderr}");
        assert!(
            stderr.starts_with("BLAST engine error: Error: Unknown error code -1"),
            "{program}: {stderr}"
        );
    }
    let _ = fs::set_permissions(&file, fs::Permissions::from_mode(0o644));
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:331-396 and
// c++/src/algo/blast/api/blast_usage_report.cpp:195-209 (quoted at `registry_search_path`).
// NCBI BLAST+ 2.17.0 (env parity S1, a malformed `.ncbirc` or `blastn.ini`, whose read NCBI
// reports): blastn run by a link on `PATH` reads one in the link's directory only with
// `BLAST_USAGE_REPORT=false` or `=0`; with an empty `NCBI_CONFIG_PATH` it reads none, not even
// one in the working directory.
#[cfg(unix)]
#[test]
fn registry_files_are_looked_for_where_ncbi_looks() {
    let dir = scratch("registry");
    let shim = dir.join("shim");
    fs::create_dir_all(&shim).unwrap();
    std::os::unix::fs::symlink(BIN, shim.join("blastn")).unwrap();
    let path = format!(
        "{}:{}",
        shim.display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let query = input("arch_q.fna");
    let subject = input("arch_s.fna");
    let args = ["-query", &query, "-subject", &subject, "-outfmt", "6"];
    let blastn = Path::new("blastn");
    let malformed = "a line without an equals sign\n";
    // A read `.ncbirc` is reported and the search goes on; a read `blastn.ini` ends it with
    // exit status 2 (`a_malformed_registry_file_is_reported_as_ncbi_reports_it`).
    let refused = |output: &Output, file: &Path| {
        let stderr = String::from_utf8_lossy(&output.stderr);
        if file.ends_with("blastn.ini") {
            output.status.code() == Some(2)
                && stderr.contains("Application's initialization failed")
        } else {
            output.status.code() == Some(0)
                && stderr.contains("Syntax error in system-wide configuration file")
        }
    };
    // The link's directory.
    let in_shim = shim.join(".ncbirc");
    fs::write(&in_shim, malformed).unwrap();
    let plain = run_in(blastn, &args, &dir, &[("PATH", &path)], b"");
    assert_eq!(plain.status.code(), Some(0), "{plain:?}");
    assert!(plain.stderr.is_empty(), "{plain:?}");
    for value in ["false", "0"] {
        let late = run_in(
            blastn,
            &args,
            &dir,
            &[("PATH", &path), ("BLAST_USAGE_REPORT", value)],
            b"",
        );
        assert!(refused(&late, &in_shim), "{value}: {late:?}");
    }
    let on = run_in(
        blastn,
        &args,
        &dir,
        &[("PATH", &path), ("BLAST_USAGE_REPORT", "true")],
        b"",
    );
    assert_eq!(on.status.code(), Some(0), "{on:?}");
    assert!(on.stderr.is_empty(), "{on:?}");
    fs::remove_file(&in_shim).unwrap();
    // `.ncbirc` turned off: the usage report looks for no file, and the path is made later
    // (ncbireg.cpp:1634-1648); `blastn.ini` in the link's directory is read.
    let ini = shim.join("blastn.ini");
    fs::write(&ini, malformed).unwrap();
    let plain = run_in(blastn, &args, &dir, &[("PATH", &path)], b"");
    assert_eq!(plain.status.code(), Some(0), "{plain:?}");
    assert!(plain.stderr.is_empty(), "{plain:?}");
    for variable in ["NCBI_DONT_USE_NCBIRC", "NCBI_CONFIG__NCBI__DONT_USE_NCBIRC"] {
        let late = run_in(blastn, &args, &dir, &[("PATH", &path), (variable, "")], b"");
        assert!(refused(&late, &ini), "{variable}: {late:?}");
    }
    fs::remove_file(&ini).unwrap();
    // The working directory, and an empty NCBI_CONFIG_PATH.
    let in_cwd = dir.join(".ncbirc");
    fs::write(&in_cwd, malformed).unwrap();
    let cwd = run_in(blastn, &args, &dir, &[("PATH", &path)], b"");
    assert!(refused(&cwd, Path::new(".ncbirc")), "{cwd:?}");
    let none = run_in(
        blastn,
        &args,
        &dir,
        &[("PATH", &path), ("NCBI_CONFIG_PATH", "")],
        b"",
    );
    assert_eq!(none.status.code(), Some(0), "{none:?}");
    assert!(none.stderr.is_empty(), "{none:?}");
    let _ = fs::remove_dir_all(&dir);
}

/// The build tree that the official Linux x64 release of NCBI BLAST+ 2.17.0 names in its `T0`
/// lines (as `ncbi_exception::NCBI_BUILD_ROOT`). The registry tests that use it and the helpers
/// below run on Unix only.
#[cfg(unix)]
const NCBI_BUILD_ROOT: &str = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";

/// NCBI's post of one read of a `.ncbirc` with the line `garbage` at `line`.
#[cfg(unix)]
fn critical(line: usize) -> String {
    format!(
        "Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'garbage' (m_Pos = {line})\n\n"
    )
}

/// NCBI's report of a `<program>.ini` of the line `garbage`, which ends the program.
#[cfg(unix)]
const INI_REPORT: &str = "Error: (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\n";

/// The name of the system-wide registry file on Unix (Windows names it `ncbi.ini`).
#[cfg(unix)]
const NCBIRC: &str = ".ncbirc";

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1650-1663 and
// c++/src/corelib/ncbiapp.cpp:669-674 (quoted at `ncbi_environment::ApplicationStart`). NCBI
// BLAST+ 2.17.0, the official Linux x64 release (env parity g_reg probes of 2026-11-11: run.sh
// `partial` with and without BLAST_USAGE_REPORT=false, x_version, x_help, x_dontuse, x_inidont,
// x_ini, x_iniUR, x_inibothbad): a `.ncbirc` that NCBI reports is posted by each of its two reads
// (the usage report's, then the configuration's), the entries before the bad line are used, and
// the search goes on; `-version` comes after the first read, `-h` after both. A malformed
// `blastn.ini` ends the program with exit status 2: its report, then every message collected
// (the usage report's and the report again).
#[cfg(unix)]
#[test]
fn a_malformed_registry_file_is_reported_as_ncbi_reports_it() {
    let dir = scratch("syntax");
    let databases = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search")
        .display()
        .to_string();
    let query = input("arch_q.fna");
    let args = [
        "blastn",
        "-query",
        &query,
        "-db",
        "arch_snucl",
        "-outfmt",
        "6",
    ];
    let losat = Path::new(BIN);
    let reference = run_in(losat, &args, &dir, &[("BLASTDB", &databases)], b"");
    assert_eq!(reference.status.code(), Some(0), "{reference:?}");
    assert!(!reference.stdout.is_empty() && reference.stderr.is_empty());
    let ncbirc = dir.join(NCBIRC);
    fs::write(
        &ncbirc,
        format!("[BLAST]\nBLASTDB={databases}\ngarbage\nBLASTDB=/nowhere\n"),
    )
    .unwrap();
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    let both = critical(3).repeat(2);
    let partial = run_in(losat, &args, &dir, &[], b"");
    assert_eq!(text(&partial.stderr), both);
    assert_eq!(partial.status.code(), Some(0));
    assert_eq!(partial.stdout, reference.stdout);
    let late = run_in(losat, &args, &dir, &[("BLAST_USAGE_REPORT", "false")], b"");
    assert_eq!(text(&late.stderr), critical(3));
    assert_eq!(late.stdout, reference.stdout);
    let off = run_in(
        losat,
        &args,
        &dir,
        &[("NCBI_DONT_USE_NCBIRC", "1"), ("BLASTDB", &databases)],
        b"",
    );
    assert!(off.stderr.is_empty(), "{off:?}");
    assert_eq!(off.stdout, reference.stdout);
    let version = run_in(losat, &["blastn", "-version"], &dir, &[], b"");
    assert_eq!(text(&version.stderr), critical(3));
    assert!(text(&version.stdout).starts_with("blastn: 2.17.0+\n"));
    assert_eq!(version.status.code(), Some(0));
    let help = run_in(losat, &["blastn", "-h"], &dir, &[], b"");
    assert_eq!(text(&help.stderr), both);
    assert!(!help.stdout.is_empty());
    // `[NCBI] DONT_USE_NCBIRC` in blastn.ini: only the usage report reads `.ncbirc`.
    let ini = dir.join("blastn.ini");
    fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
    let ini_off = run_in(losat, &args, &dir, &[("BLASTDB", &databases)], b"");
    assert_eq!(text(&ini_off.stderr), critical(3));
    assert_eq!(ini_off.stdout, reference.stdout);
    // A malformed blastn.ini, also for `-h`; without the usage report's read, its report twice.
    fs::write(&ini, "garbage\n").unwrap();
    let report = format!("{}{INI_REPORT}", critical(3));
    for words in [&args[..], &["blastn", "-h"]] {
        let failed = run_in(losat, words, &dir, &[], b"");
        assert_eq!(failed.status.code(), Some(2), "{failed:?}");
        assert!(failed.stdout.is_empty());
        assert_eq!(text(&failed.stderr), format!("{report}{report}"));
    }
    let failed = run_in(losat, &args, &dir, &[("BLAST_USAGE_REPORT", "false")], b"");
    assert_eq!(failed.status.code(), Some(2));
    assert_eq!(text(&failed.stderr), INI_REPORT.repeat(2));
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-871 (`-logfile` is read and the
// diagnostics set up before `-version`) and c++/src/corelib/ncbidiag.cpp:3863-3874 (the messages
// collected so far go to a new log). NCBI BLAST+ 2.17.0, the official Linux x64 release (env
// parity g_reg probes of 2026-11-11: o1, o2, l, work_o/log.txt, h_log, ver_log, ver_log_cmd,
// inibad_null): with a `-logfile` file, stderr gets the usage report's post and the file gets it
// with the configuration's (and, when the program's `.ini` is malformed, the report and every
// message again); `/dev/null` gets nothing; `-version` makes the file with the usage report's
// post.
#[cfg(unix)]
#[test]
fn the_registry_messages_go_to_the_log_of_blastdbcmd_and_makeblastdb() {
    let dir = scratch("logfile");
    let databases = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search")
        .display()
        .to_string();
    let losat = Path::new(BIN);
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    let log = dir.join("log.txt");
    let logged = || {
        let content = fs::read_to_string(&log).unwrap_or_default();
        let _ = fs::remove_file(&log);
        content
    };
    fs::write(
        dir.join(NCBIRC),
        format!("[BLAST]\nBLASTDB={databases}\ngarbage\n"),
    )
    .unwrap();
    let info = ["blastdbcmd", "-info", "-db", "arch_snucl"];
    let plain = run_in(losat, &info, &dir, &[], b"");
    assert_eq!(text(&plain.stderr), critical(3).repeat(2));
    assert_eq!(plain.status.code(), Some(0));
    let mut with_log = info.to_vec();
    with_log.extend(["-logfile", "log.txt"]);
    let to_file = run_in(losat, &with_log, &dir, &[], b"");
    assert_eq!(text(&to_file.stderr), critical(3));
    assert_eq!(to_file.stdout, plain.stdout);
    assert_eq!(logged(), critical(3).repeat(2));
    let mut discarded = info.to_vec();
    discarded.extend(["-logfile", "/dev/null"]);
    let nowhere = run_in(losat, &discarded, &dir, &[], b"");
    assert_eq!(text(&nowhere.stderr), critical(3));
    assert_eq!(nowhere.stdout, plain.stdout);
    let version = run_in(
        losat,
        &["blastdbcmd", "-version", "-logfile", "log.txt"],
        &dir,
        &[],
        b"",
    );
    assert_eq!(text(&version.stderr), critical(3));
    assert_eq!(logged(), critical(3));
    // A malformed blastdbcmd.ini.
    let ini = dir.join("blastdbcmd.ini");
    fs::write(&ini, "[x\n").unwrap();
    let section = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";
    let failed = run_in(losat, &with_log, &dir, &[], b"");
    assert_eq!(failed.status.code(), Some(2), "{failed:?}");
    assert!(failed.stdout.is_empty());
    assert_eq!(text(&failed.stderr), critical(3));
    let collected = format!("{}{section}", critical(3));
    assert_eq!(logged(), collected.repeat(2));
    let failed = run_in(losat, &discarded, &dir, &[], b"");
    assert_eq!(failed.status.code(), Some(2));
    assert_eq!(text(&failed.stderr), critical(3));
    fs::remove_file(&ini).unwrap();
    // makeblastdb: the log of `-logfile` gets the messages first, then the build's lines.
    fs::write(dir.join(NCBIRC), "garbage\n").unwrap();
    let fasta = input("arch_s.fna");
    let build = run_in(
        losat,
        &[
            "makeblastdb",
            "-in",
            &fasta,
            "-dbtype",
            "nucl",
            "-out",
            "built",
            "-logfile",
            "log.txt",
        ],
        &dir,
        &[],
        b"",
    );
    assert_eq!(build.status.code(), Some(0), "{build:?}");
    assert_eq!(text(&build.stderr), critical(1));
    assert!(build.stdout.is_empty());
    let content = logged();
    let head = critical(1).repeat(2);
    assert!(
        content.starts_with(&format!("{head}\n\nBuilding a new DB, current time: ")),
        "{content}"
    );
    let version = run_in(
        losat,
        &["makeblastdb", "-version", "-logfile", "log.txt"],
        &dir,
        &[],
        b"",
    );
    assert_eq!(text(&version.stderr), critical(1));
    assert!(text(&version.stdout).starts_with("makeblastdb: 2.17.0+\n"));
    assert_eq!(logged(), critical(1));
    fs::remove_file(dir.join(NCBIRC)).unwrap();
    // Without a message, `-version` still makes the file.
    let version = run_in(
        losat,
        &["makeblastdb", "-version", "-logfile", "log.txt"],
        &dir,
        &[],
        b"",
    );
    assert!(version.stderr.is_empty());
    assert!(log.is_file());
    assert_eq!(logged(), "");
    // A search program answers `-version` after the file is made too (NCBI 2.17.0 `blastn
    // -version -logfile log.txt`: the file is made, with the usage report's message when
    // `.ncbirc` is malformed); LOSAT refuses `-logfile` only when it reads the arguments.
    // (`logged` has removed the file.)
    assert!(!log.exists());
    fs::write(dir.join(NCBIRC), "garbage\n").unwrap();
    let version = run_in(
        losat,
        &["blastn", "-version", "-logfile", "log.txt"],
        &dir,
        &[],
        b"",
    );
    assert_eq!(text(&version.stderr), critical(1));
    assert!(text(&version.stdout).starts_with("blastn: 2.17.0+\n"));
    assert_eq!(logged(), critical(1));
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1181-1187 (quoted at
// `ncbi_environment::replays_to_stderr`). NCBI BLAST+ 2.17.0, the official Linux x64 release
// (env parity review G, F1; g_reg probe_h of 2026-11-12): with a malformed `<program>.ini` and a
// `-logfile` target, blastn, blastp, blastx, tblastn and psiblast post the report to the target
// and give the collected messages to stderr when they end; tblastx (like rpsblast, makeblastdb and
// blastdbcmd) gives them to the target.
#[cfg(unix)]
#[test]
fn a_malformed_ini_with_a_logfile_is_replayed_where_ncbi_replays_it() {
    let dir = scratch("replay");
    let losat = Path::new(BIN);
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    let log = dir.join("log.txt");
    let logged = || {
        let content = fs::read_to_string(&log).unwrap_or_default();
        let _ = fs::remove_file(&log);
        content
    };
    let section = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";
    let u = critical(1);
    fs::write(dir.join(NCBIRC), "garbage\n").unwrap();
    for program in [
        "blastn", "blastp", "blastx", "tblastn", "psiblast", "tblastx",
    ] {
        let ini = dir.join(format!("{program}.ini"));
        fs::write(&ini, "[x\n").unwrap();
        let to_stderr = program != "tblastx";
        let run = |words: &[&str], envs: &[(&str, &str)]| {
            let mut args = vec![program];
            args.extend_from_slice(words);
            run_in(losat, &args, &dir, envs, b"")
        };
        let file = run(&["-query", "q.fa", "-logfile", "log.txt"], &[]);
        assert_eq!(file.status.code(), Some(2), "{program}: {file:?}");
        assert!(file.stdout.is_empty());
        if to_stderr {
            assert_eq!(text(&file.stderr), format!("{u}{u}{section}"), "{program}");
            assert_eq!(logged(), format!("{u}{section}"), "{program}");
        } else {
            assert_eq!(text(&file.stderr), u, "{program}");
            assert_eq!(logged(), format!("{u}{section}").repeat(2), "{program}");
        }
        let nowhere = run(&["-query", "q.fa", "-logfile", "/dev/null"], &[]);
        assert_eq!(nowhere.status.code(), Some(2));
        let expected = if to_stderr {
            format!("{u}{u}{section}")
        } else {
            u.clone()
        };
        assert_eq!(text(&nowhere.stderr), expected, "{program}");
        let quiet = run(
            &["-query", "q.fa", "-logfile", "log.txt"],
            &[("BLAST_USAGE_REPORT", "0")],
        );
        assert_eq!(quiet.status.code(), Some(2));
        if to_stderr {
            assert_eq!(text(&quiet.stderr), section, "{program}");
            assert_eq!(logged(), section, "{program}");
        } else {
            assert!(quiet.stderr.is_empty(), "{program}");
            assert_eq!(logged(), section.repeat(2), "{program}");
        }
        // Without `-logfile` every program writes the same.
        let plain = run(&["-query", "q.fa"], &[]);
        assert_eq!(
            text(&plain.stderr),
            format!("{u}{section}").repeat(2),
            "{program}"
        );
        fs::remove_file(&ini).unwrap();
    }
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1022 (quoted at
// `ncbi_environment::preparse` and `cli::ncbi_toolkit_word_error`) and
// c++/src/algo/blast/api/blast_usage_report.cpp:195-206 (quoted at
// `ncbi_environment::usage_report_reads_ncbirc`). NCBI BLAST+ 2.17.0 (env parity review G, F2-F4):
// `-conffile other.ini` loads that file in place of a malformed `blastn.ini` (LOSAT refuses it
// before the configuration load); `BLAST_USAGE_REPORT=abc` aborts before `-version` (LOSAT
// refuses there); `-logfile -version` logs to the file `-version` and prints no version.
#[cfg(unix)]
#[test]
fn the_words_before_the_arguments_are_read_as_ncbi_reads_them() {
    let dir = scratch("prescan");
    let losat = Path::new(BIN);
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    fs::write(dir.join("blastn.ini"), "[x\n").unwrap();
    fs::write(dir.join("other.ini"), "[BLAST]\n").unwrap();
    for conffile in [&["-conffile", "other.ini"][..], &["-conffile=other.ini"]] {
        let mut args = vec!["blastn"];
        args.extend_from_slice(conffile);
        args.extend(["-query", "q.fa"]);
        let refused = run_in(losat, &args, &dir, &[], b"");
        let stderr = text(&refused.stderr);
        assert_ne!(refused.status.code(), Some(0));
        assert!(
            stderr.contains(
                "the NCBI C++ Toolkit option -conffile is not supported by LOSAT's BLASTN"
            ) && !stderr.contains("Application's initialization failed"),
            "{stderr}"
        );
    }
    fs::remove_file(dir.join("blastn.ini")).unwrap();
    for words in [&["blastn", "-version"][..], &["blastn", "-h"]] {
        let aborts = run_in(losat, words, &dir, &[("BLAST_USAGE_REPORT", "abc")], b"");
        assert_eq!(aborts.status.code(), Some(1), "{aborts:?}");
        assert!(aborts.stdout.is_empty());
        assert!(text(&aborts.stderr)
            .contains("the environment variable BLAST_USAGE_REPORT has the value \"abc\""));
    }
    fs::write(dir.join(NCBIRC), "garbage\n").unwrap();
    let logged = run_in(losat, &["blastn", "-logfile", "-version"], &dir, &[], b"");
    assert!(logged.stdout.is_empty(), "{logged:?}");
    assert_ne!(logged.status.code(), Some(0));
    assert_eq!(text(&logged.stderr[..critical(1).len()]), critical(1));
    assert_eq!(
        fs::read_to_string(dir.join("-version")).unwrap(),
        critical(1).repeat(2)
    );
    let _ = fs::remove_dir_all(&dir);
}
