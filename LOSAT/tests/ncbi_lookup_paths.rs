//! Where NCBI BLAST+ looks for files of its own around a search (env parity S1 and S4): the
//! registry files (`.ncbirc`, `<program>.ini`) on its search path, which depends on when NCBI
//! makes the path (`ncbi_environment::registry_search_path`), and the file of a `-matrix` that
//! its tables do not have (`ApplicationSettings::matrix_file`), which NCBI 2.17.0 reads and then
//! aborts on (rc 134, `free(): invalid size`; tblastx, and blastx `-ungapped -comp_based_stats 0`)
//! when it is a matrix.

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
// NCBI BLAST+ 2.17.0 (env parity S1, a malformed `.ncbirc`, which NCBI reports and LOSAT
// refuses): blastn run by a link on `PATH` reads one in the link's directory only with
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
    let refused = |output: &Output, file: &Path| {
        output.status.code() == Some(1)
            && String::from_utf8_lossy(&output.stderr).contains(&file.display().to_string())
    };
    // The link's directory.
    let in_shim = shim.join(".ncbirc");
    fs::write(&in_shim, malformed).unwrap();
    let plain = run_in(blastn, &args, &dir, &[("PATH", &path)], b"");
    assert_eq!(plain.status.code(), Some(0), "{plain:?}");
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
    fs::remove_file(&in_shim).unwrap();
    // `.ncbirc` turned off: the usage report looks for no file, and the path is made later
    // (ncbireg.cpp:1634-1648); `blastn.ini` in the link's directory is read.
    let ini = shim.join("blastn.ini");
    fs::write(&ini, malformed).unwrap();
    let plain = run_in(blastn, &args, &dir, &[("PATH", &path)], b"");
    assert_eq!(plain.status.code(), Some(0), "{plain:?}");
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
    let _ = fs::remove_dir_all(&dir);
}
