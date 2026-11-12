//! psiblast's PSSM files when a run fails (line PS, PS3, the review's B1, B2 and S1). NCBI writes a
//! query's `-out_pssm` file after its report (`SavePssmToFile`), so a run that fails before its
//! reports are out leaves every PSSM file as it was, and a file that cannot be written ends the run
//! after that query's report. LOSAT stages the files before the reports and writes them after
//! (`algorithm::psiblast::pssm::PssmFiles`). The runs need Linux's `/dev/full` (every write fails
//! with ENOSPC) and `/dev/stdout`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:319-321
//! ```c++
//!        if (m_CmdLineArgs->GetSaveLastPssm()) {
//!            SavePssmToFile(pssm, &itr);
//!        }
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;
#[cfg(target_os = "linux")]
use std::process::Output;

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// One `losat psiblast -db tests/db_search/paja -query tests/fasta/psiblast/<query>
/// -save_pssm_after_last_round <args>` from `LOSAT/`, with `HOME` the empty directory `home` and
/// without the variables that change a search (`LOSAT_BLASTX_WORKER_LOG` stays: an
/// `--all-features` binary writes its probe there).
#[cfg(target_os = "linux")]
fn psiblast(query: &str, args: &[&str], home: &Path) -> Output {
    let mut command = losat_psiblast(home);
    command
        .args(["-db", "tests/db_search/paja"])
        .args(["-query", &format!("tests/fasta/psiblast/{query}")])
        .args(["-num_threads", "1", "-save_pssm_after_last_round"])
        .args(args);
    command.output().expect("run LOSAT CLI")
}

/// `losat psiblast` from `LOSAT/`, with `HOME` the directory `home` and without the variables that
/// change a search.
fn losat_psiblast(home: &Path) -> Command {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("psiblast");
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || name == "LOSAT_TIMING"
            || name.starts_with("NCBI")
            || name.starts_with("BLAST")
            || matches!(
                name.as_str(),
                "BATCH_SIZE" | "BL2SEQ_LEGACY" | "CHUNK_SIZE" | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
    command.env("HOME", home).env("USERPROFILE", home);
    command
}

/// An empty directory of the test's own for the PSSM files, and one for HOME.
fn scratch(tag: &str) -> (PathBuf, PathBuf) {
    let base =
        std::env::temp_dir().join(format!("losat_psiblast_pssm_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    let (files, home) = (base.join("files"), base.join("home"));
    std::fs::create_dir_all(&files).unwrap();
    std::fs::create_dir_all(&home).unwrap();
    (files, home)
}

#[cfg(target_os = "linux")]
fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    names
}

fn path(path: &Path) -> &str {
    path.to_str().unwrap()
}

// NCBI 2.17.0 (the review of PS3, `pre_full`, `full2`): `-out /dev/full` fails with status 6 and
// "BLAST failed to write output" before `SavePssmToFile`; no PSSM file is made and an existing
// one keeps its bytes.
#[cfg(target_os = "linux")]
#[test]
fn a_report_write_error_leaves_the_pssm_files_as_they_were() {
    let (dir, home) = scratch("full");
    let file = dir.join("p.asn");
    let failed = |args: &[&str]| {
        let output = psiblast("two.faa", args, &home);
        assert_eq!(output.status.code(), Some(6), "{args:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("BLAST failed to write output"),
            "{args:?}"
        );
    };
    failed(&["-out_pssm", path(&file), "-out", "/dev/full"]);
    failed(&[
        "-out_pssm",
        path(&file),
        "-save_each_pssm",
        "-out_ascii_pssm",
        path(&dir.join("p.txt")),
        "-out",
        "/dev/full",
    ]);
    assert!(names(&dir).is_empty());
    let previous = b"previous file\n";
    for name in ["p.asn", "p.asn.1"] {
        std::fs::write(dir.join(name), previous).unwrap();
    }
    failed(&["-out_pssm", path(&file), "-out", "/dev/full"]);
    failed(&[
        "-out_pssm",
        path(&file),
        "-save_each_pssm",
        "-out",
        "/dev/full",
    ]);
    assert_eq!(names(&dir), ["p.asn", "p.asn.1"]);
    for name in ["p.asn", "p.asn.1"] {
        assert_eq!(std::fs::read(dir.join(name)).unwrap(), previous, "{name}");
    }
    // The same run with a report that can be written replaces them.
    let output = psiblast(
        "two.faa",
        &["-out_pssm", path(&file), "-outfmt", "6"],
        &home,
    );
    assert_eq!(output.status.code(), Some(0));
    assert!(std::fs::read(&file)
        .unwrap()
        .starts_with(b"PssmWithParameters ::= {"));
    assert_eq!(names(&dir), ["p.asn", "p.asn.1"]);
}

// NCBI 2.17.0 (the review of PS3; the rejection row db.ps.asn.rej.nodir.7 has the same stdout): an
// -out_pssm file that cannot be written ends the run after the first query's report, without the
// epilog; NCBI then writes its serial-stream error with status 0, LOSAT its refusal.
#[cfg(target_os = "linux")]
#[test]
fn an_unwritable_pssm_file_stops_after_the_querys_report() {
    let (dir, home) = scratch("devfull");
    let output = psiblast(
        "two.faa",
        &["-out_pssm", "/dev/full", "-outfmt", "7"],
        &home,
    );
    let ncbi = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/db_search/db.ps.asn.rej.nodir.7.out"),
    )
    .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(output.stdout, ncbi);
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("a -out_pssm file that cannot be written (/dev/full: "));
    assert!(names(&dir).is_empty());
}

// The review of PS3 (B2): NCBI writes a PSSM file that is also a report stream among the reports'
// bytes in an order of their buffers (`devout`, `devfile`, `outalias`), which LOSAT refuses before
// the search.
#[cfg(target_os = "linux")]
#[test]
fn a_pssm_file_that_is_a_report_stream_is_refused_before_the_search() {
    let (dir, home) = scratch("alias");
    let report = dir.join("r.txt");
    for (args, text) in [
        (
            vec!["-out_pssm", "/dev/stdout"],
            "an -out_pssm file that is also the standard output (/dev/stdout) is not supported by LOSAT's PSIBLAST",
        ),
        (
            vec!["-out_ascii_pssm", "/dev/stdout"],
            "an -out_ascii_pssm file that is also the standard output (/dev/stdout) is not supported by LOSAT's PSIBLAST",
        ),
        (
            vec!["-out_pssm", path(&report), "-out", path(&report)],
            "an -out_pssm file that is also the -out file (",
        ),
    ] {
        let output = psiblast("two.faa", &args, &home);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(text),
            "{args:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

// NCBI 2.17.0 (the review of PS4a, S3; probes s3.file.6, s3.stdin.6, s3.stdin.trunc and
// s3.stdin.empty of `PS/tmp/fix4a`): `-in_pssm -` reads the PSSM from the standard input, as
// `CArg_InputFile` reads `-` (ncbiargs.cpp:717-721); the search gives the bytes of the run that names
// the file, and NCBI's errors name the file `-`.
#[test]
fn in_pssm_dash_reads_the_standard_input() {
    let (_, home) = scratch("stdin");
    let pssm =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search/db.ps.asn.rnd.6.asn");
    let run = |in_pssm: &str, stdin: std::process::Stdio| {
        let mut command = losat_psiblast(&home);
        command
            .args(["-in_pssm", in_pssm, "-db", "tests/db_search/paja"])
            .args(["-comp_based_stats", "0", "-evalue", "1000", "-outfmt", "6"])
            .args(["-num_threads", "1"])
            .stdin(stdin);
        command.output().expect("run LOSAT CLI")
    };
    let named = run(path(&pssm), std::process::Stdio::null());
    let piped = run("-", std::fs::File::open(&pssm).unwrap().into());
    assert_eq!(named.status.code(), Some(0));
    assert!(!named.stdout.is_empty());
    assert_eq!(
        (piped.status.code(), &piped.stdout, &piped.stderr),
        (named.status.code(), &named.stdout, &named.stderr)
    );
    let truncated = std::fs::read(&pssm).unwrap();
    let half = std::env::temp_dir().join(format!(
        "losat_psiblast_pssm_stdin_cut_{}.asn",
        std::process::id()
    ));
    std::fs::write(&half, &truncated[..truncated.len() / 2]).unwrap();
    let cut = run("-", std::fs::File::open(&half).unwrap().into());
    let _ = std::fs::remove_file(&half);
    assert_eq!(cut.status.code(), Some(1));
    assert_eq!(
        cut.stderr,
        b"BLAST query/options error: Unrecognized format for PSSM in - (must be PssmWithParameters)\nPlease refer to the BLAST+ user manual.\n"
    );
    let empty = run("-", std::process::Stdio::null());
    assert_eq!(empty.status.code(), Some(1));
    assert_eq!(
        empty.stderr,
        b"BLAST query/options error: Unsupported format for PSSM\nPlease refer to the BLAST+ user manual.\n"
    );
}

// Line PS, PS5a (DECISIONS 15; NCBI 2.17.0 step-1 cases f.bad.asn, f.bad.ascii): with more rounds
// the PSSM files are written between the rounds, the first at the start of round 2, after round 1's
// report. An `-out_pssm` file that cannot be written ends the run there (NCBI: its serial-stream
// write error, status 0; LOSAT: its refusal), without round 2 and the epilog; an `-out_ascii_pssm`
// file that cannot be opened is ignored, as NCBI's `ofstream`.
#[cfg(target_os = "linux")]
#[test]
fn an_unwritable_pssm_file_stops_the_rounds_after_the_round_before() {
    let (dir, home) = scratch("rounds");
    let run = |iterations: &str, args: &[&str]| {
        let mut command = losat_psiblast(&home);
        command
            .args(["-db", "tests/db_search/paja"])
            .args(["-query", "tests/fasta/psiblast/tk.faa"])
            .args([
                "-num_threads",
                "1",
                "-num_iterations",
                iterations,
                "-outfmt",
                "6",
            ])
            .args(args);
        command.output().expect("run LOSAT CLI")
    };
    let rounds = |args: &[&str]| run("3", args);
    let one_round = run("1", &[]);
    assert_eq!(one_round.status.code(), Some(0));
    let stopped = rounds(&["-out_pssm", "/dev/full"]);
    assert_eq!(stopped.status.code(), Some(1));
    assert_eq!(stopped.stdout, one_round.stdout);
    assert!(String::from_utf8_lossy(&stopped.stderr)
        .contains("a -out_pssm file that cannot be written (/dev/full: "));
    let missing = dir.join("nodir").join("f.txt");
    let ignored = rounds(&["-out_ascii_pssm", path(&missing)]);
    assert_eq!(ignored.status.code(), Some(0));
    assert!(ignored.stdout.ends_with(b"\nSearch has CONVERGED!\n"));
    assert!(names(&dir).is_empty());
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248 (quoted at
// `algorithm::psiblast::in_pssm_files`): with `-in_pssm` the first `SavePssmToFile` writes the
// input PSSM before round 1, after the prolog (NCBI 2.17.0, `-in_pssm tk.asn -num_iterations 2
// -out_pssm nodir/f.asn -outfmt 0`: the prolog, then its serial-stream write error, status 0).
// LOSAT writes the same prolog and refuses there (DECISIONS 10 of line PS); an `-out_ascii_pssm`
// file that cannot be opened is ignored.
#[cfg(target_os = "linux")]
#[test]
fn an_unwritable_pssm_file_of_in_pssm_stops_before_round_1() {
    let (dir, home) = scratch("in_pssm_rounds");
    let run = |outfmt: &str, args: &[&str]| {
        let mut command = losat_psiblast(&home);
        command
            .args(["-db", "tests/db_search/paja"])
            .args(["-in_pssm", "tests/fixtures/db_search/db.ps.asn.rnd.6.asn"])
            .args([
                "-num_threads",
                "1",
                "-num_iterations",
                "2",
                "-outfmt",
                outfmt,
            ])
            .args(args);
        command.output().expect("run LOSAT CLI")
    };
    let full = run("0", &[]);
    assert_eq!(full.status.code(), Some(0));
    let stopped = run("0", &["-out_pssm", "/dev/full"]);
    assert_eq!(stopped.status.code(), Some(1));
    let prolog_end = full
        .stdout
        .windows(6)
        .position(|window| window == b"Query=")
        .expect("a Query= line");
    let prolog = full.stdout[..prolog_end]
        .split(|&byte| byte == b'\n')
        .take_while(|line| !line.starts_with(b"Results from round"))
        .map(|line| line.len() + 1)
        .sum::<usize>();
    assert_eq!(stopped.stdout, full.stdout[..prolog]);
    assert!(String::from_utf8_lossy(&stopped.stderr)
        .contains("a -out_pssm file that cannot be written (/dev/full: "));
    let tabular = run("6", &["-out_pssm", "/dev/full"]);
    assert_eq!(tabular.status.code(), Some(1));
    assert!(tabular.stdout.is_empty());
    let missing = dir.join("nodir").join("f.txt");
    let ignored = run("6", &["-out_ascii_pssm", path(&missing)]);
    assert_eq!(ignored.status.code(), Some(0));
    assert_eq!(ignored.stdout, run("6", &[]).stdout);
    assert!(names(&dir).is_empty());
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:123-129
// ```c++
//     if (! m_File.Open(CSeqDB_Path(m_FileName))) {
//         //m_Atlas.Unlock(locked);
//
//         string msg = string("Error: File (") + m_FileName + ") not found.";
//
//         NCBI_THROW(CSeqDBException, eFileErr, msg);
//     }
// ```
// psiblast's preliminary search of a valid query reads a database whose sequence file is
// missing: the reader's error reaches the application (no `CLocalBlast` catch), after the prolog,
// exit 2, with the file's absolute path (NCBI 2.17.0 from `LOSAT/` with this relative `-db`;
// probe `y.g.ps.tk.nopsq.0` of the PS5c review; line PS DECISIONS 31). The path depends on the
// checkout, so this is a test here rather than a frozen row (Unix paths).
#[cfg(unix)]
#[test]
fn a_missing_sequence_file_is_the_readers_database_error() {
    let (_, home) = scratch("nopsq");
    let output = losat_psiblast(&home)
        .args(["-query", "tests/fasta/psiblast/tk.faa"])
        .args(["-db", "tests/db_search/arch_sprot_nopsq"])
        .args(["-outfmt", "0", "-num_threads", "1"])
        .output()
        .expect("run LOSAT CLI");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output.status.code(), Some(2), "{stdout}");
    let psq = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/arch_sprot_nopsq.psq");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        format!(
            "BLAST Database error: Error: File ({}) not found.\n",
            psq.display()
        )
    );
    assert!(
        stdout.starts_with("PSIBLAST 2.17.0+")
            && stdout
                .ends_with("Database: arch_sprot\n           2 sequences; 279 total letters\n\n"),
        "{stdout}"
    );
}
