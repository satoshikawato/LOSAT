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
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// One `losat psiblast -db tests/db_search/paja -query tests/fasta/psiblast/<query>
/// -save_pssm_after_last_round <args>` from `LOSAT/`, with `HOME` the empty directory `home` and
/// without the variables that change a search (`LOSAT_BLASTX_WORKER_LOG` stays: an
/// `--all-features` binary writes its probe there).
fn psiblast(query: &str, args: &[&str], home: &Path) -> Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("psiblast")
        .args(["-db", "tests/db_search/paja"])
        .args(["-query", &format!("tests/fasta/psiblast/{query}")])
        .args(["-num_threads", "1", "-save_pssm_after_last_round"])
        .args(args);
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
    command.output().expect("run LOSAT CLI")
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
