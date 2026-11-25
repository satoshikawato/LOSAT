//! TBLASTN and TBLASTX `-db_hard_mask` on database sequences whose masks overlap or nest in the
//! order of their start (series NP, M4b-3c part 2):
//!
//! 1. Masks overlapping or nesting in start order (`s2`, `nest`, `fb`) match NCBI BLAST+ 2.17.0's
//!    `-num_threads 1` output byte for byte at both 1 and 4 threads with empty stderr.
//! 2. Refusals (exit status 1, empty stdout): NCBI's endless gapped interval-tree split (`g3`
//!    tblastn default), masks stored out of start order (`uns`), and a sequence whose word hits
//!    reach entries left by sequences searched before it (`lk`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276,804-855 and
//! c++/src/algo/blast/core/aa_ungapped.c:495-614

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const OUTFMT: &str = "6 std sframe";
const DBS: &[&str] = &["s2", "nest", "g3", "uns", "lk", "fb"];

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/t_db_hard_overlap")
}

/// Runs `losat` with `LOSAT_*` and NCBI search-changing env vars stripped, preserving
/// `LOSAT_BLASTX_WORKER_LOG` for `blastx-worker-probe` builds.
fn losat(dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE"
                    | "CHUNK_SIZE"
                    | "OVERLAP_CHUNK_SIZE"
                    | "RAYON_NUM_THREADS"
                    | "BLASTDB"
                    | "ARCHIVE_FORMAT"
                    | "WINDOW_MASKER_PATH"
                    | "BLAST_USAGE_REPORT"
                    | "BL2SEQ_LEGACY"
                    | "CTOOLKIT_COMPATIBLE"
                    | "OLD_FSC"
                    | "ADAPTIVE_CBS"
                    | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .current_dir(dir)
        .env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"));
    command
}

fn build_db(dir: &Path, name: &str) {
    let fasta_name = format!("{name}.fna");
    let asnt_name = format!("{name}.asnt");
    let output = losat(dir)
        .args([
            "makeblastdb",
            "-in",
            &fasta_name,
            "-dbtype",
            "nucl",
            "-out",
            name,
            "-mask_data",
            &asnt_name,
        ])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb {name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_search(dir: &Path, program: &str, db: &str, extra: &[&str], threads: u32) -> Output {
    let query = if program == "tblastn" {
        "q.faa"
    } else {
        "qx.fna"
    };
    losat(dir)
        .arg(program)
        .args([
            "-query",
            query,
            "-db",
            db,
            "-db_hard_mask",
            "11",
            "-seg",
            "no",
        ])
        .args(extra)
        .args(["-outfmt", OUTFMT, "-num_threads", &threads.to_string()])
        .output()
        .expect("run losat search")
}

fn check_ncbi_case(
    dir: &Path,
    program: &str,
    db: &str,
    extra: &[&str],
    expected_file: &str,
    failed: &mut Vec<String>,
) {
    let expected_path = fixture_dir().join(expected_file);
    let want = std::fs::read(&expected_path)
        .unwrap_or_else(|err| panic!("{}: {err}", expected_path.display()));
    for threads in [1, 4] {
        let output = run_search(dir, program, db, extra, threads);
        if !output.status.success() || !output.stderr.is_empty() || output.stdout != want {
            failed.push(format!(
                "{expected_file} at {threads} threads: status={} stderr={} stdout={}",
                output.status,
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&output.stdout),
            ));
        }
    }
}

fn check_refusal_case(
    dir: &Path,
    program: &str,
    db: &str,
    extra: &[&str],
    want_fragment: &str,
    failed: &mut Vec<String>,
) {
    for threads in [1, 4] {
        let result = run_search(dir, program, db, extra, threads);
        let stderr = String::from_utf8_lossy(&result.stderr);
        let stdout = String::from_utf8_lossy(&result.stdout);
        if result.status.code() != Some(1)
            || !result.stdout.is_empty()
            || !stderr.contains(want_fragment)
        {
            failed.push(format!(
                "refusal {program} {db} at {threads} threads: status={} stdout={stdout} stderr={stderr}",
                result.status,
            ));
        }
    }
}

/// Searches `s2`, `nest`, and `fb` against frozen NCBI BLAST+ 2.17.0 `-num_threads 1` outputs
/// and checks refusals on `g3` (interval tree), `uns` (out of start order), and `lk` (reaches
/// earlier sequence entries), at 1 and 4 threads.
#[test]
fn overlapping_hard_masks_match_ncbi_and_refuse_defects() {
    let dir = std::env::temp_dir().join(format!("losat-t-db-hard-overlap-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let src = fixture_dir();
    for file in ["q.faa", "qx.fna"] {
        std::fs::copy(src.join(file), dir.join(file)).unwrap();
    }
    for db in DBS {
        for ext in ["fna", "asnt"] {
            let name = format!("{db}.{ext}");
            std::fs::copy(src.join(&name), dir.join(&name)).unwrap();
        }
        build_db(&dir, db);
    }

    let w0_args = &["-window_size", "0", "-evalue", "1000"];
    let ung_args = &["-ungapped", "-comp_based_stats", "0", "-evalue", "1000"];
    let xdef_args = &["-evalue", "1000"];
    let xw0_args = &["-evalue", "1000", "-window_size", "0"];

    let mut failed = Vec::new();
    let ncbi_cases: &[(&str, &str, &[&str], &str)] = &[
        ("tblastn", "s2", &[], "tblastn.s2.def.1a.6.out"),
        ("tblastn", "s2", w0_args, "tblastn.s2.w0.1a.6.out"),
        ("tblastn", "s2", ung_args, "tblastn.s2.ung.1a.6.out"),
        ("tblastx", "s2", xdef_args, "tblastx.s2.xdef.1a.6.out"),
        ("tblastx", "s2", xw0_args, "tblastx.s2.xw0.1a.6.out"),
        ("tblastn", "nest", &[], "tblastn.nest.def.1a.6.out"),
        ("tblastn", "nest", w0_args, "tblastn.nest.w0.1a.6.out"),
        ("tblastx", "nest", xw0_args, "tblastx.nest.xw0.1a.6.out"),
        ("tblastn", "fb", w0_args, "tblastn.fb.w0.1a.6.out"),
        ("tblastx", "fb", xw0_args, "tblastx.fb.xw0.1a.6.out"),
    ];
    for &(program, db, extra, expected) in ncbi_cases {
        check_ncbi_case(&dir, program, db, extra, expected, &mut failed);
    }

    const REACH_MSG: &str = "reach the table of diagonals of the sequences searched";
    let refusals: &[(&str, &str, &[&str], &str)] = &[
        ("tblastn", "g3", &[], "splits the tree's nodes without end"),
        ("tblastn", "uns", &[], "has masks out of order"),
        ("tblastx", "uns", xdef_args, "has masks out of order"),
        ("tblastn", "lk", w0_args, REACH_MSG),
        ("tblastx", "lk", xw0_args, REACH_MSG),
    ];
    for &(program, db, extra, fragment) in refusals {
        check_refusal_case(&dir, program, db, extra, fragment, &mut failed);
    }

    std::fs::remove_dir_all(&dir).unwrap();
    assert!(failed.is_empty(), "{failed:#?}");
}
