//! BLASTX's diagonal entries are NCBI's `DiagStruct` bit-field: a 31-bit `last_hit` (truncated on
//! write, sign-extended on read) and a 1-bit flag.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
//! ```c
//! typedef struct DiagStruct {
//!    signed int last_hit   : 31; /**< Offset of the last hit */
//!    unsigned int flag      : 1 ; /**< Reset the next extension? */
//! } DiagStruct;
//! ```
//!
//! The diagonal offset grows by the subject's length plus the window after each subject until it
//! reaches INT4_MAX / 4 (blast_extend.c:167-174). At `-window_size 536870000` the second
//! subject's offset is 2 x 536870000 + 4,910, past 2^30, so its `s + diag_offset` wraps in the
//! 31-bit field, and NCBI 2.17.0 reports none of its rows. The fixtures
//! (`tests/blastx_diag_bitfield/README.md`) are NCBI's outputs with `-subject` and with `-db`;
//! LOSAT gives the same bytes. At that window the diagonal table has 2^30 entries (4 GiB, NCBI's
//! size too), so that test is ignored by default; the same inputs at the default window run
//! always.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/blastx_diag_bitfield")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// A temporary directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-diag-bitfield-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The protein database of the fixture's subjects, made by `losat makeblastdb` in `dir`.
fn database(dir: &Path) -> PathBuf {
    let db = dir.join("s2db");
    let made = losat()
        .args(["makeblastdb", "-in"])
        .arg(fixture("subjects.faa"))
        .args(["-dbtype", "prot", "-out"])
        .arg(&db)
        .output()
        .unwrap();
    assert!(
        made.status.success(),
        "makeblastdb failed: {}",
        String::from_utf8_lossy(&made.stderr)
    );
    db
}

/// `blastx -outfmt 6` at one thread of the fixture's query against `target` (`-subject FILE`
/// or `-db NAME`) with `-window_size window`.
fn blastx(target: (&str, &Path), window: i32) -> String {
    let out = losat()
        .args(["blastx", "-query"])
        .arg(fixture("query.fna"))
        .arg(target.0)
        .arg(target.1)
        .args(["-window_size", &window.to_string()])
        .args(["-outfmt", "6", "-num_threads", "1"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "blastx failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn expected(name: &str) -> String {
    std::fs::read_to_string(fixture(name)).unwrap()
}

#[test]
fn the_fixture_at_the_default_window_matches_ncbi() {
    let dir = temp_dir("w40");
    let subjects = fixture("subjects.faa");
    assert_eq!(
        blastx(("-subject", &subjects), 40),
        expected("ncbi.subject.w40.out6")
    );
    assert_eq!(
        blastx(("-db", &database(&dir)), 40),
        expected("ncbi.db.w40.out6")
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
#[ignore = "a 4 GiB diagonal table (2^30 entries, NCBI's size too); run with --ignored"]
fn the_second_subject_wraps_in_the_31_bit_last_hit_as_in_ncbi() {
    let dir = temp_dir("w536870000");
    let subjects = fixture("subjects.faa");
    let rows = blastx(("-subject", &subjects), 536_870_000);
    assert!(!rows.contains("YBR142W"), "{rows}");
    assert_eq!(rows, expected("ncbi.subject.w536870000.out6"));
    assert_eq!(
        blastx(("-db", &database(&dir)), 536_870_000),
        expected("ncbi.db.w536870000.out6")
    );
    std::fs::remove_dir_all(dir).ok();
}
