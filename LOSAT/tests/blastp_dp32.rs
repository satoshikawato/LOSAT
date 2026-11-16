//! The 32-bit lanes of the vector X-drop DP and the X-drop clamp (`LOSAT_X_DP32`,
//! `src/utils/xdrop_simd.rs`) on whole blastp searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! An extension whose X-drop the 16-bit lanes cannot hold runs, with the switch, with the X-drop
//! clamped to the value past which NCBI's loops drop no cell, on the 16-bit or 8-bit lanes when
//! that fits them and on 32-bit lanes when it does not; with `=0` it runs on the scalar ports of
//! `Blast_SemiGappedAlign` and `ALIGN_EX` as before. Both give the same result, so the property
//! under test is that the report does not depend on the switch: the default run and
//! `LOSAT_X_DP32=0` give the same bytes and print nothing on standard error, and the default run
//! takes the switch's path (`LOSAT_X_REACH`). That the extensions agree cell by cell is the
//! business of the unit tests; that the bytes equal NCBI's, of the fixture suites.
//!
//! The searches: the proteins of the fixture `xdp.blastp.f1e8` (`outfmt0_manifest.tsv`) with a final
//! X-drop of 1e8 bits and without composition-based statistics, whose traceback then runs with
//! that X-drop (the composition-based redo has its own; every extension is clamped onto the 16-bit
//! lanes), and a protein of 5,000 random letters
//! against a copy with every ninth letter changed, without composition-based statistics and with
//! preliminary and final X-drops of 1e5 bits (extensions of about 2,500 letters each way, whose
//! clamped X-drop stays past the 16-bit lanes: the 32-bit lanes).

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DP32", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-bpdp32-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

/// The random protein and its copy with every ninth letter changed.
fn pair(dir: &TestDir) -> (PathBuf, PathBuf) {
    const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = 0x2545_F491u32;
    let mut next = move || {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        (state >> 8) as usize
    };
    let query: Vec<u8> = (0..5_000).map(|_| LETTERS[next() % 20]).collect();
    let subject: Vec<u8> = query
        .iter()
        .enumerate()
        .map(|(i, &c)| if i % 9 == 4 { LETTERS[next() % 20] } else { c })
        .collect();
    let q = dir.0.join("q.faa");
    std::fs::write(&q, format!(">q\n{}\n", String::from_utf8(query).unwrap())).expect("write");
    let s = dir.0.join("s.faa");
    std::fs::write(&s, format!(">s\n{}\n", String::from_utf8(subject).unwrap())).expect("write");
    (q, s)
}

/// One search (`-outfmt 0` on standard output) with the given variables as its only `LOSAT_*`
/// variables that choose a code path. Returns the report and the standard error.
fn search(
    query: &Path,
    subject: &Path,
    args: &[&str],
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(args)
        .args(["-outfmt", "0", "-num_threads", "1"])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{args:?} {switches:?}: {stderr}");
    (output.stdout, stderr)
}

// Each default search takes the switch's path and gives the bytes of `LOSAT_X_DP32=0`, which runs
// the scalar ports; no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_the_32_bit_lanes() {
    let dir = TestDir::new("report");
    let (query, subject) = pair(&dir);
    let fixture = (
        fasta("outfmt0/e2e_protein_query.faa"),
        fasta("outfmt0/e2e_protein_subject.faa"),
    );
    let searches: [(&Path, &Path, &[&str]); 2] = [
        (
            &fixture.0,
            &fixture.1,
            &["-comp_based_stats", "0", "-xdrop_gap_final", "1e8"],
        ),
        (
            &query,
            &subject,
            &[
                "-comp_based_stats",
                "0",
                "-xdrop_gap",
                "1e5",
                "-xdrop_gap_final",
                "1e5",
            ],
        ),
    ];
    for (q, s, args) in searches {
        let reach = reach_support::Reach::new("bpdp32");
        let (default, default_stderr) = search(q, s, args, &[reach.var()]);
        assert!(
            String::from_utf8_lossy(&default).contains("Identities"),
            "{args:?}: no alignment"
        );
        assert_eq!(default_stderr, "", "{args:?}");
        if cfg!(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128"),
            target_arch = "aarch64"
        )) {
            assert!(
                reach.sides("LOSAT_X_DP32").0,
                "{args:?}: no extension took the switch's path"
            );
        }
        let (report, stderr) = search(q, s, args, &[SWITCH_OFF]);
        assert_eq!(report, default, "{args:?} {SWITCH_OFF:?}");
        assert_eq!(stderr, "", "{args:?} {SWITCH_OFF:?}");
    }
}
