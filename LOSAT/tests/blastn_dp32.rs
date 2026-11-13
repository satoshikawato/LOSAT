//! The 32-bit lanes of the vector X-drop DP (`LOSAT_X_DP32`, `src/utils/xdrop_simd.rs`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! A call whose X-drop the 16-bit lanes cannot hold runs, with the switch, with the X-drop clamped
//! to the value past which NCBI's loops drop no cell and, where that is still past the 16-bit
//! lanes, on 32-bit lanes; with `=0` it runs on the scalar ports of `s_BlastAlignPackedNucl` and
//! `ALIGN_EX` as before. Both give the same result, so the property under test is that the report
//! does not depend on the switch: the default run and `LOSAT_X_DP32=0` give the same bytes and
//! print nothing on standard error, and the default run takes the switch's path (`LOSAT_X_REACH`:
//! an extension whose X-drop is past the 16-bit lanes).
//! That the extensions agree cell by cell is the business of the unit tests (`wide_lanes_*` in
//! `xdrop_simd.rs`, `wide_lanes_give_the_scalar_ports_results` in the blastn and blastp modules);
//! that the bytes equal NCBI's, of the fixture suites.
//!
//! The search: 6,000 letters of LvMJNV against a copy with a substitution every 13 letters, a
//! deletion and an insertion, `-task blastn` with reward 5, penalty -4 and gap costs 10 / 6 (a
//! supported set of NCBI's), and preliminary and final X-drops of 1e5 bits (raw about 425,000):
//! the extensions from a seed reach thousands of letters each way, so their clamped X-drop stays
//! past the 16-bit lanes and both the packed score-only extension and the traceback take the
//! 32-bit lanes.

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
        let dir = std::env::temp_dir().join(format!("losat-dp32-{test}-{}", std::process::id()));
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

fn letters(path: &Path) -> String {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// The query (letters 3,000..9,000 of LvMJNV) and the subject (the query with every 13th letter
/// changed, letters 2,000..2,003 deleted and two letters inserted at 4,000).
fn pair(dir: &TestDir) -> (PathBuf, PathBuf) {
    let genome = letters(&fasta("LvMJNV.fasta"));
    let query = &genome.as_bytes()[3_000..9_000];
    let mut subject: Vec<u8> = Vec::with_capacity(query.len());
    for (i, &c) in query.iter().enumerate() {
        if (2_000..2_003).contains(&i) {
            continue;
        }
        if i == 4_000 {
            subject.extend_from_slice(b"GA");
        }
        subject.push(if i % 13 == 6 {
            match c.to_ascii_uppercase() {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            }
        } else {
            c
        });
    }
    let q = dir.0.join("q.fna");
    std::fs::write(
        &q,
        format!(">q\n{}\n", String::from_utf8(query.to_vec()).unwrap()),
    )
    .expect("write the query");
    let s = dir.0.join("s.fna");
    std::fs::write(&s, format!(">s\n{}\n", String::from_utf8(subject).unwrap()))
        .expect("write the subject");
    (q, s)
}

/// One search (`-outfmt 0` on standard output) with the given variables as its only `LOSAT_*`
/// variables that choose a code path. Returns the report and the standard error.
fn search(query: &Path, subject: &Path, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-task", "blastn", "-reward", "5", "-penalty", "-4"])
        .args(["-gapopen", "10", "-gapextend", "6"])
        .args(["-xdrop_gap", "1e5", "-xdrop_gap_final", "1e5"])
        .args(["-outfmt", "0", "-num_threads", "1"])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{switches:?}: {stderr}");
    (output.stdout, stderr)
}

// The default search takes the switch's path (here the 32-bit lanes) and gives the bytes of
// `LOSAT_X_DP32=0`, which runs the scalar ports; neither prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_the_32_bit_lanes() {
    let dir = TestDir::new("report");
    let (query, subject) = pair(&dir);
    let reach = reach_support::Reach::new("dp32");
    let (default, default_stderr) = search(&query, &subject, &[reach.var()]);
    assert!(
        String::from_utf8_lossy(&default).contains("Identities"),
        "no alignment"
    );
    assert_eq!(default_stderr, "");
    if cfg!(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    )) {
        assert!(
            reach.sides("LOSAT_X_DP32").0,
            "no extension took the switch's path"
        );
    }
    let (report, stderr) = search(&query, &subject, &[SWITCH_OFF]);
    assert_eq!(report, default, "{SWITCH_OFF:?}");
    assert_eq!(stderr, "", "{SWITCH_OFF:?}");
}
