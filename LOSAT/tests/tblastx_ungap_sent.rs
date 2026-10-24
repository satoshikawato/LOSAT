//! TBLASTX's two-hit ungapped extension with walks that the sentinels stop
//! (`LOSAT_X_TBXUNGAPSENT`, `src/algorithm/tblastx/extension/x_two_hit_fast.rs`, module header) on
//! whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155
//! ```c
//! static Int4 s_BlastAaExtendTwoHit(...)
//! {
//!     ...
//!     left_score = s_BlastAaExtendLeft(matrix, subject, query, ...);
//!     ...
//!         right_score = s_BlastAaExtendRight(matrix, subject, query, ...);
//! ```
//! The switch is in the strict set (on unless it is `0`). The path gives the reference
//! extension's result (the module's argument), so the property under test is that the report does
//! not depend on it: in the unit path at 1 and 4 threads (`LOSAT_X_PAIRPARGATE=0`: the gate would
//! send these short batches to the reference frame loop) and in the reference frame loop
//! (`LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1`). With `LOSAT_TBLASTX_TEST_CHUNK_SIZE` (debug builds) every
//! subject frame is searched in chunks, which stay off the path. That reports equal NCBI's is the
//! business of the TBLASTX fixtures. `LOSAT_X_REACH` shows which side each run took.
//!
//! The query is the first 12,000 bases of AP027131, the subject the first 12,000 bases of the
//! related AP027133 (genetic code 4 for both), so there are long walks and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: &str = "LOSAT_X_TBXUNGAPSENT";
const OFF: (&str, &str) = ("LOSAT_X_TBXUNGAPSENT", "0");
const REFERENCE_LOOP: (&str, &str) = ("LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_tblastx_ungap_sent_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The first `n` bases of the one record of a FASTA file, as FASTA text.
fn bases(name: &str, n: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut lines = text.lines();
    let header = lines.next().expect("a header");
    let sequence: String = lines.flat_map(str::chars).take(n).collect();
    assert_eq!(sequence.len(), n, "{name} is shorter than {n} bases");
    format!("{header}\n{sequence}\n")
}

/// The query and subject files, removed when dropped.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("query.fna");
        let subject = dir.join("subject.fna");
        std::fs::write(&query, bases("AP027131.fasta", 12_000)).unwrap();
        std::fs::write(&subject, bases("AP027133.fasta", 12_000)).unwrap();
        Inputs {
            dir,
            query,
            subject,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The binary with `env` as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` names a file of its own: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. `LOSAT_X_PAIRPARGATE=0`
/// comes first (as in `tblastx_kmerge.rs`).
fn losat(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .env("LOSAT_X_PAIRPARGATE", "0")
        .envs(env.iter().copied());
    command
}

/// The sides of `SWITCH` that a run with `LOSAT_X_REACH=dir` recorded: (new, ref).
fn reached(dir: &Path) -> (bool, bool) {
    let (mut new, mut reference) = (false, false);
    for file in std::fs::read_dir(dir).expect("the reach directory") {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for row in text.lines().skip(1) {
            let f: Vec<&str> = row.split('\t').collect();
            if f[0] == SWITCH {
                new |= f[1] == "1";
                reference |= f[2] == "1";
            }
        }
    }
    (new, reference)
}

/// `tblastx` of the inputs at `threads` with `env` and `LOSAT_X_REACH` set to a directory of its
/// own: the report and the sides recorded.
fn tblastx(inputs: &Inputs, threads: u32, env: &[(&str, &str)]) -> (Vec<u8>, (bool, bool)) {
    let dir = temp_path("reach");
    std::fs::create_dir_all(&dir).unwrap();
    let output = losat(env)
        .env("LOSAT_X_REACH", &dir)
        .arg("tblastx")
        .arg("-query")
        .arg(&inputs.query)
        .arg("-subject")
        .arg(&inputs.subject)
        .args(["-query_gencode", "4", "-db_gencode", "4", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "{env:?} threads {threads}: {stderr}"
    );
    let sides = reached(&dir);
    let _ = std::fs::remove_dir_all(&dir);
    (output.stdout, sides)
}

/// The default (the unit path at 1 and 4 threads, the reference frame loop) gives the report of
/// `=0` and records only the switch's own side; `=0` records nothing.
#[test]
fn reports_do_not_depend_on_the_sentinel_walks() {
    let inputs = Inputs::new();
    let (expected, sides) = tblastx(&inputs, 1, &[OFF]);
    assert!(!expected.is_empty(), "no hits");
    assert_eq!(sides, (false, false), "=0 records nothing");
    for (name, threads, env) in [
        ("unit path, 1 thread", 1, &[][..]),
        ("unit path, 4 threads", 4, &[][..]),
        ("reference frame loop", 1, &[REFERENCE_LOOP][..]),
    ] {
        let (report, sides) = tblastx(&inputs, threads, env);
        assert_eq!(report, expected, "{name}");
        assert_eq!(sides, (true, false), "{name}: the sentinel walks ran");
    }
    let (report, _) = tblastx(&inputs, 1, &[OFF, REFERENCE_LOOP]);
    assert_eq!(report, expected, "reference frame loop, =0");
}

/// Subject frames searched in chunks (`LOSAT_TBLASTX_TEST_CHUNK_SIZE`, read by debug builds): a
/// chunk is not the whole frame, so the extensions take the walks with a bound, in the unit path
/// and in the reference frame loop; the report is that of `=0` with the same chunks.
#[cfg(debug_assertions)]
#[test]
fn chunks_of_a_frame_stay_off_the_path() {
    let inputs = Inputs::new();
    let chunks = ("LOSAT_TBLASTX_TEST_CHUNK_SIZE", "1500");
    let (expected, _) = tblastx(&inputs, 1, &[chunks, OFF]);
    assert!(!expected.is_empty(), "no hits with chunks");
    for (name, threads, env) in [
        ("unit path, 1 thread", 1, &[chunks][..]),
        ("unit path, 4 threads", 4, &[chunks][..]),
        ("reference frame loop", 1, &[chunks, REFERENCE_LOOP][..]),
    ] {
        let (report, sides) = tblastx(&inputs, threads, env);
        assert_eq!(report, expected, "{name}");
        assert_eq!(
            sides,
            (false, true),
            "{name}: chunks take the bounded walks"
        );
    }
}
