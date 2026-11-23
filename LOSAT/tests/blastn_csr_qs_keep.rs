//! The `qs` of a CSR lookup table kept from batch to batch of blastn's batch loop
//! (`LOSAT_X_CSRQSKEEP`, `src/algorithm/blastn/x_mb_lut_keep.rs` module comment) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It chooses where a CSR build's `qs` comes from (the keep's, from an earlier batch's CSR, or a
//! fresh allocation), not what the build writes there, so the property under test is that the
//! report and the standard error do not depend on it, and that the kept `qs` is taken where it has
//! room: with `LOSAT_TIMING`, `[TIMING] x_mb_lut_keep: ... qs_kept= qs_fresh=` counts the CSR
//! builds that took it and those that allocated theirs, which add up to the search's CSR builds;
//! both are 0 with `=0`. That those bytes equal NCBI's is the business of the fixture suites.
//!
//! The subject is LC738868; `BATCH_SIZE=4000` puts each query record in a batch of its own (a
//! batch holds whole records). dc-megablast: three records of two genomes each (LC738869 and
//! LC738870, 594 kb; LC738873 and LC738874, 579 kb; LC738871 and LC738875, 720 kb), each above
//! the size rule of `LOSAT_X_DCCSR` (a quarter as many positions as the 4^11 cells), so three CSR
//! builds: the first finds no `qs` kept, the second has fewer words (and positions) than the first
//! and takes its `qs`, the third has more and allocates its own. At one thread the serial build
//! (`x_dc_csr::build`, room for the words) runs, at four the pool's (`x_disc_lookup_par`, room for
//! the positions). Megablast: four pieces of 100,000 letters from the middles of LC738869,
//! LC738870, LC738873 and LC738874, whose tables of 4^11 cells are packed CSR heads
//! (`LOSAT_X_BNCSRPACK`, `x_bn_csr::x_serial_packed` at one thread, `x_mb_lookup_once::csr_on_pool`
//! at four), each with room for a quarter of the cells (the counts' bytes), so the last three take
//! the first one's `qs`.

use std::path::{Path, PathBuf};
use std::process::Command;

mod temp_support;
use temp_support::TempFile;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_CSRQSKEEP", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");
const BATCH: (&str, &str) = ("BATCH_SIZE", "4000");

const KEEP: &str = "[TIMING] x_mb_lut_keep: ";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of a one-record FASTA fixture.
fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(&format!("{name}.fasta")))
        .expect("read fixture")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// A query file of these records in a temporary directory of its own, removed with the file when
/// the value is dropped.
fn query_file(name: &str, records: &[(String, String)]) -> TempFile {
    let text: String = records
        .iter()
        .map(|(title, seq)| format!(">{title}\n{seq}\n"))
        .collect();
    let path = TempFile::new("losat_csrqskeep", &format!("{name}.fasta"));
    std::fs::write(&path, text).expect("write the query file");
    path
}

/// One search of `query` against LC738868, `-outfmt 6` on standard output, with the given
/// variables as its only `LOSAT_*` variables (and `BATCH_SIZE`). `LOSAT_BLASTX_WORKER_LOG` stays, or
/// is set to a file of this test when the environment has none: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. Returns the report and the
/// standard error.
fn search(task: &str, query: &Path, threads: u32, vars: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    let log = std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
        .is_none()
        .then(|| {
            let log = TempFile::keep_on_failure("losat_csrqskeep_worker", "worker.log");
            command.env("LOSAT_BLASTX_WORKER_LOG", &log);
            log
        });
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs([BATCH].into_iter().chain(vars.iter().copied()))
        .output()
        .expect("run LOSAT");
    drop(log);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{task} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// `name=<n>` of the line of `stderr` that starts with `prefix`.
fn field(stderr: &str, prefix: &str, name: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(prefix))
        .unwrap_or_else(|| panic!("no {prefix} line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|word| word.strip_prefix(name))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
}

/// The standard error without its `[TIMING]` lines.
fn untimed(stderr: &str) -> Vec<&str> {
    stderr
        .lines()
        .filter(|line| !line.starts_with("[TIMING]"))
        .collect()
}

/// At one thread and four, with LOSAT_TIMING: the default and `=0` give the same report and the
/// same standard error but for the `[TIMING]` lines; the default's CSR builds (`[TIMING] <csr>
/// builds=`, all of them through the keep) took the kept `qs` `qs_kept` times and allocated the
/// others; `=0` counts neither.
fn check(task: &str, query: &Path, csr: &str, builds: u64, qs_kept: u64) {
    for threads in [1, 4] {
        let what = format!("{task} threads {threads}");
        let (expected, expected_stderr) = search(task, query, threads, &[SWITCH_OFF, TIMING]);
        assert!(!expected.is_empty(), "{what}: no hits");
        let (report, stderr) = search(task, query, threads, &[TIMING]);
        assert!(report == expected, "{what}: the report differs");
        assert_eq!(untimed(&stderr), untimed(&expected_stderr), "{what}");
        assert_eq!(field(&stderr, csr, "builds="), builds, "{what}:\n{stderr}");
        assert_eq!(
            field(&expected_stderr, csr, "builds="),
            builds,
            "{what}:\n{expected_stderr}"
        );
        assert_eq!(
            field(&stderr, KEEP, "qs_kept="),
            qs_kept,
            "{what}:\n{stderr}"
        );
        assert_eq!(
            field(&stderr, KEEP, "qs_fresh="),
            builds - qs_kept,
            "{what}:\n{stderr}"
        );
        for name in ["qs_kept=", "qs_fresh="] {
            assert_eq!(
                field(&expected_stderr, KEEP, name),
                0,
                "{what} =0:\n{expected_stderr}"
            );
        }
    }
}

// dc-megablast: three CSR builds of 4^11 cells; the second takes the first one's `qs`, the third
// (more words) allocates its own.
#[test]
fn dc_megablast_csr_builds_take_the_kept_qs() {
    let joined = |a: &str, b: &str| (format!("{a}_{b}"), letters(a) + &letters(b));
    let query = query_file(
        "dc_joined",
        &[
            joined("LC738869", "LC738870"),
            joined("LC738873", "LC738874"),
            joined("LC738871", "LC738875"),
        ],
    );
    check("dc-megablast", &query, "[TIMING] disc_lookup_csr: ", 3, 1);
}

// Megablast: four packed CSR heads of 4^11 cells; the last three take the first one's `qs`.
#[test]
fn megablast_csr_builds_take_the_kept_qs() {
    let records: Vec<(String, String)> = ["LC738869", "LC738870", "LC738873", "LC738874"]
        .iter()
        .map(|name| {
            let genome = letters(name);
            let start = genome.len() / 2 - 50_000;
            (
                format!("{name}_middle"),
                genome[start..start + 100_000].to_string(),
            )
        })
        .collect();
    let query = query_file("mb_middles", &records);
    check("megablast", &query, "[TIMING] x_bn_csr: ", 4, 3);
}
