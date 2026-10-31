//! `LOSAT_X_TBXBAND` on whole runs of the binary: the unit path takes a query
//! batch whose longest context and longest subject chunk are both longer than the slack of the
//! reference's diagonal table (`guard` in `src/algorithm/tblastx/x_pair_par.rs`, which the switch no
//! longer applies; module doc item 3 for why the units stay exact), when the subject chunks have at
//! most as many letters as that table has cells (`x_band_fits`, a memory rule), or, with
//! `LOSAT_X_TBXUNITMASK` (also in the strict set), when the units whose chunks are longer take
//! masked tables. The query is the first 1,900 bases of MjeNMV (six frames of 632-633 letters, a
//! table of 4,096 cells, a slack of about 290); the subject is MjeNMV's first 9,000 bases (frames of
//! about 3,000 letters: the band) or the PemoMJNVA genome (frames of about 99,000 letters: the
//! memory rule keeps the guard without `LOSAT_X_TBXUNITMASK`, masked units with it).
//! `LOSAT_X_PAIRPARGATE` still decides after the guards: at one thread the batch (about 3,800
//! letters, under the 5,000 of `LOSAT_X_PAIRPARN1`) keeps the reference loop ("fallback: cost"), at
//! several the units split the one subject over the pool.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-63
//! ```c
//! diag_array_length = 1;
//! /* What power of 2 is just longer than the query? */
//! while (diag_array_length < (qlen+window_size))
//! {
//!         diag_array_length = diag_array_length << 1;
//! }
//! ...
//! diag_table->diag_mask = diag_array_length-1;
//! ```
//! The switch is in the strict set (on unless `LOSAT_X_TBXBAND=0`). The tests compare the reports
//! of the default and of `LOSAT_X_TBXBAND=0`, read `[X_STATS] PAIRPAR ... units= fallback_guard=
//! fallback_cost= band=`.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path, then `switches`.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    for (name, value) in switches {
        command.env(name, value);
    }
    command
}

/// A new file of one record, the first `length` bases of `name`. Every call writes its own file: the
/// tests of this binary run in parallel, and a shared name let one test read the file while another
/// was writing it again (an empty or cut query).
fn query(name: &str, length: usize) -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!("losat-tblastx-band-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let letters: Vec<u8> = std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect();
    let mut text = String::from(">q\n");
    for line in letters[..length].chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    let call = CALLS.fetch_add(1, Ordering::Relaxed);
    let file = dir.join(format!("q{length}-{call}.fna"));
    std::fs::write(&file, text).unwrap();
    file
}

/// One run: the report (`-outfmt 6`) and the whole stderr.
fn tblastx(
    query: &Path,
    subject: &Path,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat(switches)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The counters `names` of the stderr line starting with `prefix`.
fn counters<const N: usize>(stderr: &str, prefix: &str, names: [&str; N]) -> [u64; N] {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(prefix))
        .unwrap_or_else(|| panic!("no {prefix} line:\n{stderr}"));
    names.map(|name| {
        let key = format!("{name}=");
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .unwrap()
    })
}

const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const BAND_OFF: (&str, &str) = ("LOSAT_X_TBXBAND", "0");
const PAIRPAR_LINE: &str = "[X_STATS] PAIRPAR ";
const PAIRPAR_FIELDS: [&str; 4] = ["units", "fallback_guard", "fallback_cost", "band"];

// No NCBI counterpart: LOSAT_X_TBXBAND chooses between two LOSAT paths with the same result.
#[test]
fn a_band_batch_takes_the_units_at_several_threads_with_the_same_report() {
    let q = query("MjeNMV.fasta", 1_900);
    let s = query("MjeNMV.fasta", 9_000);
    for threads in [2, 4] {
        let (reference, off) = tblastx(&q, &s, threads, &[STATS, BAND_OFF]);
        assert!(
            reference.iter().filter(|&&b| b == b'\n').count() > 5,
            "the search finds hits"
        );
        assert_eq!(
            counters(&off, PAIRPAR_LINE, PAIRPAR_FIELDS),
            [0, 1, 0, 0],
            "threads {threads}: without the switch the guard sends the batch to the reference loop"
        );
        let (report, on) = tblastx(&q, &s, threads, &[STATS]);
        assert!(
            report == reference,
            "threads {threads}: the report differs with the switch"
        );
        let [units, guard, cost, band] = counters(&on, PAIRPAR_LINE, PAIRPAR_FIELDS);
        assert!(units > 0, "threads {threads}: the units run\n{on}");
        assert_eq!(
            (guard, cost, band),
            (0, 0, 1),
            "threads {threads}: the batch of the band"
        );
    }
}

// At one thread `LOSAT_X_PAIRPARGATE` keeps the reference loop for this batch (under the 5,000
// letters of `LOSAT_X_PAIRPARN1`), with the switch as without it.
#[test]
fn the_cost_gate_still_decides_at_one_thread() {
    let q = query("MjeNMV.fasta", 1_900);
    let s = query("MjeNMV.fasta", 9_000);
    let (reference, off) = tblastx(&q, &s, 1, &[STATS, BAND_OFF]);
    assert_eq!(counters(&off, PAIRPAR_LINE, PAIRPAR_FIELDS), [0, 1, 0, 0]);
    let (report, on) = tblastx(&q, &s, 1, &[STATS]);
    assert!(report == reference, "the report differs with the switch");
    assert_eq!(counters(&on, PAIRPAR_LINE, PAIRPAR_FIELDS), [0, 0, 1, 0]);
}

const MASK_OFF: (&str, &str) = ("LOSAT_X_TBXUNITMASK", "0");
/// The size rule's tables of at most 2 KiB (512 cells): below this query's group tables (1,024
/// cells or more), so only the band batch's rule masks them.
const MASK_BAND_ONLY: (&str, &str) = ("LOSAT_X_TBXUNITMASK", "2");

// A subject longer than the reference table's cells (the PemoMJNVA genome): without
// `LOSAT_X_TBXUNITMASK` the memory rule keeps the guard (the reference loop at several threads);
// with it (the default) the units take the batch, every one on a masked table, with the same
// report.
#[test]
fn a_subject_longer_than_the_table_takes_masked_units() {
    let q = query("MjeNMV.fasta", 1_900);
    let genome = fasta("PemoMJNVA.fasta");
    let (reference, off) = tblastx(&q, &genome, 2, &[STATS, BAND_OFF]);
    assert!(
        reference.iter().filter(|&&b| b == b'\n').count() > 5,
        "the search finds hits"
    );
    assert_eq!(counters(&off, PAIRPAR_LINE, PAIRPAR_FIELDS), [0, 1, 0, 0]);
    let (report, unmasked) = tblastx(&q, &genome, 2, &[STATS, MASK_OFF]);
    assert!(
        report == reference,
        "the report differs without the masked tables"
    );
    assert_eq!(
        counters(&unmasked, PAIRPAR_LINE, PAIRPAR_FIELDS),
        [0, 1, 0, 0]
    );
    for threads in [2, 4] {
        let (report, on) = tblastx(&q, &genome, threads, &[STATS]);
        assert!(
            report == reference,
            "threads {threads}: the report differs with the masked tables"
        );
        let [units, guard, cost, band] = counters(&on, PAIRPAR_LINE, PAIRPAR_FIELDS);
        assert!(units > 0, "threads {threads}: the units run\n{on}");
        assert_eq!((guard, cost, band), (0, 0, 1), "threads {threads}");
        assert_eq!(
            counters(&on, "[X_STATS] TBXUNITMASK ", ["units_masked"]),
            [units],
            "threads {threads}: every unit's chunk is longer than the reference table"
        );
        // The band batch's rule alone (chunks longer than the reference table) masks every unit.
        let (report, band_only) = tblastx(&q, &genome, threads, &[STATS, MASK_BAND_ONLY]);
        assert!(
            report == reference,
            "threads {threads}: the report differs with the band rule alone"
        );
        let [units, ..] = counters(&band_only, PAIRPAR_LINE, PAIRPAR_FIELDS);
        assert_eq!(
            counters(
                &band_only,
                "[X_STATS] TBXUNITMASK ",
                ["kib", "units_masked"]
            ),
            [2, units],
            "threads {threads}"
        );
    }
}
