//! TBLASTX's report stage on the search's pool (`LOSAT_X_TBXRPTPAR`,
//! `src/algorithm/tblastx/x_rpt_par.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_TBXRPTPAR` is
//! `0`, or `LOSAT_X_DEFAULT=0` and it is not set. With it on and two or more threads, the displayed
//! identities, the outfmt 0 rows and the tabular rows of the final list are computed by the pool's
//! threads; `LOSAT_X_TBXRPTPARSHADOW=1` also runs the serial stage and compares the two.
//!
//! The property under test is that the reports, the messages and the exit status do not depend on
//! the switch or on the thread count: the searches are compared with the same LOSAT search with the
//! switch off at one thread. The first query below has more than 2,048 rows (AP027131 x AP027133,
//! 120 kb each, genetic code 4 for both, which a local subject honors: the approved deviation of
//! `docs/release/tblastx_v0.1.0_certification.md`), so the identities, the outfmt 0 rows and the
//! tabular rows are cut into several chunks on the pool.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const OFF: (&str, &str) = ("LOSAT_X_TBXRPTPAR", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBXRPTPARSHADOW", "1");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const THREADS: [u32; 5] = [1, 2, 4, 8, 16];

fn letters(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    std::fs::read(path)
        .expect("test FASTA")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// A FASTA file of `records` (title, letters) in this test's temporary directory.
fn write_fasta(file: &str, records: &[(&str, &[u8])]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-tbxrptpar-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    let path = dir.join(file);
    let mut text = Vec::new();
    for (title, seq) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in seq.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).expect("write FASTA");
    path
}

/// Two queries from AP027131 (the second with a title of several words) and a record without
/// letters between them, which has a warning before its report; written once.
fn queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let a = letters("AP027131.fasta");
        write_fasta(
            "queries.fa",
            &[
                ("a", &a[..120_000]),
                ("empty record", b""),
                ("a second part", &a[400_000..430_000]),
            ],
        )
    })
}

/// A subject from AP027133 (a title of several words), written once (the tests run in parallel).
fn subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let b = letters("AP027133.fasta");
        write_fasta("subjects.fa", &[("b first part", &b[..120_000])])
    })
}

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path
/// (`LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature needs it
/// as soon as a search has more than one thread, and it selects nothing).
fn tblastx(threads: u32, outfmt: &str, switches: &[(&str, &str)]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(switches.iter().copied());
    command
        .arg("tblastx")
        .arg("-query")
        .arg(queries())
        .arg("-subject")
        .arg(subjects())
        .args(["-query_gencode", "4", "-db_gencode", "4"])
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The standard error without the shadow's lines.
fn messages(output: &Output) -> String {
    stderr(output)
        .lines()
        .filter(|line| !line.starts_with("[X_SHADOW]"))
        .map(|line| format!("{line}\n"))
        .collect()
}

/// The counters of the `[X_SHADOW] TBXRPTPAR` line: HSPs, outfmt 0 rows, tabular bytes.
fn shadow_counts(output: &Output) -> [u64; 3] {
    let text = stderr(output);
    let line = text
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] TBXRPTPAR "))
        .unwrap_or_else(|| panic!("no TBXRPTPAR shadow line:\n{text}"));
    assert!(line.ends_with(" identical"), "{line}");
    ["hsps", "pairwise", "bytes"].map(|name| {
        let key = format!("{name}=");
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .expect("a count")
    })
}

fn check_format(outfmt: &str) {
    let off = tblastx(1, outfmt, &[OFF]);
    assert!(off.status.success(), "{}", stderr(&off));
    assert!(!off.stdout.is_empty(), "the pairs have hits");
    assert!(
        stderr(&off).contains("empty record"),
        "the empty record has a warning: {}",
        stderr(&off)
    );
    if outfmt != "0" {
        // The first query's rows are formatted in windows of chunks (more than 1,024 rows).
        let rows = String::from_utf8_lossy(&off.stdout)
            .lines()
            .filter(|line| line.starts_with("a\t"))
            .count();
        assert!(
            rows > 2_048,
            "outfmt {outfmt}: {rows} rows of the first query"
        );
    }
    for threads in THREADS {
        let on = tblastx(threads, outfmt, &[]);
        assert!(on.status.success(), "threads {threads}: {}", stderr(&on));
        assert_eq!(on.stdout, off.stdout, "outfmt {outfmt}, threads {threads}");
        assert_eq!(
            messages(&on),
            messages(&off),
            "outfmt {outfmt}, threads {threads}"
        );
        // The shadow at one thread (where the stage is serial) and at 2 and 16.
        if ![1, 2, 16].contains(&threads) {
            continue;
        }
        let shadow = tblastx(threads, outfmt, &[SHADOW]);
        assert!(
            shadow.status.success(),
            "outfmt {outfmt}, threads {threads}: {}",
            stderr(&shadow)
        );
        assert_eq!(
            shadow.stdout, off.stdout,
            "outfmt {outfmt}, threads {threads}, shadow"
        );
        assert_eq!(
            messages(&shadow),
            messages(&off),
            "outfmt {outfmt}, threads {threads}"
        );
        if threads == 1 {
            // One thread: the serial stage, as before.
            assert!(
                !stderr(&shadow).contains("TBXRPTPAR"),
                "{}",
                stderr(&shadow)
            );
            continue;
        }
        let [hsps, pairwise, bytes] = shadow_counts(&shadow);
        assert!(
            hsps > 2_048,
            "outfmt {outfmt}, threads {threads}: {hsps} HSPs"
        );
        if outfmt == "0" {
            assert_eq!((pairwise, bytes), (hsps, 0), "outfmt {outfmt}");
        } else {
            assert_eq!(pairwise, 0, "outfmt {outfmt}");
            assert!(bytes > 0, "outfmt {outfmt}");
        }
    }
}

#[test]
fn the_tabular_report_does_not_depend_on_the_switch_at_any_thread_count() {
    check_format("6");
}

#[test]
fn the_commented_tabular_report_does_not_depend_on_the_switch() {
    check_format("7");
}

#[test]
fn the_pairwise_report_does_not_depend_on_the_switch() {
    check_format("0");
}

#[test]
fn the_previous_defaults_write_the_same_report() {
    let off = tblastx(4, "6", &[OFF]);
    let previous = tblastx(4, "6", &[PREVIOUS]);
    assert!(previous.status.success(), "{}", stderr(&previous));
    assert_eq!(previous.stdout, off.stdout);
}
