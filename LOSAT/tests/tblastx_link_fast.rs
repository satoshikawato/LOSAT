//! `LOSAT_LINK_FAST`: the index-backed TBLASTX linking kernel on sequences that
//! reach each of its branches, against the output of NCBI BLAST+ 2.17.0.
//!
//! The kernel reports what it did with `LOSAT_LINK_STATS=1`; every test checks
//! that the branch it is about was taken, not only that the output is equal.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

/// One TBLASTX run (`-outfmt 6` on standard output) with only the given
/// `LOSAT_*` variables set. Returns the report and the standard error.
fn tblastx(query: &Path, subject: &Path, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6"])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{switches:?}: {stderr}");
    (output.stdout, stderr)
}

/// Sum of one counter over the `[LINK_FAST_STATS]` lines (one per group).
fn counter(stderr: &str, name: &str) -> u64 {
    let key = format!("{name}=");
    let mut groups = 0;
    let total = stderr
        .lines()
        .filter(|line| line.starts_with("[LINK_FAST_STATS]"))
        .map(|line| {
            groups += 1;
            line.split_whitespace()
                .find_map(|field| field.strip_prefix(key.as_str()))
                .unwrap_or_else(|| panic!("no {name} in {line}"))
                .parse::<u64>()
                .expect("counter")
        })
        .sum();
    assert!(groups > 0, "no [LINK_FAST_STATS] line:\n{stderr}");
    total
}

/// Every check of the kernel on: each choice against the plain scan, each
/// group against the NCBI kernel.
const CHECKED: [(&str, &str); 4] = [
    ("LOSAT_LINK_FAST", "1"),
    ("LOSAT_LINK_FAST_VERIFY", "1"),
    ("LOSAT_LINK_FAST_SHADOW", "1"),
    ("LOSAT_LINK_STATS", "1"),
];

// A 4-residue HSP whose large-gap search finds, in the tree, a long HSP that
// follows it in list order (`fixtures/tblastx_link_fast/make_short_hsp.py`).
// NCBI does not link the two.
#[test]
fn a_short_hsp_followed_by_its_tree_answer_matches_ncbi() {
    let query = fixture("short_hsp.query.fna");
    let subject = fixture("short_hsp.subject.fna");
    let expected = std::fs::read(fixture("short_hsp.ncbi_2.17.0.fmt6.out")).expect("NCBI output");
    let report = String::from_utf8(expected.clone()).expect("outfmt 6");
    // The two HSPs of the construction, each with an E-value of its own.
    assert!(
        report.contains("\t100.000\t4\t0\t0\t313\t324\t313\t324\t0.065\t20.8\n"),
        "{report}"
    );
    assert!(
        report.contains("\t26\t1\t0\t311\t388\t508\t585\t1.83e-14\t62.5\n"),
        "{report}"
    );

    let (ncbi_kernel, _) = tblastx(&query, &subject, &[]);
    assert_eq!(ncbi_kernel, expected, "NCBI kernel");

    let (index_kernel, stderr) = tblastx(&query, &subject, &CHECKED);
    assert_eq!(index_kernel, expected, "index-backed kernel");
    assert!(
        counter(&stderr, "fallbacks") >= 1,
        "the scan fallback was not reached:\n{stderr}"
    );
    assert!(counter(&stderr, "verified") >= 1, "{stderr}");

    let (searching, _) = tblastx(
        &query,
        &subject,
        &[("LOSAT_LINK_FAST", "1"), ("LOSAT_LINK_FAST_REUSE0", "0")],
    );
    assert_eq!(
        searching, expected,
        "index-backed kernel, index-0 reuse off"
    );
}

fn fasta_sequence(path: &Path) -> String {
    std::fs::read_to_string(path)
        .expect("FASTA")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

struct TempFasta(PathBuf);

impl TempFasta {
    fn new(stem: &str, name: &str, sequence: &str) -> Self {
        let path =
            std::env::temp_dir().join(format!("losat_link_fast_{}_{stem}.fna", std::process::id()));
        std::fs::write(&path, format!(">{name}\n{sequence}\n")).expect("write FASTA");
        Self(path)
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

// Two overlapping 9 kb windows of one genome: eight groups of up to 65 HSPs
// that are linked over several passes, so that choices are kept from one pass
// to the next under both ordering methods.
#[test]
fn choices_kept_between_passes_match_ncbi() {
    let genome =
        fasta_sequence(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/LC738884.fasta"));
    let query = TempFasta::new("kept_query", "kept_query", &genome[0..9_000]);
    let subject = TempFasta::new("kept_subject", "kept_subject", &genome[600..9_600]);
    let expected =
        std::fs::read(fixture("kept_choices.ncbi_2.17.0.fmt6.out")).expect("NCBI output");

    let (ncbi_kernel, _) = tblastx(&query.0, &subject.0, &[]);
    assert_eq!(ncbi_kernel, expected, "NCBI kernel");

    let (index_kernel, stderr) = tblastx(&query.0, &subject.0, &CHECKED);
    assert_eq!(index_kernel, expected, "index-backed kernel");
    for name in [
        "idx0_kept_link",
        "idx0_kept_none",
        "idx1_kept_link",
        "idx1_kept_none",
    ] {
        assert!(
            counter(&stderr, name) >= 1,
            "{name} was not reached:\n{stderr}"
        );
    }
    assert!(
        counter(&stderr, "recompute_rounds") > counter(&stderr, "n") / 20,
        "{stderr}"
    );
    assert!(counter(&stderr, "verified") >= 1, "{stderr}");

    // With NCBI's index-0 rule (a search in every pass) nothing is kept there.
    let (searching, stderr) = tblastx(
        &query.0,
        &subject.0,
        &[
            ("LOSAT_LINK_FAST", "1"),
            ("LOSAT_LINK_FAST_REUSE0", "0"),
            ("LOSAT_LINK_STATS", "1"),
        ],
    );
    assert_eq!(
        searching, expected,
        "index-backed kernel, index-0 reuse off"
    );
    assert_eq!(
        counter(&stderr, "idx0_kept_link") + counter(&stderr, "idx0_kept_none"),
        0
    );
    assert!(counter(&stderr, "idx1_kept_link") >= 1, "{stderr}");
}
