//! The default TBLASTX linking kernel (the incremental kernel of
//! `sum_stats_linking/linking_incr.rs`) on sequences that reach each of its
//! branches, against the output of NCBI BLAST+ 2.17.0.
//!
//! The kernel reports what it did with `LOSAT_LINK_STATS=1`; every test checks
//! that the branch it is about was taken, not only that the output is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-861
//! ```c
//! b0 = sum <= H_hsp_sum;
//! ...
//! b1 = q_off_t <= H_query_etrim;
//! b2 = s_off_t <= H_sub_etrim;
//! ...
//! if (!(b0|b1|b2) )
//! ```
//! The expected outputs in `tests/fixtures/tblastx_link_fast` were produced by NCBI BLAST+ 2.17.0,
//! which runs s_BlastEvenGapLinkHSPs (link_hsps.c:414-1091). The tests compare with them the output
//! of the default kernel and of the literal port (`LOSAT_LINKING_LEGACY`).

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

/// Sum of one counter over the `[LINK_INCR_STATS]` lines (one per group).
fn counter(stderr: &str, name: &str) -> u64 {
    let key = format!("{name}=");
    let mut groups = 0;
    let total = stderr
        .lines()
        .filter(|line| line.starts_with("[LINK_INCR_STATS]"))
        .map(|line| {
            groups += 1;
            line.split_whitespace()
                .find_map(|field| field.strip_prefix(key.as_str()))
                .unwrap_or_else(|| panic!("no {name} in {line}"))
                .parse::<u64>()
                .expect("counter")
        })
        .sum();
    assert!(groups > 0, "no [LINK_INCR_STATS] line:\n{stderr}");
    total
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419
// ```c
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
// ```
// The checks compare the kernel with the port of s_BlastEvenGapLinkHSPs (SHADOW) and every pass
// with a pass computed from scratch by NCBI's scans (VERIFY).
/// Every check of the default kernel on: each pass against a pass from
/// scratch, each group against the literal port.
const CHECKED: [(&str, &str); 3] = [
    ("LOSAT_LINK_FAST_VERIFY", "1"),
    ("LOSAT_LINK_FAST_SHADOW", "1"),
    ("LOSAT_LINK_STATS", "1"),
];

/// The literal port of s_BlastEvenGapLinkHSPs instead of the default kernel.
const LITERAL_PORT: [(&str, &str); 1] = [("LOSAT_LINKING_LEGACY", "1")];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-861
// ```c
// for (H2_index=H_index-1; H2_index>1;)
// ...
//    b0 = sum <= H_hsp_sum;
// ...
//    b1 = q_off_t <= H_query_etrim;
//    b2 = s_off_t <= H_sub_etrim;
// ...
//    if (!(b0|b1|b2) )
// ```
// NCBI scans only the HSPs before H in list order, so it does not offer the later long HSP to H.
// The test reaches the fallback that does the same scan.
// A 4-residue HSP whose large-gap search finds, in the tree of the sweep, a
// long HSP that follows it in list order
// (`fixtures/tblastx_link_fast/make_short_hsp.py`). NCBI does not link the two.
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

    let (literal_port, _) = tblastx(&query, &subject, &LITERAL_PORT);
    assert_eq!(literal_port, expected, "literal port");

    let (default_kernel, _) = tblastx(&query, &subject, &[]);
    assert_eq!(default_kernel, expected, "default kernel");

    let (checked, stderr) = tblastx(&query, &subject, &CHECKED);
    assert_eq!(checked, expected, "default kernel, checked");
    assert!(
        counter(&stderr, "fallbacks") >= 1,
        "the scan fallback was not reached:\n{stderr}"
    );
    assert!(counter(&stderr, "verified_passes") >= 1, "{stderr}");
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:781-795
// ```c
// H2 = H->hsp_link.link[index];
// if ((!first_pass) && ((H2==0) || (H2->hsp_link.changed==0)))
// ...
//    if(H2){
//       H_hsp_num=H2->hsp_link.num[index];
//       H_hsp_sum=H2->hsp_link.sum[index];
//       H_hsp_xsum=H2->hsp_link.xsum[index];
//    }
//    H_hsp_link=H2;
//    H->hsp_link.changed=0;
// ```
// NCBI keeps a choice at index 1 when the HSP it selected kept its own; the incremental kernel
// computes again, at both indexes, only the HSPs whose choice a removal can change, and the
// other HSPs keep theirs. The test reaches those later passes, including searches that select
// the previous choice again (same0, same1); the unit tests of linking.rs reach the choices kept
// without a search (kept0, kept1).
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

    let (literal_port, _) = tblastx(&query.0, &subject.0, &LITERAL_PORT);
    assert_eq!(literal_port, expected, "literal port");

    let (default_kernel, _) = tblastx(&query.0, &subject.0, &[]);
    assert_eq!(default_kernel, expected, "default kernel");

    let (checked, stderr) = tblastx(&query.0, &subject.0, &CHECKED);
    assert_eq!(checked, expected, "default kernel, checked");
    for name in [
        "searched0",
        "searched1",
        "same0",
        "same1",
        "removed",
        "swept",
        "verified_passes",
    ] {
        assert!(
            counter(&stderr, name) >= 1,
            "{name} was not reached:\n{stderr}"
        );
    }
    assert!(
        counter(&stderr, "passes") > counter(&stderr, "n") / 20,
        "{stderr}"
    );
}
