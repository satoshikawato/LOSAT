//! TBLASTN `-ungapped` (sum statistics with the even-gap linking) of a short query against one short
//! subject (line RG item X-short-tbx-link). NCBI's expected length of the linking cutoffs is then 0
//! (the logarithm of K times 3, NCBI's average query length of `WWWW`, times 4, the subject's 12
//! bases in protein letters, over H, rounds to 0), and with one subject the database's length is
//! the subject's, so `y_variable` and `x_variable` are 0, the
//! logarithm is -inf and the big-gap cutoff is one more than `INT_MIN`; the linking's `Int4` sums
//! (each score less the cutoff) overflow
//! and NCBI BLAST+ 2.17.0 crashes in `BLAST_LinkHsps` (SIGSEGV, exit status 139, no output), as in
//! approved exception 20 of PD-LOSAT-NCBI-DEFECTS. LOSAT computes the same cutoff
//! (`core::link_hsps::calculate_link_hsp_cutoffs`) and adds the sums in 64 bits. The expected file
//! is NCBI's own output for the nearest search it completes: the subject with a second one of 12 N,
//! which has no hits and makes the database longer than the subject (a finite cutoff), and
//! `-dbsize 12` (the subject's length), frozen from the official 2.17.0 binary
//! (`tblastn -query p4.faa -subject s12_pad.fna -ungapped -outfmt 6 -seg no -evalue 1e10
//! -comp_based_stats 0 -word_size 2 -threshold 1 -dbsize 12`, from `tests/fixtures/tblastn_short_link`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1036-1072
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:868,931-932

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastn_short_link")
        .join(name)
}

/// `losat tblastn` of `p4.faa` against the one subject `s12.fna`, with no `LOSAT_*`, `NCBI*`,
/// `BATCH_SIZE` or `BLASTDB` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with the
/// `blastx-worker-probe` feature needs that file name once a search has more than one thread).
fn tblastn(threads: &str) -> Vec<u8> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(text.as_ref(), "BATCH_SIZE" | "BLASTDB")
        {
            command.env_remove(name);
        }
    }
    let worker_log = std::env::temp_dir().join(format!(
        "losat-tblastn-ungapped-short-link-{}-{threads}.log",
        std::process::id()
    ));
    let output = command
        .env("LOSAT_BLASTX_WORKER_LOG", &worker_log)
        .arg("tblastn")
        .arg("-query")
        .arg(fixture("p4.faa"))
        .arg("-subject")
        .arg(fixture("s12.fna"))
        .args([
            "-ungapped",
            "-outfmt",
            "6",
            "-seg",
            "no",
            "-evalue",
            "1e10",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-num_threads",
            threads,
        ])
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_file(&worker_log);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success() && stderr.is_empty(),
        "threads {threads}: {stderr}"
    );
    output.stdout
}

#[test]
fn a_short_query_against_one_short_subject_gives_ncbis_nearest_completed_output() {
    let expected = std::fs::read(fixture("p4_s12.ncbi_pad_dbsize12.out")).expect("expected output");
    assert!(!expected.is_empty());
    for threads in ["1", "4"] {
        assert_eq!(
            String::from_utf8_lossy(&tblastn(threads)),
            String::from_utf8_lossy(&expected),
            "threads {threads}"
        );
    }
}
