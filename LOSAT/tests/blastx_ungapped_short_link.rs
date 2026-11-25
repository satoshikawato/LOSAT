//! BLASTX `-ungapped` with the one-hit word finder (`-window_size 0`; sum statistics with the
//! even-gap linking) of a short query against one short subject (line RG item X-short-tbx-link).
//! NCBI's expected length of the linking cutoffs is then 0 or less, and with one subject the
//! database is no longer than the subject, so the big-gap cutoff is one more than `INT_MIN` (gdb on
//! NCBI BLAST+ 2.17.0 for `TGGTGGTG` against `WWWW`); the linking's `Int4` sums (each score less
//! the cutoff) overflow and NCBI crashes in `BLAST_LinkHsps` (SIGSEGV, exit status 139, no output),
//! as in approved exception 20 of PD-LOSAT-NCBI-DEFECTS. LOSAT computes the same cutoff
//! (`core::link_hsps::calculate_link_hsp_cutoffs`) and adds the sums in 64 bits. The expected file
//! is NCBI's own output for the nearest search it completes: the subject with a second one of 4 X
//! (the subject's length), which has no hits and makes the database longer than the subject (a
//! finite cutoff), and `-dbsize 4` (the subject's length), frozen from the official 2.17.0 binary
//! (`blastx -query q8.fna -subject p4_pad.faa -ungapped -window_size 0 -outfmt 6 -seg no -evalue
//! 1e10 -comp_based_stats 0 -word_size 2 -threshold 1 -dbsize 4`, from
//! `tests/fixtures/blastx_short_link`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1036-1072
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:868,931-932

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/blastx_short_link")
        .join(name)
}

/// `losat blastx` of `q8.fna` against the one subject `p4.faa`, with no `LOSAT_*`, `NCBI*`,
/// `BATCH_SIZE` or `BLASTDB` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with the
/// `blastx-worker-probe` feature needs that file name once a search has more than one thread).
fn blastx(threads: &str) -> Vec<u8> {
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
        "losat-blastx-ungapped-short-link-{}-{threads}.log",
        std::process::id()
    ));
    let output = command
        .env("LOSAT_BLASTX_WORKER_LOG", &worker_log)
        .arg("blastx")
        .arg("-query")
        .arg(fixture("q8.fna"))
        .arg("-subject")
        .arg(fixture("p4.faa"))
        .args([
            "-ungapped",
            "-window_size",
            "0",
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
    let expected = std::fs::read(fixture("q8_p4.ncbi_pad_dbsize4.out")).expect("expected output");
    assert!(!expected.is_empty());
    for threads in ["1", "4"] {
        assert_eq!(
            String::from_utf8_lossy(&blastx(threads)),
            String::from_utf8_lossy(&expected),
            "threads {threads}"
        );
    }
}
