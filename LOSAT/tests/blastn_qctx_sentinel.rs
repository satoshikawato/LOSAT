//! Whole `-task blastn` searches whose approximate ungapped extension starts left of its query
//! context, run by the binary of the test profile, which keeps debug assertions (M-32).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:286,294
//! ```c
//!     len = MIN(q_ext, s_ext) / COMPRESSION_RATIO;
//!         Uint1 q_byte = (q[-4] << 6) | (q[-3] << 4) | (q[-2] << 2) | q[-1];
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2390-2391
//! ```c
//!     ASSERT(init_hsp->ungapped_data == NULL ||
//!            init_hsp->ungapped_data->q_start >= 0);
//! ```
//! The left loop of the approximate extension runs over the concatenated query and reads the
//! sentinel between contexts as bases, so a result kept below the reduced cutoff (and then
//! rejected) may start in the previous context; NCBI checks the start against the context only
//! for saved hits. In both searches below such a result starts two letters left of the minus
//! strand; a build that asserted it before the cutoff test stopped there: `qctx.assert` at one
//! thread in `search_query_batch`, `qctx.pairpar` (a query longer than 8000) at four threads in
//! the search pool's seed stage (`x_pair_seed.rs`). Each report must equal the NCBI BLAST+
//! 2.17.0 output frozen in `tests/fixtures/blastn_regression/` (`blastn_regression_fixtures.py`).

use std::path::{Path, PathBuf};
use std::process::Command;

fn engine(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

/// One search of a `qctx_*` query against LC738868 at `threads` threads; the report must be the
/// frozen NCBI output of `case`, with nothing on standard error. `LOSAT_*` variables other than
/// `LOSAT_BLASTX_WORKER_LOG` are dropped. A binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs
/// `LOSAT_BLASTX_WORKER_LOG` as soon as a search has more than one thread; it selects nothing, so
/// a file in the temporary directory is given when it is not set.
fn check(case: &str, query: &str, threads: u32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join(format!("losat-qctx-worker-{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(engine("tests/fixtures/blastn_regression/inputs").join(query))
        .arg("-subject")
        .arg(engine("tests/fasta/LC738868.fasta"))
        .args(["-task", "blastn", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{case}: {stderr}");
    assert_eq!(stderr, "", "{case}");
    let expected =
        std::fs::read(engine("tests/fixtures/blastn_regression").join(format!("{case}.out")))
            .expect("read the frozen NCBI output");
    assert!(!expected.is_empty(), "{case}: the frozen output has hits");
    assert!(
        output.stdout == expected,
        "{case}: the report differs from the frozen NCBI output"
    );
}

#[test]
fn a_rejected_result_left_of_the_minus_strand_at_one_thread() {
    check("qctx.assert", "qctx_assert_query.fa", 1);
}

#[test]
fn a_rejected_result_left_of_the_minus_strand_in_the_search_pool() {
    check("qctx.pairpar", "qctx_pairpar_query.fa", 4);
}
