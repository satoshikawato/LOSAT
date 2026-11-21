//! BLASTP's refusals of a step where NCBI's search of a later query batch stops (line RG item X4,
//! fix round 2; PD-LOSAT-NCBI-DEFECTS 1.53, decision D13): NCBI prints each batch's reports before
//! it searches the next, so where its search of batch 2 stops with an exception, batch 1's reports
//! are printed (without the epilog). LOSAT refuses after them. `x4/chunk_second.faa`: a protein of
//! PajaWSV (92 letters, `x4/chunk_first.faa`) and the TK-TMK protein of psiblast's rows (398
//! letters); with `BATCH_SIZE=50` each is a batch, and with `CHUNK_SIZE=150` (the overlap 100) the
//! first is not split (92 / 50 makes one chunk) and the second is split into chunks of 143 letters
//! that NCBI would split again, which NCBI stops on with a null reference (`db.ps.it.rej.split`'s
//! case). `x4/chunk_empty.faa`: 9 and 31 letters of that protein; with `BATCH_SIZE=1 CHUNK_SIZE=100
//! OVERLAP_CHUNK_SIZE=95` the first is not split (9 / 5 makes one chunk) and the second makes 6
//! chunks of an evened size of 85, below the overlap: chunk 0 covers the batch (and would be split
//! again) and chunks 1 to 5 hold no query, on which NCBI's `CQuerySplitter::Split` stops first
//! (split_query_cxx.cpp:234,876-879). A negative pair (fix round 4): `outfmt0/xsb_ql.faa` (qa, 31
//! letters, and ql, qa ten times) with `BATCH_SIZE=1 CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200` (NCBI's
//! `size_t` sizes 100 apart): qa is one chunk (31 / 100), ql makes 3 chunks whose evened size wraps
//! past its 310 letters, so chunk 0 covers ql and chunks 1 and 2 hold no query: NCBI stops there
//! after qa's report (NCBI BLAST+ 2.17.0 n1: exit status 3, rows `db.x4.chunk.neg.*`). Where the wrapped
//! ranges leave gaps between chunks that all hold a query (`x4/chunk_second.faa` with `CHUNK_SIZE=-1
//! OVERLAP_CHUNK_SIZE=-200`: the 398 letters make chunks [0, 100) and [300, 398)), NCBI searches on
//! without the hits in the gap; the maintainer kept that setting refused (PD-LOSAT-NCBI-DEFECTS
//! version 1.1), and LOSAT refuses where that batch's search would start, after batch 1's report
//! (review 3 of X4; row `db.x4.chunk.neg.gap.6`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276

use std::path::{Path, PathBuf};
use std::process::Command;

fn path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(name)
}

struct Run {
    stdout: String,
    stderr: String,
    code: Option<i32>,
}

/// `losat blastp -query QUERY -db tests/db_search/paja -outfmt OUTFMT -num_threads THREADS` with
/// `env` as the only variables that choose a code path or NCBI's batches and chunks.
fn search(query: &str, outfmt: &str, threads: usize, env: &[(&str, &str)]) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        // `LOSAT_BLASTX_WORKER_LOG` is the log file of the build with `--all-features` (its worker
        // probe), not a code path.
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE"
                    | "CHUNK_SIZE"
                    | "OVERLAP_CHUNK_SIZE"
                    | "RAYON_NUM_THREADS"
                    | "BLASTDB"
            )
        {
            command.env_remove(name);
        }
    }
    let output = command
        .envs(env.iter().copied())
        .arg("blastp")
        .arg("-query")
        .arg(path(query))
        .arg("-db")
        .arg(path("tests/db_search/paja"))
        .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

#[test]
fn a_chunk_split_again_in_batch_2_is_refused_after_batch_1() {
    let env = [("BATCH_SIZE", "50"), ("CHUNK_SIZE", "150")];
    let alone = search("tests/fasta/x4/chunk_first.faa", "6", 1, &env);
    let alone7 = search("tests/fasta/x4/chunk_first.faa", "7", 1, &env);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert!(!alone.stdout.is_empty());
    let epilog = alone7
        .stdout
        .rfind("# BLAST processed")
        .expect("outfmt 7's last line");
    for threads in [1, 4] {
        let run = search("tests/fasta/x4/chunk_second.faa", "6", threads, &env);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert!(
            run.stderr.contains("would split a query chunk again"),
            "n{threads}: {}",
            run.stderr
        );
        let commented = search("tests/fasta/x4/chunk_second.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}");
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        // In one batch (the default batch size) the refusal comes before the first byte.
        let one_batch = search(
            "tests/fasta/x4/chunk_second.faa",
            "6",
            threads,
            &[("CHUNK_SIZE", "150")],
        );
        assert_eq!(one_batch.code, Some(1), "n{threads}");
        assert_eq!(one_batch.stdout, "", "n{threads}");
    }
}

#[test]
fn a_chunk_without_a_query_in_batch_2_is_refused_before_a_split_again() {
    let env = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "100"),
        ("OVERLAP_CHUNK_SIZE", "95"),
    ];
    let alone7 = search("tests/fasta/x4/chunk_empty.faa", "7", 1, &[]);
    assert_eq!(alone7.code, Some(0), "{}", alone7.stderr);
    // Batch 1's comment lines and rows: the part of the run without the batches' epilog before
    // q31's `# Query:` line.
    let batch1 = alone7
        .stdout
        .find(concat!(
            "# BLASTP (LOSAT ",
            env!("CARGO_PKG_VERSION"),
            ")\n# Query: q31"
        ))
        .expect("q31's report");
    for threads in [1, 4] {
        let commented = search("tests/fasta/x4/chunk_empty.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}: {}", commented.stderr);
        assert!(
            commented
                .stderr
                .contains("a query chunk without a query is not supported by LOSAT's BLASTP"),
            "n{threads}: {}",
            commented.stderr
        );
        assert_eq!(commented.stdout, alone7.stdout[..batch1], "n{threads}");
    }
}

#[test]
fn a_negative_pair_that_leaves_chunks_without_a_query_in_batch_2_is_refused_after_batch_1() {
    let env = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-100"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ];
    // qa's report in one batch of both queries without the variables (not split): its rows, and
    // its comment lines and rows before ql's `# Query:` line.
    let whole = search("tests/fasta/outfmt0/xsb_ql.faa", "6", 1, &[]);
    let whole7 = search("tests/fasta/outfmt0/xsb_ql.faa", "7", 1, &[]);
    assert_eq!(whole.code, Some(0), "{}", whole.stderr);
    let qa: String = whole
        .stdout
        .lines()
        .filter(|line| line.starts_with("qa\t"))
        .map(|line| format!("{line}\n"))
        .collect();
    assert!(!qa.is_empty());
    let batch1 = whole7
        .stdout
        .find(concat!(
            "# BLASTP (LOSAT ",
            env!("CARGO_PKG_VERSION"),
            ")\n# Query: ql"
        ))
        .expect("ql's report");
    for threads in [1, 4] {
        let run = search("tests/fasta/outfmt0/xsb_ql.faa", "6", threads, &env);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert!(
            run.stderr
                .contains("a query chunk without a query is not supported by LOSAT's BLASTP"),
            "n{threads}: {}",
            run.stderr
        );
        assert_eq!(run.stdout, qa, "n{threads}");
        let commented = search("tests/fasta/outfmt0/xsb_ql.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}: {}", commented.stderr);
        assert_eq!(commented.stdout, whole7.stdout[..batch1], "n{threads}");
    }
}

#[test]
fn a_negative_pair_that_leaves_gaps_between_its_chunks_in_batch_2_is_refused_after_batch_1() {
    let env = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-1"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ];
    // Batch 1 (qsok2, 92 letters: one chunk with this pair) alone, without the variables.
    let alone = search("tests/fasta/x4/chunk_first.faa", "6", 1, &[]);
    let alone7 = search("tests/fasta/x4/chunk_first.faa", "7", 1, &[]);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert!(!alone.stdout.is_empty());
    let epilog = alone7
        .stdout
        .rfind("# BLAST processed")
        .expect("outfmt 7's last line");
    let message =
        "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch";
    for threads in [1, 4] {
        let run = search("tests/fasta/x4/chunk_second.faa", "6", threads, &env);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert!(run.stderr.contains(message), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        let commented = search("tests/fasta/x4/chunk_second.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}: {}", commented.stderr);
        assert!(
            commented.stderr.contains(message),
            "n{threads}: {}",
            commented.stderr
        );
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        // In one batch (the default batch size) the refusal comes before the first byte.
        let one_batch = search(
            "tests/fasta/x4/chunk_second.faa",
            "6",
            threads,
            &[("CHUNK_SIZE", "-1"), ("OVERLAP_CHUNK_SIZE", "-200")],
        );
        assert_eq!(one_batch.code, Some(1), "n{threads}: {}", one_batch.stderr);
        assert_eq!(one_batch.stdout, "", "n{threads}");
    }
}
