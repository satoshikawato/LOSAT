//! BLASTP's refusals of a step where NCBI's search of a later query batch stops (line RG item X4,
//! fix round 2; PD-LOSAT-NCBI-DEFECTS 1.51, decision D13): NCBI prints each batch's reports before
//! it searches the next, so where its search of batch 2 stops with an exception, batch 1's reports
//! are printed (without the epilog). LOSAT refuses after them. `x4/chunk_second.faa`: a protein of
//! PajaWSV (92 letters, `x4/chunk_first.faa`) and the TK-TMK protein of psiblast's rows (398
//! letters); with `BATCH_SIZE=50` each is a batch, and with `CHUNK_SIZE=150` (the overlap 100) the
//! first is not split (92 / 50 makes one chunk) and the second is split into chunks of 143 letters
//! that NCBI would split again, which NCBI stops on with a null reference (`db.ps.it.rej.split`'s
//! case).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             ...
//!                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//!                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//!                 results = lcl_blast.Run();
//! ```

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
