//! Where psiblast writes the warnings of its PSSM's reads (line DBID pssmread, review 2 H4): NCBI
//! makes each query's PSSM after the query's report, inside the query loop, and writes the epilog
//! after the loop, so the last query's `Failed to retrieve sequence` warnings come before outfmt 7's
//! `# BLAST processed N queries` and before outfmt 0's database statistics. NCBI's diagnostics go
//! to `cerr`, which flushes `cout` first: the order shows where both streams go to one file (a
//! terminal, `2>&1`). `tests/db_search/ids/prd_mb.pal` holds Y (`gb|BBB00001.1|`, 100 letters) and
//! A (`gb|AAA00001.1|`, the best hit), B (`gb|BBB00001|`, whose Seq-align id the lookup gives back as
//! Y, too short for B's HSP) and C: with `-max_target_seqs 1` no report reads B, and every PSSM of
//! `prd_q3.faa`'s three queries warns once for it (`tests/db_search/README.md`, review 2 of
//! pssmread; NCBI 2.17.0's merged stream of the same runs is byte-identical to LOSAT's).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-308 and 562
//! ```c++
//!                         pssm =
//!                          ComputePssmForNextIteration(*seq, aln, psi_opts, scope,
//!                                      results_1st_query.GetAncillaryData());
//!     ...
//!         formatter.PrintEpilog(opt);
//! ```

use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// `losat psiblast -db tests/db_search/ids/prd_mb -query tests/db_search/ids/prd_q3.faa` with one
/// round, its ASCII PSSM, `-max_target_seqs 1` and `-outfmt outfmt`, standard output and standard
/// error into one file
/// (`name` in a fresh directory of the target's temporary space); returns that file's bytes.
fn merged(name: &str, outfmt: &str) -> String {
    let dir = std::env::temp_dir().join(format!(
        "losat-pssm-read-order-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("merged");
    let file = std::fs::File::create(&path).unwrap();
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("psiblast")
        .args(["-db", "tests/db_search/ids/prd_mb"])
        .args(["-query", "tests/db_search/ids/prd_q3.faa"])
        .args(["-num_iterations", "1", "-save_pssm_after_last_round"])
        .arg("-out_ascii_pssm")
        .arg(dir.join("pssm"))
        .args([
            "-max_target_seqs",
            "1",
            "-outfmt",
            outfmt,
            "-num_threads",
            "1",
        ])
        .stdout(Stdio::from(file.try_clone().unwrap()))
        .stderr(Stdio::from(file));
    for (variable, _) in std::env::vars() {
        // `LOSAT_BLASTX_WORKER_LOG` stays: an `--all-features` binary writes its probe there.
        if variable.starts_with("LOSAT_X_")
            || variable.starts_with("NCBI")
            || matches!(
                variable.as_str(),
                "BATCH_SIZE" | "BLASTDB" | "CHUNK_SIZE" | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(variable);
        }
    }
    command.env("HOME", &dir).env("USERPROFILE", &dir);
    let status = command.status().expect("run LOSAT CLI");
    let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(status.success(), "{text}");
    text
}

/// The byte offsets of the lines of `text` that start with `prefix`.
fn lines_at(text: &str, prefix: &str) -> Vec<usize> {
    let mut at = Vec::new();
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        if line.starts_with(prefix) {
            at.push(offset);
        }
        offset += line.len();
    }
    at
}

const WARNING: &str = "Warning: [psiblast] Failed to retrieve sequence gb|BBB00001|";

#[test]
fn the_last_query_s_pssm_warning_comes_before_outfmt_7_s_count() {
    let text = merged("f7", "7 qseqid evalue");
    let warnings = lines_at(&text, WARNING);
    let count = lines_at(&text, "# BLAST processed 3 queries");
    let headers = lines_at(&text, "# PSIBLAST");
    assert_eq!(
        (warnings.len(), count.len(), headers.len()),
        (3, 1, 3),
        "{text}"
    );
    // Each query's warning follows its report and precedes the next query's header; the last
    // one precedes the epilog.
    assert!(
        headers[0] < warnings[0] && warnings[0] < headers[1],
        "{text}"
    );
    assert!(
        headers[1] < warnings[1] && warnings[1] < headers[2],
        "{text}"
    );
    assert!(headers[2] < warnings[2] && warnings[2] < count[0], "{text}");
}

#[test]
fn the_last_query_s_pssm_warning_comes_before_outfmt_0_s_database_statistics() {
    let text = merged("f0", "0");
    let warnings = lines_at(&text, WARNING);
    let footers = lines_at(&text, "Effective search space used:");
    let database = lines_at(&text, "  Database: ");
    assert_eq!((warnings.len(), footers.len()), (3, 3), "{text}");
    // The epilog's database block is the last `  Database: ` line.
    let epilog = *database.last().expect("the epilog's database block");
    assert!(footers[2] < warnings[2] && warnings[2] < epilog, "{text}");
    // Its two blank lines come after the warning (blast_format.cpp:2249).
    assert!(
        text[warnings[2]..epilog].starts_with(&format!("{WARNING}\n\n\n")),
        "{text}"
    );
}
