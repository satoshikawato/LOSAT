//! psiblast's search in blocks of query batches (`LOSAT_X_BPQBLOCK`, on by default) of a database
//! with a taxonomy database on `BLASTDB` (line DBID P1). psiblast's batches of one query each
//! prefetch their subjects, which looks their taxids up there (`ListedTitles::prefetch_messages`):
//! the checks before the first block also read every sequence's taxid and ask whether a lookup can
//! be refused (`x_titles_allow_blocks`), so that a refusal in a later block does not come after
//! the reports of the earlier ones; where none can be, the blocks write the messages of the search
//! of all queries. `tests/db_search/tm_psi_u` is a database without parsed Seq-ids (the block
//! search's) of the 8 subjects of `tm_psi`, all with taxid 9606 (`tests/db_search/README.md`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262-264
//! ```c++
//!                    results = psiblast->Run();
//!                    BlastFormatter_PreFetchSequenceData(*results, scope,
//!                            m_CmdLineArgs->GetFormattingArgs()->GetFormattedOutputChoice());
//! ```

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// One `losat psiblast -db tests/db_search/tm_psi_u` from `LOSAT/` with
/// `BLASTDB=tests/seqdb/taxdb/<taxdb>`, `LOSAT_X_BPQBLOCK=<block>` (`1`: a block per query batch,
/// `0`: the search of all queries) and `LOSAT_TIMING` (a `[TIMING] blastp_total` line per search
/// of a query set), without the other variables that change a search; `home` is an empty
/// directory for HOME (no `~/.ncbirc`).
fn psiblast(query: &str, taxdb: &str, block: &str, args: &[&str], home: &str) -> Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("psiblast")
        .args(["-db", "tests/db_search/tm_psi_u"])
        .args(["-query", &format!("tests/db_search/{query}")])
        .args(["-num_threads", "1"])
        .args(args);
    for (name, _) in std::env::vars() {
        // `LOSAT_BLASTX_WORKER_LOG` stays: an `--all-features` binary writes its probe there.
        if name.starts_with("LOSAT_X_")
            || name.starts_with("NCBI")
            || matches!(
                name.as_str(),
                "BATCH_SIZE" | "BL2SEQ_LEGACY" | "CHUNK_SIZE" | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("BLASTDB", format!("tests/seqdb/taxdb/{taxdb}"))
        .env("LOSAT_X_BPQBLOCK", block)
        .env("LOSAT_TIMING", "1")
        .env("HOME", home)
        .env("USERPROFILE", home);
    command.output().expect("run LOSAT CLI")
}

/// The searches of query sets that `output` timed.
fn searches(output: &Output) -> usize {
    String::from_utf8_lossy(&output.stderr)
        .lines()
        .filter(|line| line.starts_with("[TIMING] blastp_total"))
        .count()
}

/// `output`'s standard error without the `[TIMING]` lines.
fn messages(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr)
        .lines()
        .filter(|line| !line.starts_with("[TIMING]"))
        .map(|line| format!("{line}\n"))
        .collect()
}

fn home(tag: &str) -> String {
    let home = std::env::temp_dir().join(format!(
        "losat_psiblast_taxdb_blocks_{tag}_{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&home).unwrap();
    home.to_str().unwrap().to_owned()
}

// `tm_psi_qb.faa` holds q3, without hits at `-evalue 1e-10`, then q1, with hits: with an empty
// `taxdb.btd` the lookup of taxid 9606 in the batch of q1 is refused (NCBI reads unmapped memory
// there). Outfmt 7 writes the lines of q3 without hits; the refusal comes before them in blocks
// as in the search of all queries.
#[test]
fn a_lookup_refused_in_a_later_block_is_refused_before_any_report() {
    let home = home("refused");
    for block in ["1", "0"] {
        let run = psiblast(
            "tm_psi_qb.faa",
            "btd_empty",
            block,
            &["-evalue", "1e-10", "-outfmt", "7"],
            &home,
        );
        let stderr = messages(&run);
        assert!(!run.status.success(), "LOSAT_X_BPQBLOCK={block}: {stderr}");
        assert!(
            run.stdout.is_empty(),
            "LOSAT_X_BPQBLOCK={block}: {}",
            String::from_utf8_lossy(&run.stdout)
        );
        assert!(
            stderr.contains("reaches past the mapped end of the file (0 bytes)")
                && stderr.contains("LOSAT's PSIBLAST"),
            "LOSAT_X_BPQBLOCK={block}: {stderr}"
        );
    }
    std::fs::remove_dir_all(&home).unwrap();
}

// A taxonomy database whose lookups write messages and cannot be refused: its open text (a wrong
// magic number) and the offset error of taxid 99990, which no sequence here has (a lookup of 9606
// writes nothing more). The three queries of `tm_psi_q.faa` make three blocks, which write the
// bytes of the search of all queries.
#[test]
fn blocks_write_the_prefetch_messages_of_the_search_of_all_queries() {
    let home = home("same");
    for taxdb in ["bad_magic", "offset_last_past_end", "good"] {
        for outfmt in ["6", "7", "0"] {
            let args = ["-outfmt", outfmt];
            let reference = psiblast("tm_psi_q.faa", taxdb, "0", &args, &home);
            let blocks = psiblast("tm_psi_q.faa", taxdb, "1", &args, &home);
            let label = format!("{taxdb} outfmt {outfmt}");
            assert!(
                reference.status.success(),
                "{label}: {}",
                messages(&reference)
            );
            assert!(blocks.status.success(), "{label}: {}", messages(&blocks));
            assert_eq!(searches(&reference), 1, "{label}");
            assert_eq!(searches(&blocks), 3, "{label}: the block search ran");
            assert_eq!(blocks.stdout, reference.stdout, "{label}");
            assert_eq!(messages(&blocks), messages(&reference), "{label}");
            if taxdb == "bad_magic" {
                assert_eq!(
                    messages(&blocks)
                        .matches(
                            "Error: [psiblast] Error: Tax database file has wrong magic number."
                        )
                        .count(),
                    1,
                    "{label}"
                );
            }
        }
    }
    std::fs::remove_dir_all(&home).unwrap();
}
