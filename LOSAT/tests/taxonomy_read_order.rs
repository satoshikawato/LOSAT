//! Where the reports that read their subjects one by one write the taxonomy database's texts
//! (line DBID taxrow; review 3 H2 of line DBID pssmread). Where `PRE_FETCH_SEQS_LIMIT` turns a
//! batch's prefetch off, NCBI's formatter loads each subject when it formats it: at a tabular row
//! that reads it, at its description line (or, for one the table does not show, at its
//! alignment's heading) and at its BLAST XML `Hit`. The first load of a blob looks its taxids up,
//! and NCBI's diagnostics go to `cerr`, which flushes `cout` first: the texts land between the
//! report's lines where both streams go to one file (a terminal, `2>&1`). Each case's expectation is
//! NCBI 2.17.0's merged stream of the same run (`tests/fixtures/taxonomy_read_order/<case>.txt`):
//! `ids/txr_p` and `ids/txr_n` (`make_ids.py --taxrow-only`) hold A (taxid 9606, the best hit) and B
//! (taxid 99990, whose lookup in `tests/seqdb/taxdb/offset_last_past_end` writes the offset error),
//! so the text comes between A's and B's lines; in `ids/prd_txself` B (taxid 0, the best hit) comes
//! before A (99990), and the plural name of B's row opens the damaged `tests/seqdb/taxdb/bad_magic`
//! first.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
//! ```c++
//!        		const CBioseq_Handle& subject_bh =
//!                 scope.GetBioseqHandle(align.GetSeq_id(1));
//! ```

use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// The cases: name, taxonomy database (under `tests/seqdb/taxdb`), program, database and query
/// (under `tests/db_search/ids`), options and `-outfmt`.
const CASES: &[(&str, &str, &str, &str, &str, &[&str], &str)] = &[
    (
        "blastp.f7",
        "offset_last_past_end",
        "blastp",
        "txr_p",
        "txr_q.faa",
        &[],
        "7 qseqid sseqid evalue",
    ),
    (
        "blastp.f0",
        "offset_last_past_end",
        "blastp",
        "txr_p",
        "txr_q.faa",
        &[],
        "0",
    ),
    (
        "blastp.f5",
        "offset_last_past_end",
        "blastp",
        "txr_p",
        "txr_q.faa",
        &[],
        "5",
    ),
    (
        "blastp.f0d1",
        "offset_last_past_end",
        "blastp",
        "txr_p",
        "txr_q.faa",
        &["-num_descriptions", "1"],
        "0",
    ),
    (
        "blastp.sci",
        "bad_magic",
        "blastp",
        "prd_txself",
        "prd_q.faa",
        &[],
        "7 qseqid sseqid sscinames",
    ),
    (
        "psiblast.f7",
        "offset_last_past_end",
        "psiblast",
        "txr_p",
        "txr_q.faa",
        &["-num_iterations", "1"],
        "7 qseqid sseqid evalue",
    ),
    (
        "psiblast.f0",
        "offset_last_past_end",
        "psiblast",
        "txr_p",
        "txr_q.faa",
        &["-num_iterations", "1"],
        "0",
    ),
    (
        "blastx.f7",
        "offset_last_past_end",
        "blastx",
        "txr_p",
        "txr_q.fna",
        &[],
        "7 qseqid sseqid evalue",
    ),
    (
        "blastx.f0",
        "offset_last_past_end",
        "blastx",
        "txr_p",
        "txr_q.fna",
        &[],
        "0",
    ),
    (
        "blastx.f5",
        "offset_last_past_end",
        "blastx",
        "txr_p",
        "txr_q.fna",
        &[],
        "5",
    ),
    (
        "tblastn.f7",
        "offset_last_past_end",
        "tblastn",
        "txr_n",
        "txr_q.faa",
        &[],
        "7 qseqid sseqid evalue",
    ),
    (
        "tblastn.f0",
        "offset_last_past_end",
        "tblastn",
        "txr_n",
        "txr_q.faa",
        &[],
        "0",
    ),
    (
        "tblastn.f5",
        "offset_last_past_end",
        "tblastn",
        "txr_n",
        "txr_q.faa",
        &[],
        "5",
    ),
    (
        "blastn.f7",
        "offset_last_past_end",
        "blastn",
        "txr_n",
        "txr_q.fna",
        &["-task", "blastn"],
        "7 qseqid sseqid evalue",
    ),
    (
        "blastn.f0",
        "offset_last_past_end",
        "blastn",
        "txr_n",
        "txr_q.fna",
        &["-task", "blastn"],
        "0",
    ),
    (
        "blastn.f5",
        "offset_last_past_end",
        "blastn",
        "txr_n",
        "txr_q.fna",
        &["-task", "blastn"],
        "5",
    ),
    (
        "tblastx.f7",
        "offset_last_past_end",
        "tblastx",
        "txr_n",
        "txr_q.fna",
        &[],
        "7 qseqid sseqid evalue",
    ),
    (
        "tblastx.f0",
        "offset_last_past_end",
        "tblastx",
        "txr_n",
        "txr_q.fna",
        &[],
        "0",
    ),
    (
        "tblastx.f5",
        "offset_last_past_end",
        "tblastx",
        "txr_n",
        "txr_q.fna",
        &[],
        "5",
    ),
];

/// Runs a case at `threads` threads with standard output and standard error into one file (in a
/// fresh directory of the target's temporary space) and returns that file's bytes.
fn merged(case: &(&str, &str, &str, &str, &str, &[&str], &str), threads: &str) -> Vec<u8> {
    let (name, taxdb, program, db, query, options, outfmt) = *case;
    let dir = std::env::temp_dir().join(format!(
        "losat-taxonomy-read-order-{}-{name}-{threads}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("merged");
    let file = std::fs::File::create(&path).unwrap();
    let mut command = Command::new(BIN);
    command
        .current_dir(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("db_search"),
        )
        .arg(program)
        .args(["-db", &format!("ids/{db}")])
        .args(["-query", &format!("ids/{query}")])
        .args(options)
        .args(["-outfmt", outfmt, "-num_threads", threads])
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
    command
        .env("BLASTDB", format!("../seqdb/taxdb/{taxdb}"))
        .env("PRE_FETCH_SEQS_LIMIT", "0")
        .env("HOME", &dir)
        .env("USERPROFILE", &dir);
    let status = command.status().expect("run LOSAT CLI");
    let bytes = std::fs::read(&path).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(
        status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&bytes)
    );
    bytes
}

#[test]
fn the_texts_of_the_reports_reads_come_where_ncbi_reads_the_subjects() {
    for case in CASES {
        let expected = std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join("taxonomy_read_order")
                .join(format!("{}.txt", case.0)),
        )
        .unwrap();
        for threads in ["1", "4"] {
            let merged = merged(case, threads);
            assert!(
                merged == expected,
                "{} at {threads} threads:\n{}",
                case.0,
                String::from_utf8_lossy(&merged)
            );
        }
    }
}
