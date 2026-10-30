//! BLASTX `-task blastx-fast` on the call ABRicate makes (bin/abricate:129-141): the query comes on
//! stdin and there is no `-query`. The expected bytes are the frozen NCBI BLAST+ 2.17.0 row
//! `db.bxfast.abr.6` of `db_search_manifest.tsv`, the same search with `-query`.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

#[test]
fn abricate_protein_database_call_with_the_query_on_stdin() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let query = std::fs::read(root.join("tests/db_search/bslice_q.fna")).unwrap();
    let expected =
        std::fs::read(root.join("tests/fixtures/db_search/db.bxfast.abr.6.out")).unwrap();
    assert!(!expected.is_empty(), "the search has hits");
    let mut child = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .args([
            "blastx",
            "-task",
            "blastx-fast",
            "-seg",
            "no",
            "-db",
            "tests/db_search/bslice",
            "-outfmt",
            "6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps \
             gapopen stitle",
            "-num_threads",
            "1",
            "-evalue",
            "1E-20",
            "-culling_limit",
            "1",
            "-max_target_seqs",
            "10000",
        ])
        .env_remove("BLASTDB")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&query).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, expected);
    assert!(output.stderr.is_empty());
}
