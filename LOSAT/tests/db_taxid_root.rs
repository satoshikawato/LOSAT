//! `-taxids 1` with `-db` (DI-8, DECISIONS 10). Taxid 1 is its own parent in
//! `taxonomy4blast.sqlite3`, so NCBI BLAST+ 2.17.0's expansion of it does not end; LOSAT skips
//! the taxids it has visited and keeps every descendant of the root. NCBI's result for the
//! nearest input that ends is the expansion of the root's children: in
//! `tests/db_search/taxfilter/root/taxonomy4blast.sqlite3` the only child is 131567, and the
//! sequences of 7227 and 8237 (rows left out of that file) are not descendants. The expected
//! bytes are NCBI's `-taxids 131567` of the manifest rows `db.tf.bp.root.6` and
//! `db.tf.tbn.root.6` (`tests/db_search/README.md`, "Taxid filtering"). NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn engine_dir() -> PathBuf {
    fs::canonicalize(Path::new(env!("CARGO_MANIFEST_DIR"))).unwrap()
}

/// One `losat PROGRAM -db DB -query QUERY words... -outfmt 6 -num_threads 1` in LOSAT/, with
/// `BLASTDB=tests/db_search/taxfilter/root`, no registry file (`NCBI_CONFIG_PATH` names that
/// directory, which has none) and no other NCBI, BLAST or LOSAT variable. Returns stdout,
/// stderr and the exit status.
fn search(program: &str, db: &str, query: &str, words: &[&str]) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    let root = "tests/db_search/taxfilter/root";
    let output = command
        .current_dir(engine_dir())
        .env("BLASTDB", root)
        .env("NCBI_CONFIG_PATH", engine_dir().join(root))
        .arg(program)
        .args(["-db", db, "-query", query])
        .args(words)
        .args(["-outfmt", "6", "-num_threads", "1"])
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

fn frozen(row: &str) -> Vec<u8> {
    fs::read(engine_dir().join(format!("tests/fixtures/db_search/{row}.out"))).unwrap()
}

#[test]
fn taxid_1_keeps_what_ncbi_keeps_for_the_roots_children() {
    let t = "tests/db_search/taxfilter";
    for (program, db, row) in [
        ("blastp", format!("{t}/tfp"), "db.tf.bp.root.6"),
        ("tblastn", format!("{t}/tfn"), "db.tf.tbn.root.6"),
    ] {
        let query = format!("{t}/q.faa");
        let (out, err, rc) = search(program, &db, &query, &["-taxids", "1"]);
        assert_eq!(
            (rc, String::from_utf8_lossy(&err).into_owned()),
            (0, String::new()),
            "{program} -taxids 1"
        );
        assert!(out == frozen(row), "{program} -taxids 1 differs from {row}");
        // The left-out taxids' sequences are gone: the whole database writes other bytes.
        let (whole, _, rc) = search(program, &db, &query, &[]);
        assert_eq!(rc, 0, "{program} without a taxid list");
        assert!(whole != out, "{program}: -taxids 1 kept the whole database");
    }
}
