//! The outfmt 0 of a database without parsed Seq-ids whose sequence has an empty title and taxid
//! 9606 (`tests/db_search/ids/nbare_t`, `pbare_t`; line DBID H3), for the programs that write the
//! head of their report before they read a reported sequence's title, so that a refusal leaves
//! output (BLASTN, TBLASTX, BLASTX; the `db.id.bare.*` rows of `db_search_manifest/` have the
//! others): NCBI BLAST+ 2.17.0 makes the title from the taxid where the taxonomy database names it
//! (`Homo sapiens`, `unnamed protein product [Homo sapiens]`), which LOSAT refuses; with a
//! taxonomy database that names no taxid (`no_btd`) it writes the stored title, as LOSAT does.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// The programs with their database and query.
const PROGRAMS: [(&str, &str, &str, &str); 3] = [
    ("blastn", "nbare_t", "ndup_q.fna", "BLASTN"),
    ("tblastx", "nbare_t", "ndup_q.fna", "TBLASTX"),
    ("blastx", "pbare_t", "ndup_q.fna", "BLASTX"),
];

/// One run of `losat <program> -db <db> -query <query> -outfmt 0` from `LOSAT/` with
/// `BLASTDB=tests/seqdb/taxdb/<taxdb>` and without the other variables that change a search;
/// `home` is an empty directory for HOME (no `~/.ncbirc`).
fn run((program, db, query, _): (&str, &str, &str, &str), taxdb: &str, home: &str) -> Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg(program)
        .args(["-db", &format!("tests/db_search/ids/{db}")])
        .args(["-query", &format!("tests/db_search/ids/{query}")])
        .args(["-outfmt", "0", "-num_threads", "1"]);
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
        .env("HOME", home)
        .env("USERPROFILE", home);
    command.output().expect("run LOSAT CLI")
}

#[test]
fn an_empty_title_with_a_named_taxid_is_refused_in_outfmt_0() {
    let home = std::env::temp_dir().join(format!("losat_db_empty_title_{}", std::process::id()));
    std::fs::create_dir_all(&home).unwrap();
    let home = home.to_str().unwrap().to_owned();
    for program in PROGRAMS {
        let refused = run(program, "good", &home);
        let stderr = String::from_utf8_lossy(&refused.stderr);
        assert!(!refused.status.success(), "{}: {stderr}", program.0);
        assert!(
            stderr.contains(
                "database sequence 0 has an empty title that NCBI BLAST+'s outfmt 0 makes from its taxid and the taxonomy database"
            ) && stderr.contains(&format!("not supported by LOSAT's {}", program.3)),
            "{}: {stderr}",
            program.0
        );
        let written = run(program, "no_btd", &home);
        assert!(
            written.status.success(),
            "{}: {}",
            program.0,
            String::from_utf8_lossy(&written.stderr)
        );
    }
    std::fs::remove_dir_all(&home).unwrap();
}
