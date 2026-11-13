//! The outfmt 0 of a database without parsed Seq-ids whose sequence has an empty title and taxid
//! 9606 (`tests/db_search/ids/nbare_t`, `pbare_t`; line DBID H3), for the programs that write the
//! head of their report before they read a reported sequence's title (BLASTN, TBLASTX, BLASTX; the
//! `db.id.bare.*` rows of `db_search_manifest/` have NCBI BLAST+ 2.17.0's bytes for all six): NCBI
//! makes the title from the taxonomy database's name of the taxid (`Homo sapiens`, `unnamed protein
//! product [Homo sapiens]`), as LOSAT does since line RG item T5; with a taxonomy database that
//! names no taxid (`no_btd`) the title stays empty; and where the lookup of the taxid reads past
//! the mapped end of `taxdb.btd` (`btd_empty`: NCBI reads unmapped memory and crashes) LOSAT
//! refuses.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// The programs with their database, query, the heading NCBI makes with the good taxonomy
/// database, and the name in LOSAT's refusals.
const PROGRAMS: [(&str, &str, &str, &str, &str); 3] = [
    (
        "blastn",
        "nbare_t",
        "ndup_q.fna",
        ">Homo sapiens\n",
        "BLASTN",
    ),
    (
        "tblastx",
        "nbare_t",
        "ndup_q.fna",
        ">Homo sapiens\n",
        "TBLASTX",
    ),
    (
        "blastx",
        "pbare_t",
        "ndup_q.fna",
        ">unnamed protein product [Homo sapiens]\n",
        "BLASTX",
    ),
];

/// One run of `losat <program> -db <db> -query <query> -outfmt 0` from `LOSAT/` with
/// `BLASTDB=tests/seqdb/taxdb/<taxdb>` and without the other variables that change a search;
/// `home` is an empty directory for HOME (no `~/.ncbirc`).
fn run(program: &str, db: &str, query: &str, taxdb: &str, home: &str) -> Output {
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
fn an_empty_title_is_made_from_the_name_of_its_taxid_in_outfmt_0() {
    let home = std::env::temp_dir().join(format!("losat_db_empty_title_{}", std::process::id()));
    std::fs::create_dir_all(&home).unwrap();
    let home = home.to_str().unwrap().to_owned();
    for (program, db, query, heading, name) in PROGRAMS {
        let named = run(program, db, query, "good", &home);
        let stdout = String::from_utf8_lossy(&named.stdout);
        assert!(
            named.status.success(),
            "{program}: {}",
            String::from_utf8_lossy(&named.stderr)
        );
        assert!(stdout.contains(heading), "{program}: {stdout}");
        let unnamed = run(program, db, query, "no_btd", &home);
        assert!(
            unnamed.status.success(),
            "{program}: {}",
            String::from_utf8_lossy(&unnamed.stderr)
        );
        assert!(
            !String::from_utf8_lossy(&unnamed.stdout).contains("Homo sapiens"),
            "{program}"
        );
        let crashes = run(program, db, query, "btd_empty", &home);
        let stderr = String::from_utf8_lossy(&crashes.stderr);
        assert!(!crashes.status.success(), "{program}: {stderr}");
        assert!(
            stderr.contains("reads unmapped memory")
                && stderr.contains(&format!("not supported by LOSAT's {name}")),
            "{program}: {stderr}"
        );
    }
    std::fs::remove_dir_all(&home).unwrap();
}
