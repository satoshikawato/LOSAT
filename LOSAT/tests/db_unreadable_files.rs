//! A database file that cannot be read (line WQA D18, 2026-11-21): a file of mode 000 or a
//! directory in place of the file, for the files of version 4 and 5 databases of both molecule
//! types, an alias file and a member of an alias, and the lookup files that `-taxids`,
//! `-negative_taxids`, `-seqidlist` and `-negative_seqidlist` open, in the search programs and in
//! blastdbcmd; and an empty lookup or LMDB file (`.pos`, `.pot`, `.ptf`, `.pdb`: review M4 of D18,
//! probe 10 of its round 2; NCBI's texts are cited at `seqdb::NCBI_FORMAT_MESSAGES`). The expected bytes and exit
//! statuses are NCBI BLAST+ 2.17.0's, run as a user that is not the superuser
//! (tests/seqdb/expected/unreadable_files.tsv, frozen from the probes of D18); NCBI is not run here.
//! A superuser reads a file of mode 000, so there only the rows of a directory are compared.
#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The rows: `id set program db query words file how fmt status stdout stderr`. A text field
/// escapes `\\`, `\n` and `\t`; `\e` is empty; `{W}` is the run's directory.
const ROWS: &str = include_str!("seqdb/expected/unreadable_files.tsv");

/// The files of each set, under tests/db_search/, copied into the run's directory.
fn set_files(set: &str) -> Vec<&'static str> {
    let nucl_v4 = [
        "ll2/nucl_v4.nhr",
        "ll2/nucl_v4.nin",
        "ll2/nucl_v4.njs",
        "ll2/nucl_v4.nnd",
        "ll2/nucl_v4.nni",
        "ll2/nucl_v4.nog",
        "ll2/nucl_v4.nsd",
        "ll2/nucl_v4.nsi",
        "ll2/nucl_v4.nsq",
        "ll2/q.fna",
    ];
    let tm_psi = [
        "tm_psi.pdb",
        "tm_psi.phr",
        "tm_psi.pin",
        "tm_psi.pjs",
        "tm_psi.pog",
        "tm_psi.pos",
        "tm_psi.pot",
        "tm_psi.psq",
        "tm_psi.ptf",
        "tm_psi.pto",
        "x16_taxid_qwarn.faa",
    ];
    match set {
        "P4" => vec!["paja.pin", "paja.phr", "paja.psq", "d1_q.faa"],
        "P5" => tm_psi.to_vec(),
        "N4" | "AN" => nucl_v4.to_vec(),
        "N5" => vec![
            "ll2/nucl_v5.ndb",
            "ll2/nucl_v5.nhr",
            "ll2/nucl_v5.nin",
            "ll2/nucl_v5.njs",
            "ll2/nucl_v5.nog",
            "ll2/nucl_v5.nos",
            "ll2/nucl_v5.not",
            "ll2/nucl_v5.nsq",
            "ll2/nucl_v5.ntf",
            "ll2/q.fna",
        ],
        "AP" => vec![
            "tm_2vol.pal",
            "tm_2vol_a.phr",
            "tm_2vol_a.pin",
            "tm_2vol_a.psq",
            "tm_2vol_b.phr",
            "tm_2vol_b.pin",
            "tm_2vol_b.psq",
            "x16_taxid_qwarn.faa",
        ],
        "TF" => vec![
            "taxfilter/tfe_p.pdb",
            "taxfilter/tfe_p.phr",
            "taxfilter/tfe_p.pin",
            "taxfilter/tfe_p.pog",
            "taxfilter/tfe_p.pos",
            "taxfilter/tfe_p.pot",
            "taxfilter/tfe_p.psq",
            "taxfilter/tfe_p.ptf",
            "taxfilter/tfe_p.pto",
            "taxfilter/q.faa",
            "taxfilter/taxonomy4blast.sqlite3",
        ],
        _ => panic!("set {set}"),
    }
}

fn unescape(field: &str) -> String {
    if field == "\\e" {
        return String::new();
    }
    let mut out = String::with_capacity(field.len());
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            other => panic!("escape {other:?} in {field}"),
        }
    }
    out
}

/// A fresh directory with the set's files (and the alias and list files of the probes), the
/// row's file made unreadable (`how`: mode 000, a directory, or empty).
fn prepare(dir: &Path, set: &str, file: &str, how: &str) {
    if dir.exists() {
        let _ = fs::set_permissions(dir, fs::Permissions::from_mode(0o755));
        fs::remove_dir_all(dir).unwrap();
    }
    fs::create_dir_all(dir.join("home")).unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    for name in set_files(set) {
        let to = dir.join(Path::new(name).file_name().unwrap());
        fs::copy(source.join(name), &to).unwrap();
        fs::set_permissions(&to, fs::Permissions::from_mode(0o644)).unwrap();
    }
    if set == "AN" {
        fs::write(dir.join("aln.nal"), "TITLE aln\nDBLIST nucl_v4\n").unwrap();
    }
    if set == "P5" {
        fs::write(dir.join("list.txt"), "PTX00001.1\nPTX00004.1\n").unwrap();
    }
    match how {
        "000" => fs::set_permissions(dir.join(file), fs::Permissions::from_mode(0o000)).unwrap(),
        "dir" => {
            fs::remove_file(dir.join(file)).unwrap();
            fs::create_dir(dir.join(file)).unwrap();
        }
        "empty" => fs::write(dir.join(file), b"").unwrap(),
        _ => assert_eq!(how, "none"),
    }
}

/// `losat PROGRAM -db DB -query QUERY words... -outfmt FMT -num_threads 1` in `dir` (blastdbcmd:
/// `losat blastdbcmd -db DB words...`), with no NCBI, BLAST or LOSAT variable but `BLASTDB=dir`
/// for the set whose taxonomy file NCBI found there, and `HOME` a directory without files.
fn run(dir: &Path, set: &str, row: &[String]) -> (String, String, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    if set == "TF" {
        command.env("BLASTDB", dir);
    }
    let output = command
        .current_dir(dir)
        .env("HOME", dir.join("home"))
        .args([&row[2], "-db", &row[3]])
        .args(
            ["-query", &row[4]]
                .into_iter()
                .filter(|_| !row[4].is_empty()),
        )
        .args(row[5].split_whitespace())
        .args(
            ["-outfmt", &row[8], "-num_threads", "1"]
                .into_iter()
                .filter(|_| !row[8].is_empty()),
        )
        .output()
        .unwrap();
    (
        String::from_utf8(output.stdout).unwrap(),
        String::from_utf8(output.stderr).unwrap(),
        output.status.code().unwrap_or(-1),
    )
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5893-5947
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:888-895
// (cited at `seqdb::SeqDbError::memory_file_exception`, `SeqDbError::ncbi_request_message`
// and `seqdb::filter::lookup_file_error`).
#[test]
fn database_files_that_cannot_be_read_give_ncbis_bytes() {
    let root =
        std::env::temp_dir().join(format!("losat_db_unreadable_files_{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let root: PathBuf = fs::canonicalize(&root).unwrap();
    let probe = root.join("probe");
    fs::write(&probe, "").unwrap();
    fs::set_permissions(&probe, fs::Permissions::from_mode(0o000)).unwrap();
    // The superuser reads files of any mode.
    let superuser = fs::File::open(&probe).is_ok();
    let mut compared = 0;
    for line in ROWS.lines().skip(1) {
        let row: Vec<String> = line.split('\t').map(unescape).collect();
        assert_eq!(row.len(), 12, "{line}");
        let (id, set, file, how) = (&row[0], &row[1], &row[6], &row[7]);
        if superuser && how == "000" {
            continue;
        }
        let dir = root.join(id);
        prepare(&dir, set, file, how);
        let place = |text: &str| text.replace("{W}", dir.to_str().unwrap());
        let (stdout, stderr, status) = run(&dir, set, &row);
        let ncbi_stderr = place(&row[11]);
        // NCBI's report of a C++ exception (exit status 255) in LOSAT's words
        // (`losat::wording`, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 12).
        let expected_stderr = if row[9] == "255" {
            losat::wording::losat_text(&ncbi_stderr).into_owned()
        } else {
            ncbi_stderr
        };
        assert_eq!(stderr, expected_stderr, "{id}: stderr");
        assert_eq!(stdout, place(&row[10]), "{id}: stdout");
        assert_eq!(status.to_string(), row[9], "{id}: exit status");
        if how == "000" {
            fs::set_permissions(dir.join(file), fs::Permissions::from_mode(0o644)).unwrap();
        }
        fs::remove_dir_all(&dir).unwrap();
        compared += 1;
    }
    assert!(
        compared >= if superuser { 29 } else { 119 },
        "{compared} rows"
    );
    let _ = fs::set_permissions(&probe, fs::Permissions::from_mode(0o644));
    fs::remove_dir_all(&root).unwrap();
}
