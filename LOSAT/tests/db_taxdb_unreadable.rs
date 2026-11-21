//! A file of the taxonomy database that cannot be read (line WQA D18b, 2026-11-21): `taxdb.bti`
//! or `taxdb.btd` of mode 000, or a directory in place of either, on `BLASTDB` beside a database
//! whose sequences have taxids, for blastp (a database with parsed Seq-ids, `tm_psi`, and one
//! without, `put`), psiblast, blastx, blastn, tblastn and tblastx by output format, and for
//! blastdbcmd. The expected bytes and exit statuses are NCBI BLAST+ 2.17.0's, run as a user that
//! is not the superuser (tests/seqdb/expected/taxdb_unreadable.tsv, frozen from probe 12 of D18b;
//! NCBI is not run here). NCBI's every load of a sequence with a taxid throws there; the rows
//! that LOSAT writes as NCBI (`same`) are compared byte for byte at 1 and 4 threads, the others
//! (`refused`) are LOSAT's refusal, exit status 1, which says what NCBI does. A superuser reads a
//! file of mode 000, so there only the rows of a directory are compared.
#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The rows: `id program db query words target status stdout stderr losat`. A text field escapes
/// `\\`, `\n` and `\t`; `\e` is empty; `{W}` is the run's directory. `words` are the arguments
/// after `-db DB -query QUERY`, separated by `|`.
const ROWS: &str = include_str!("seqdb/expected/taxdb_unreadable.tsv");

/// The files of the database `db` and of the query `query`, under tests/.
fn files(db: &str, query: &str) -> Vec<String> {
    let (dir, exts): (&str, &[&str]) = match db {
        "tm_psi" => (
            "db_search",
            &[
                "pdb", "phr", "pin", "psq", "ptf", "pto", "pot", "pjs", "pog", "pos",
            ],
        ),
        "put" => (
            "db_search/multidb",
            &["pdb", "phr", "pin", "psq", "ptf", "pto", "pot"],
        ),
        "nut" => (
            "db_search/multidb",
            &["ndb", "nhr", "nin", "nsq", "ntf", "nto", "not"],
        ),
        other => panic!("db {other}"),
    };
    let mut files: Vec<String> = exts.iter().map(|ext| format!("{dir}/{db}.{ext}")).collect();
    if !query.is_empty() {
        files.push(format!("{dir}/{query}"));
    }
    files.push("seqdb/taxdb/good/taxdb.bti".to_owned());
    files.push("seqdb/taxdb/good/taxdb.btd".to_owned());
    files
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

/// A fresh directory with the row's files, the taxonomy file of `target` made unreadable
/// (`bti`, `btd`: mode 000; `btidir`, `btddir`: a directory in its place).
fn prepare(dir: &Path, db: &str, query: &str, target: &str) {
    if dir.exists() {
        let _ = fs::set_permissions(dir, fs::Permissions::from_mode(0o755));
        fs::remove_dir_all(dir).unwrap();
    }
    fs::create_dir_all(dir.join("home")).unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    for name in files(db, query) {
        let to = dir.join(Path::new(&name).file_name().unwrap());
        fs::copy(source.join(&name), &to).unwrap();
        fs::set_permissions(&to, fs::Permissions::from_mode(0o644)).unwrap();
    }
    let file = dir.join(format!("taxdb.{}", &target[..3]));
    match target {
        "bti" | "btd" => fs::set_permissions(&file, fs::Permissions::from_mode(0o000)).unwrap(),
        "btidir" | "btddir" => {
            fs::remove_file(&file).unwrap();
            fs::create_dir(&file).unwrap();
        }
        other => panic!("target {other}"),
    }
}

/// `losat PROGRAM -db DB [-query QUERY] words... [-num_threads N]` in `dir`, with no NCBI, BLAST or
/// LOSAT variable but `BLASTDB=dir`, and `HOME` a directory without files; the output files of
/// outfmt 13 and 14 (`-out o`) are appended to standard output in name order.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`) writes its worker diagnostics there and needs it at more than one thread.
fn run(dir: &Path, row: &[String], threads: Option<&str>) -> (String, String, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || text.starts_with("NCBI")
            || text.starts_with("BLAST")
        {
            command.env_remove(&name);
        }
    }
    let output = command
        .current_dir(dir)
        .env("HOME", dir.join("home"))
        .env("BLASTDB", dir)
        .args([&row[1], "-db", &row[2]])
        .args(
            ["-query", &row[3]]
                .into_iter()
                .filter(|_| !row[3].is_empty()),
        )
        .args(row[4].split('|').filter(|word| !word.is_empty()))
        .args(threads.map(|n| ["-num_threads", n]).into_iter().flatten())
        .output()
        .unwrap();
    let mut stdout = String::from_utf8(output.stdout).unwrap();
    let mut outputs: Vec<PathBuf> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            let name = path.file_name().unwrap().to_string_lossy();
            name == "o" || name.starts_with("o_")
        })
        .collect();
    outputs.sort();
    for path in outputs {
        stdout.push_str(&fs::read_to_string(&path).unwrap());
        fs::remove_file(&path).unwrap();
    }
    (
        stdout,
        String::from_utf8(output.stderr).unwrap(),
        output.status.code().unwrap_or(-1),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:154-199
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-790
// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:873-886
// (cited at `seqdb::taxdb::TaxDb::open`, `common::database::batch_prefetch` and
// `blastdbcmd::retrieve::Source::tax_info`).
#[test]
fn taxonomy_files_that_cannot_be_read_give_ncbis_bytes_or_losats_refusal() {
    let root = std::env::temp_dir().join(format!("losat_taxdb_unreadable_{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let root: PathBuf = fs::canonicalize(&root).unwrap();
    let probe = root.join("probe");
    fs::write(&probe, "").unwrap();
    fs::set_permissions(&probe, fs::Permissions::from_mode(0o000)).unwrap();
    // The superuser reads files of any mode.
    let superuser = fs::File::open(&probe).is_ok();
    let (mut same, mut refused) = (0, 0);
    for line in ROWS.lines().skip(1) {
        let row: Vec<String> = line.split('\t').map(unescape).collect();
        assert_eq!(row.len(), 10, "{line}");
        let (id, target, losat) = (&row[0], &row[5], &row[9]);
        if superuser && matches!(target.as_str(), "bti" | "btd") {
            continue;
        }
        let dir = root.join(id);
        prepare(&dir, &row[2], &row[3], target);
        let place = |text: &str| text.replace("{W}", dir.to_str().unwrap());
        let blastdbcmd = row[1] == "blastdbcmd";
        let threads: &[Option<&str>] = if blastdbcmd {
            &[None]
        } else {
            &[Some("1"), Some("4")]
        };
        for &threads in threads {
            let (stdout, stderr, status) = run(&dir, &row, threads);
            let at = format!("{id} -num_threads {}", threads.unwrap_or("-"));
            match losat.as_str() {
                "same" => {
                    let ncbi_stderr = place(&row[8]);
                    // NCBI's report of a C++ exception (exit status 255) in LOSAT's words
                    // (`losat::wording`, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 12).
                    let expected_stderr = if row[6] == "255" {
                        losat::wording::losat_text(&ncbi_stderr).into_owned()
                    } else {
                        ncbi_stderr
                    };
                    assert_eq!(stderr, expected_stderr, "{at}: stderr");
                    assert_eq!(stdout, place(&row[7]), "{at}: stdout");
                    assert_eq!(status.to_string(), row[6], "{at}: exit status");
                }
                "refused" => {
                    let file = dir.join(format!("taxdb.{}", &target[..3]));
                    let words = format!(
                        "{} cannot be opened: Permission denied, where NCBI BLAST+ 2.17.0 fails to \
                         load each database sequence with a taxid",
                        file.display()
                    );
                    assert!(stderr.contains(&words), "{at}: {stderr}");
                    assert_eq!(status, 1, "{at}: exit status");
                }
                other => panic!("{id}: losat {other}"),
            }
        }
        if matches!(target.as_str(), "bti" | "btd") {
            let file = dir.join(format!("taxdb.{target}"));
            fs::set_permissions(file, fs::Permissions::from_mode(0o644)).unwrap();
        }
        fs::remove_dir_all(&dir).unwrap();
        if losat == "same" {
            same += 1;
        } else {
            refused += 1;
        }
    }
    let _ = fs::set_permissions(&probe, fs::Permissions::from_mode(0o644));
    fs::remove_dir_all(&root).unwrap();
    eprintln!("{same} rows as NCBI, {refused} refused");
    assert!(
        same + refused >= if superuser { 30 } else { 162 },
        "{same} {refused}"
    );
}
