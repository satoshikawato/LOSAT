//! `losat makeblastdb` on a full disk: the index file of a volume (and every other output file
//! that passes the limit) is cut short, against the outcomes of makeblastdb 2.17.0 (the official
//! Linux x64 release) under the same limit, frozen in `fixtures/makeblastdb_failure/fsize_rows.txt`.
//!
//! makeblastdb opens the database it built (`CSeqDB`) before it writes the metadata file, and
//! NCBI's reader takes the index file as it finds it: a file that ends early reads as zeros to the
//! end of its last page, so a file of 1 or 2 bytes is "Not a valid version 4 or 5 database.", one
//! of 7 bytes holds the wrong sequence type, one of 8 or 9 bytes is a database with no title, no
//! date and no sequences, and one of no bytes is a fault of NCBI's process (killed by SIGSEGV,
//! which LOSAT refuses: "NCBI crashes on it"). An alias file that is cut short names no database
//! ("No database names were found in alias file [...]", without the file's path).
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1296-1318,
//! c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:134-236 (`CSeqDBIdxFile`).
//!
//! The limit is `prlimit --fsize` with SIGXFSZ ignored (the child of `sh` keeps the ignored signal);
//! the test is skipped where `prlimit` is missing. The cases that cut the date string of the index
//! (the outcome of which depends on the hour of the day) are in the date table of
//! `seqdb::date`'s tests.

#![cfg(target_os = "linux")]

use std::collections::BTreeMap;
use std::process::Command;

const ROWS: &str = include_str!("fixtures/makeblastdb_failure/fsize_rows.txt");

const RICH_P: &str = ">gi|1001|ref|NP_000001.1| alpha\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n>gnl|ti|12345 trace one\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA\n>sp|P55555|CCC_HUMAN Gamma\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARR\n";

/// 130 proteins of 60 residues.
fn many() -> String {
    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut out = String::new();
    for i in 0..130usize {
        out.push_str(&format!(">sp|P{i:05}|S{i} desc {i}\n"));
        out.extend((0..60usize).map(|j| AA[(i * 7 + j * 3 + j * j % 11) % 20] as char));
        out.push('\n');
    }
    out
}

/// A row of the fixture: what NCBI's makeblastdb did.
#[derive(Default)]
struct Row {
    kind: String,
    limit: u64,
    rc: i32,
    err: String,
    out: String,
    files: String,
    pjs: String,
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

fn rows() -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();
    for line in ROWS.lines().filter(|line| !line.starts_with('#')) {
        let (field, value) = line.split_once('\t').unwrap();
        if field == "case" {
            let (kind, limit) = value.split_once('\t').unwrap();
            rows.push(Row {
                kind: kind.to_string(),
                limit: limit.parse().unwrap(),
                ..Row::default()
            });
            continue;
        }
        let row = rows.last_mut().unwrap();
        match field {
            "rc" => row.rc = value.parse().unwrap(),
            "err" => row.err = unescape(value),
            "out" => row.out = unescape(value),
            "files" => row.files = value.to_string(),
            "pjs" => row.pjs = unescape(value),
            other => panic!("unknown field {other}"),
        }
    }
    rows
}

fn has_prlimit() -> bool {
    Command::new("prlimit")
        .arg("--version")
        .output()
        .is_ok_and(|output| output.status.success())
}

/// The time and the duration of the log and the build's date of the metadata masked.
fn masked(text: &str, d: &str) -> String {
    let mut out = text.replace(d, "<D>");
    for (start, end, fill) in [
        ("current time: ", "\n", "<T>"),
        (" sequences in ", " seconds.", "<S>"),
        ("\"last-updated\": \"", "\"", "<ISO>"),
    ] {
        let mut at = 0;
        while let Some(i) = out[at..].find(start) {
            let from = at + i + start.len();
            let Some(j) = out[from..].find(end) else {
                break;
            };
            out.replace_range(from..from + j, fill);
            at = from + fill.len();
        }
    }
    out
}

fn dated(name: &str) -> bool {
    [".pin", ".nin", ".pjs", ".njs", ".pal", ".nal"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

/// The number of files of the output directory, then `name:size` of its index, alias and metadata
/// files (`*` for the size of a file that holds the build's date), as the fixture holds them.
fn listing(dir: &std::path::Path) -> String {
    let mut entries = BTreeMap::new();
    for entry in std::fs::read_dir(dir).unwrap().map_while(Result::ok) {
        let name = entry.file_name().to_string_lossy().into_owned();
        let size = entry.metadata().unwrap().len();
        entries.insert(
            name.clone(),
            if dated(&name) {
                "*".to_string()
            } else {
                size.to_string()
            },
        );
    }
    let chosen: Vec<String> = entries
        .iter()
        .filter(|(name, _)| {
            [".pin", ".nin", ".pal", ".nal", ".pjs", ".njs"]
                .iter()
                .any(|ext| name.ends_with(ext))
        })
        .map(|(n, s)| format!("{n}:{s}"))
        .collect();
    format!("{} {}", entries.len(), chosen.join(" "))
}

#[test]
fn a_full_disk_leaves_the_index_as_ncbi_reads_it() {
    if !has_prlimit() {
        eprintln!("prlimit is not installed: skipped");
        return;
    }
    let mut compared = 0;
    for row in rows() {
        let tmp = std::fs::canonicalize(std::env::temp_dir())
            .unwrap()
            .join(format!(
                "losat-mkdb-failed-index-{}-{}-{}",
                row.kind,
                row.limit,
                std::process::id()
            ));
        let _ = std::fs::remove_dir_all(&tmp);
        let out_dir = tmp.join("o");
        std::fs::create_dir_all(&out_dir).unwrap();
        let (input, args): (String, Vec<&str>) = match row.kind.as_str() {
            "p4" => (
                RICH_P.to_string(),
                vec!["-hash_index", "-blastdb_version", "4"],
            ),
            "m6" => (
                many(),
                vec!["-max_file_sz", "2000B", "-blastdb_version", "4"],
            ),
            other => panic!("unknown kind {other}"),
        };
        std::fs::write(tmp.join("in.fa"), input).unwrap();
        let mut command = Command::new("sh");
        command
            .current_dir(&tmp)
            .env("HOME", &tmp)
            .env_remove("BLASTDB")
            .arg("-c")
            .arg(format!(
                "trap '' XFSZ; exec prlimit --fsize={} -- \"$@\"",
                row.limit
            ))
            .arg("sh")
            .arg(env!("CARGO_BIN_EXE_losat"))
            .args([
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                "prot",
                "-parse_seqids",
            ])
            .args(&args)
            .args(["-out", "o/db"]);
        let output = command.output().expect("run losat makeblastdb");
        let d = tmp.to_str().unwrap();
        let status = output.status.code().unwrap_or(-1);
        let stderr = masked(&String::from_utf8_lossy(&output.stderr), d);
        let stdout = masked(&String::from_utf8_lossy(&output.stdout), d);
        let files = listing(&out_dir);
        let name = format!("{} limit {}", row.kind, row.limit);
        if row.rc == -11 {
            // NCBI's process is killed by SIGSEGV where it reads past its mapping: LOSAT refuses
            // (explicit rejection, exit status 1) and leaves the files NCBI leaves.
            assert_eq!(status, 1, "{name}");
            assert!(
                stderr.starts_with("losat makeblastdb: <D>/o/db.pin: ")
                    && stderr.ends_with("(NCBI crashes on it)\n"),
                "{name}: {stderr}"
            );
        } else {
            assert_eq!(status, row.rc, "{name}");
            let expected_err =
                String::from_utf8(losat::wording::losat_bytes(row.err.as_bytes()).into_owned())
                    .unwrap();
            assert_eq!(stderr, expected_err, "{name}");
        }
        assert_eq!(stdout, row.out, "{name}");
        let pjs: String = ["db.pjs", "db.njs"]
            .iter()
            .filter_map(|n| std::fs::read_to_string(out_dir.join(n)).ok())
            .collect();
        assert_eq!(masked(&pjs, d), row.pjs, "{name}");
        assert_eq!(files, row.files, "{name}");
        let _ = std::fs::remove_dir_all(&tmp);
        compared += 1;
    }
    assert_eq!(compared, 13);
}
