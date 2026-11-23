//! `losat blastdbcmd -metadata` where its output cannot be written: the rows of
//! `tests/blastdbcmd/metadata_failures.tsv`, NCBI BLAST+ 2.17.0 blastdbcmd writing the JSON of
//! `db/p4` (with prefixes that make documents of 401 to 21,408 bytes) to standard output or to
//! `-out`, on `/dev/full` or on a file under a size limit (`RLIMIT_FSIZE`, with SIGXFSZ
//! ignored: the write that reaches the limit is short and the next one fails). LOSAT must exit
//! with NCBI's status, post NCBI's lines (in LOSAT's words) and leave NCBI's bytes in the file.
//! The rows need `/dev/full`, `sh` and util-linux's `prlimit`: Linux only.
#![cfg(target_os = "linux")]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::ncbi_windows_paths::canonical;

fn run_dir() -> PathBuf {
    canonical(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastdbcmd"))
}

/// The 64-bit FNV-1a hash of `bytes` (the fixture's `file_fnv1a64`).
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

/// One row: blastdbcmd on `db/p4` with a prefix of `prefix` bytes, its JSON to standard output
/// (`std`) or `-out` in `dir`, which is `/dev/full` (`limit` "full") or a file of at most `limit`
/// bytes; the exit status, the standard error and the file's bytes.
fn run(where_: &str, prefix: usize, limit: &str, dir: &Path) -> (i32, Vec<u8>, Option<Vec<u8>>) {
    let target = if limit == "full" {
        PathBuf::from("/dev/full")
    } else {
        dir.join(format!("{where_}_{prefix}_{limit}.json"))
    };
    let mut words = vec![
        "blastdbcmd".to_string(),
        "-db".to_string(),
        "db/p4".to_string(),
        "-metadata".to_string(),
    ];
    if prefix > 0 {
        words.push("-metadata_output_prefix".to_string());
        words.push("x".repeat(prefix));
    }
    let stdout_target = if where_ == "std" {
        target.clone()
    } else {
        words.push("-out".to_string());
        words.push(target.to_str().unwrap().to_string());
        PathBuf::from("/dev/null")
    };
    let mut command = Command::new("sh");
    command
        .arg("-c")
        .arg(r#"trap '' XFSZ; t=$1; l=$2; shift 2; if [ "$l" = full ]; then exec "$@" > "$t"; else exec prlimit --fsize="$l" -- "$@" > "$t"; fi"#)
        .arg("sh")
        .arg(&stdout_target)
        .arg(limit)
        .arg(env!("CARGO_BIN_EXE_losat"))
        .args(&words);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    let output = command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(run_dir())
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .output()
        .unwrap();
    let file = (limit != "full").then(|| fs::read(&target).unwrap());
    (output.status.code().unwrap(), output.stderr, file)
}

#[test]
fn a_metadata_write_that_fails_is_reported_as_ncbi_reports_it() {
    let dir = std::env::temp_dir().join(format!("losat-metadata-failures-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let rows = fs::read_to_string(run_dir().join("metadata_failures.tsv")).unwrap();
    let mut checked = 0;
    for row in rows.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = row.split('\t').collect();
        let [where_, prefix, limit, status, size, hash, stderr] = fields[..] else {
            panic!("row {row:?}");
        };
        let (code, err, file) = run(where_, prefix.parse().unwrap(), limit, &dir);
        let want_err =
            losat::wording::losat_output_bytes(stderr.replace("\\n", "\n").as_bytes()).into_owned();
        assert_eq!(
            String::from_utf8_lossy(&err),
            String::from_utf8_lossy(&want_err),
            "{row}: stderr"
        );
        assert_eq!(code.to_string(), status, "{row}: exit status");
        if let Some(file) = file {
            assert_eq!(file.len().to_string(), size, "{row}: file size");
            assert_eq!(
                format!("{:016x}", fnv1a64(&file)),
                hash,
                "{row}: file bytes"
            );
        }
        checked += 1;
    }
    assert_eq!(checked, 300);
    let _ = fs::remove_dir_all(&dir);
}
