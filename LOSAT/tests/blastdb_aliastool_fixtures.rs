//! `losat blastdb_aliastool` against NCBI BLAST+ 2.17.0 blastdb_aliastool: each row of
//! `tests/fixtures/blastdb_aliastool/cases.tsv` runs in its own copy of `in/` and gives NCBI's
//! standard output, standard error (in LOSAT's words, `losat::wording`; the run's directory
//! replaced by `<RUN>`) and exit status (`ncbi/<row>.{out,err,rc}`) and, for each file the run made,
//! changed or removed (in `in/` too), its size and FNV-1a digest (`ncbi/<row>.files`; the date of an
//! alias file and of a binary seqid list replaced by `<DATE>`, its size `*`), frozen from the oracle
//! of row RG-D2a (2026-11-24; README.md of the directory). NCBI is not run here.
//!
//! The rows of alias files that name a seqid or taxid list over a version 5 database
//! (`AWAITING_V5_LISTS`) are those that LOSAT's reader refuses until line RG-D (D1c-2) reads such
//! lists: there LOSAT must refuse and leave no alias file.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

mod temp_support;

/// `tests/fixtures/blastdb_aliastool`.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdb_aliastool")
}

/// The rows whose alias file names a seqid or taxid list over a version 5 database, which LOSAT's
/// reader refuses to open until line RG-D reads them (D1c-2, `seqdb::alias::refuse_v5_lists`).
const AWAITING_V5_LISTS: [&str; 10] = [
    "S02", "S03", "S05", "S07", "S09", "S10", "T01", "T02", "T05", "T06",
];

/// The rows of the table: the name and the command-line words.
fn cases() -> Vec<(String, Vec<String>)> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let (name, arguments) = line.split_once('\t').unwrap_or((line, ""));
            (name.to_string(), words(arguments))
        })
        .collect()
}

/// The words of `text` as bash splits them: spaces separate words, single quotes keep their text
/// (the table's only quoting).
fn words(text: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut word: Option<String> = None;
    let mut quoted = false;
    for c in text.chars() {
        match c {
            '\'' => {
                quoted = !quoted;
                word.get_or_insert_with(String::new);
            }
            ' ' if !quoted => {
                if let Some(word) = word.take() {
                    words.push(word);
                }
            }
            _ => word.get_or_insert_with(String::new).push(c),
        }
    }
    words.extend(word);
    words
}

/// FNV-1a 64 of `bytes`.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// The bytes of a file that holds the run's date with the date replaced by `<DATE>`: the
/// `MM/DD/YYYY hh:mm:ss` of an alias file's header, and in a binary seqid list its date
/// (`temp_support::undated_seqid_list`: its size field zeroed only when it holds the file's
/// length); `None` for another file.
fn masked(name: &str, bytes: &[u8]) -> Option<Vec<u8>> {
    if [".pal", ".nal", ".pvl", ".nvl"]
        .iter()
        .any(|extension| name.ends_with(extension))
    {
        let marker = b"# Alias file created ";
        let at = bytes.windows(marker.len()).position(|w| w == marker)? + marker.len();
        let mut out = bytes[..at].to_vec();
        out.extend_from_slice(b"<DATE>");
        out.extend_from_slice(&bytes[at + "MM/DD/YYYY hh:mm:ss".len()..]);
        return Some(out);
    }
    if name.ends_with(".bsl") {
        return temp_support::undated_seqid_list(bytes);
    }
    None
}

/// A command of `binary` with no LOSAT, NCBI or BLAST variable, run in `dir` in UTC.
fn command(binary: &Path, dir: &Path) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .env("HOME", dir)
        .env("TZ", "UTC")
        .current_dir(dir);
    command
}

fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// Every file below `dir` with its bytes, by its `/`-separated path relative to `dir`.
fn listing(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut found = BTreeMap::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else {
                let relative = path.strip_prefix(dir).unwrap();
                let parts: Vec<String> = relative
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                found.insert(parts.join("/"), fs::read(&path).unwrap());
            }
        }
    }
    found
}

/// What a run gave: stdout, stderr, exit status, and the files it made, changed or removed
/// (`path\tsize\tdigest`).
type Run = (Vec<u8>, Vec<u8>, i32, String);

/// `bytes` with the run's directory `dir` replaced by `<RUN>`.
fn with_run(bytes: &[u8], dir: &Path) -> Vec<u8> {
    let run = dir.as_os_str().as_encoded_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i..].starts_with(run) {
            out.extend_from_slice(b"<RUN>");
            i += run.len();
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    out
}

/// The run of `losat blastdb_aliastool words` in a fresh copy of `in/`, in the directory `dir`.
fn run(dir: &Path, words: &[String]) -> Run {
    let _ = fs::remove_dir_all(dir);
    copy_tree(&fixtures().join("in"), &dir.join("in"));
    let dir = dir.canonicalize().unwrap();
    let before = listing(&dir);
    let output = command(Path::new(env!("CARGO_BIN_EXE_losat")), &dir)
        .arg("blastdb_aliastool")
        .args(words)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let after = listing(&dir);
    let mut files = String::new();
    for (path, bytes) in &after {
        if before.get(path) == Some(bytes) {
            continue;
        }
        match masked(path, bytes) {
            Some(masked) => files.push_str(&format!("{path}\t*\t{:016x}\n", fnv1a64(&masked))),
            None => files.push_str(&format!(
                "{path}\t{}\t{:016x}\n",
                bytes.len(),
                fnv1a64(bytes)
            )),
        }
    }
    for path in before.keys().filter(|path| !after.contains_key(*path)) {
        files.push_str(&format!("{path}\tREMOVED\t-\n"));
    }
    (
        with_run(&output.stdout, &dir),
        with_run(&output.stderr, &dir),
        output.status.code().unwrap(),
        files,
    )
}

/// NCBI's run of the row `name`, its standard error in LOSAT's words (`losat::wording`).
fn expected(name: &str) -> Run {
    let file = |extension: &str| fixtures().join("ncbi").join(format!("{name}.{extension}"));
    let err = losat::wording::losat_bytes(&fs::read(file("err")).unwrap()).into_owned();
    let out = losat::wording::losat_output_bytes(&fs::read(file("out")).unwrap()).into_owned();
    (
        out,
        err,
        fs::read_to_string(file("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap(),
        fs::read_to_string(file("files")).unwrap(),
    )
}

fn scratch(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("blastdb_aliastool-{name}-{}", std::process::id()))
}

/// The text of `bytes` for a report, cut to 800 characters.
fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).chars().take(800).collect()
}

#[test]
fn blastdb_aliastool_gives_ncbis_bytes() {
    let mut differ = Vec::new();
    for (name, words) in cases() {
        let dir = scratch(&name);
        let got = run(&dir, &words);
        if AWAITING_V5_LISTS.contains(&name.as_str()) {
            // LOSAT's refusal, before the alias file is kept.
            let refused = got.2 == 1
                && got.1.starts_with(b"Error: ")
                && text(&got.1).contains("not supported by LOSAT")
                && !got.3.contains(".pal\t");
            if !refused {
                eprintln!(
                    "--- {name}: expected LOSAT's refusal, got status {}\nstderr: {:?}\nfiles:\n{}",
                    got.2,
                    text(&got.1),
                    got.3
                );
                differ.push(name);
            } else {
                let _ = fs::remove_dir_all(&dir);
            }
            continue;
        }
        let want = expected(&name);
        if got != want {
            eprintln!(
                "--- {name}: LOSAT status {}, NCBI {}\nLOSAT stdout: {:?}\nNCBI stdout:  {:?}\n\
                 LOSAT stderr: {:?}\nNCBI stderr:  {:?}\nLOSAT files:\n{}NCBI files:\n{}",
                got.2,
                want.2,
                text(&got.0),
                text(&want.0),
                text(&got.1),
                text(&want.1),
                got.3,
                want.3,
            );
            differ.push(name);
        } else {
            let _ = fs::remove_dir_all(&dir);
        }
    }
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}
