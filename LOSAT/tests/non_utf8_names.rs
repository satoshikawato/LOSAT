//! Names and words that are not UTF-8 on the command line (line LIB-B): a file name of
//! `-query`, `-subject` or `-out` of BLASTX and RPSTBLASTN, a word of a text seqid list, the names
//! that `blastdbcmd -list` finds, a directory of `NCBI_CONFIG_PATH` and the LMDB file named in a
//! version 5 volume (missing, and listed by `blastdbcmd -metadata`), each with the byte 0xE9
//! (`é` of Latin-1), and a member name of a text `-mask_data` file of makeblastdb that NCBI's
//! exception report quotes (`R4`); and the values of dustmasker, segmasker, windowmasker,
//! makembindex, makeprofiledb, convert2blastmask and blastdbcmd (`K10_*`, line RG item K10a: names
//! of files and databases, `BLASTDB` of the environment and of `.ncbirc`, and words that name a
//! format, an algorithm, an entry or a range). NCBI BLAST+ 2.17.0 keeps such a name's bytes (`std::string`) where it opens, lists and writes it; LOSAT wrote U+FFFD, found
//! nothing or failed before. The expected bytes are NCBI's, frozen by `non_utf8_names/make_fixtures.py`
//! (`non_utf8_names/README.md`) from the cases of `non_utf8_names/cases.tsv`, which this test runs
//! the same way: `{D}` in them is the case's directory, and the time of makeblastdb's line
//! `Building a new DB, current time: ` is not compared (as `dropin_corpus.py`); the names of the
//! files and directories the run made are (`expected/<case>.files`, `made`). `R9` (a `-matrix` name that NCBI's 1,024-byte buffer cuts inside a character) is a case the
//! survey left open.
//!
//! Unix only, and not macOS: a file name there is the system's bytes, and APFS refuses a name that
//! is not UTF-8 (Windows and WASI keep a name as Unicode text; `LOSAT/src/utils/os_bytes.rs`).
#![cfg(all(unix, not(target_os = "macos")))]

use std::collections::BTreeSet;
use std::ffi::OsString;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};
use std::process::Command;

mod temp_support;
use temp_support::{expected_losat_bytes, home_and_worker_log, unique_dir};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/non_utf8_names")
}

/// The bytes of a word of the manifest: `\xHH`, `\n`, `\s` (a space), `\\`; `{D}` is `d`.
fn words(text: &str, d: &[u8]) -> Vec<u8> {
    let raw = text.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        if raw[i] == b'\\' {
            match raw[i + 1] {
                b'x' => {
                    let hex = std::str::from_utf8(&raw[i + 2..i + 4]).unwrap();
                    out.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 4;
                    continue;
                }
                b'n' => out.push(b'\n'),
                b's' => out.push(b' '),
                b'\\' => out.push(b'\\'),
                other => panic!("escape \\{}", other as char),
            }
            i += 2;
            continue;
        }
        out.push(raw[i]);
        i += 1;
    }
    replace(&out, b"{D}", d)
}

fn replace(data: &[u8], from: &[u8], to: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len());
    let mut i = 0;
    while i < data.len() {
        if data[i..].starts_with(from) {
            out.extend_from_slice(to);
            i += from.len();
        } else {
            out.push(data[i]);
            i += 1;
        }
    }
    out
}

fn path_in(d: &Path, name: Vec<u8>) -> PathBuf {
    d.join(OsString::from_vec(name))
}

/// The setup steps of a case (`make_fixtures.py`'s text): `copy`, `copydb`, `mkdir`, `write`,
/// `patch`.
fn set_up(setup: &str, d: &Path) {
    let losat = Path::new(env!("CARGO_MANIFEST_DIR"));
    let db = d.as_os_str().as_bytes();
    for step in setup.split(" ; ") {
        let parts: Vec<&str> = step.split(' ').collect();
        match parts[0] {
            "copy" => {
                std::fs::copy(losat.join(parts[1]), path_in(d, words(parts[2], db))).unwrap();
            }
            "copydb" => {
                let source = losat.join(parts[1]);
                let prefix = format!("{}.", source.file_name().unwrap().to_str().unwrap());
                let mut names: Vec<_> = std::fs::read_dir(source.parent().unwrap())
                    .unwrap()
                    .map(|entry| entry.unwrap().file_name())
                    .filter(|name| name.to_str().is_some_and(|n| n.starts_with(&prefix)))
                    .collect();
                names.sort();
                let target = path_in(d, words(parts[2], db));
                for name in names {
                    std::fs::copy(source.parent().unwrap().join(&name), target.join(&name))
                        .unwrap();
                }
            }
            "mkdir" => std::fs::create_dir(path_in(d, words(parts[1], db))).unwrap(),
            "write" => {
                std::fs::write(path_in(d, words(parts[1], db)), words(parts[2], db)).unwrap()
            }
            "patch" => {
                let path = path_in(d, words(parts[1], db));
                let (old, new) = (words(parts[2], db), words(parts[3], db));
                assert_eq!(old.len(), new.len());
                let data = std::fs::read(&path).unwrap();
                let at = data
                    .windows(old.len())
                    .position(|w| w == old.as_slice())
                    .expect("the bytes to patch");
                let mut patched = data;
                patched[at..at + old.len()].copy_from_slice(&new);
                std::fs::write(&path, patched).unwrap();
            }
            other => panic!("unknown step {other}"),
        }
    }
}

/// The frozen bytes of `case` (`expected/<case>.<what>`), `{D}` as `d`, in LOSAT's words.
fn expected(case: &str, what: &str, d: &[u8]) -> Vec<u8> {
    let frozen = std::fs::read(fixtures().join(format!("expected/{case}.{what}"))).unwrap();
    without_time(&expected_losat_bytes(replace(&frozen, b"{D}", d)))
}

/// `output` with the time of makeblastdb's `Building a new DB, current time: ` line as `{TIME}`.
fn without_time(output: &[u8]) -> Vec<u8> {
    const LINE: &[u8] = b"Building a new DB, current time: ";
    let Some(at) = output.windows(LINE.len()).position(|w| w == LINE) else {
        return output.to_vec();
    };
    let start = at + LINE.len();
    let end = output[start..]
        .iter()
        .position(|&b| b == b'\n')
        .map_or(output.len(), |n| start + n);
    [&output[..start], b"{TIME}".as_slice(), &output[end..]].concat()
}

/// The paths below `d`, `/`-separated, a directory's with a final `/` (`make_fixtures.py`'s
/// `paths`).
fn paths(d: &Path) -> BTreeSet<Vec<u8>> {
    let mut found = BTreeSet::new();
    let mut stack = vec![(d.to_path_buf(), Vec::new())];
    while let Some((dir, prefix)) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap() {
            let entry = entry.unwrap();
            let mut name = prefix.clone();
            name.extend_from_slice(entry.file_name().as_bytes());
            if entry.file_type().unwrap().is_dir() {
                name.push(b'/');
                stack.push((entry.path(), name.clone()));
            }
            found.insert(name);
        }
    }
    found
}

/// Each path of `after` not in `before`, one per line in byte order, a byte that is not printable
/// ASCII or is the backslash as `\xHH` (`make_fixtures.py`'s `made`, `expected/<case>.files`).
fn made(before: &BTreeSet<Vec<u8>>, after: &BTreeSet<Vec<u8>>) -> Vec<u8> {
    let mut out = Vec::new();
    for path in after.difference(before) {
        for &b in path {
            if b.is_ascii_graphic() && b != b'\\' {
                out.push(b);
            } else {
                out.extend_from_slice(format!("\\x{b:02x}").as_bytes());
            }
        }
        out.push(b'\n');
    }
    out
}

#[test]
fn names_that_are_not_utf8_keep_their_bytes_as_in_ncbi() {
    let manifest = std::fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut failures = Vec::new();
    let mut ran = 0;
    for line in manifest.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let [case, kind, env, setup, argv] = fields[..] else {
            panic!("manifest line {line:?}");
        };
        let dir = unique_dir(&format!("non_utf8_{case}"));
        let d = std::fs::canonicalize(dir.path()).unwrap();
        let db = d.as_os_str().as_bytes().to_vec();
        set_up(setup, &d);
        let mut words_of_argv = argv.split(' ');
        let program = words_of_argv.next().unwrap();
        let mut command = Command::new(BIN);
        command.current_dir(&d).arg(program);
        for word in words_of_argv {
            command.arg(OsString::from_vec(words(word, &db)));
        }
        for name in [
            "BLASTDB",
            "NCBI_CONFIG_PATH",
            "NCBI_DONT_USE_LOCAL_CONFIG",
            "NCBI",
        ] {
            command.env_remove(name);
        }
        if env != "-" {
            let (name, value) = env.split_once('=').unwrap();
            command.env(name, OsString::from_vec(words(value, &db)));
        }
        let _files = home_and_worker_log(&mut command, &format!("non_utf8_{case}"));
        let before = paths(&d);
        let output = command.output().expect("run LOSAT");
        let made = made(&before, &paths(&d));
        ran += 1;
        let rc: i32 = std::fs::read_to_string(fixtures().join(format!("expected/{case}.rc")))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        let ok = match kind {
            "ncbi" => {
                output.status.code() == Some(rc)
                    && without_time(&output.stdout) == expected(case, "stdout", &db)
                    && output.stderr == expected(case, "stderr", &db)
                    && made
                        == std::fs::read(fixtures().join(format!("expected/{case}.files"))).unwrap()
            }
            other => panic!("kind {other}"),
        };
        if !ok {
            failures.push(format!(
                "{case}: rc {:?} (NCBI {rc})\nstdout {:?}\nstderr {:?}\nmade {:?}",
                output.status.code(),
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&made)
            ));
        }
    }
    assert_eq!(ran, 82, "the cases of cases.tsv");
    assert!(failures.is_empty(), "{}", failures.join("\n\n"));
}
