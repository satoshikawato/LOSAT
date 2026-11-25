//! Names and words that are not UTF-8 on the command line (line LIB-B): a file name of
//! `-query`, `-subject` or `-out` of BLASTX and RPSTBLASTN, a word of a text seqid list, the names
//! that `blastdbcmd -list` finds, a directory of `NCBI_CONFIG_PATH` and the LMDB file named in a
//! version 5 volume (missing, and listed by `blastdbcmd -metadata`), each with the byte 0xE9
//! (`é` of Latin-1), and a member name of a text `-mask_data` file of makeblastdb that NCBI's
//! exception report quotes (`R4`); and the values of dustmasker, segmasker, windowmasker,
//! makembindex, makeprofiledb, convert2blastmask and blastdbcmd (`K10_*`, line RG item K10a: names
//! of files and databases, `BLASTDB` of the environment and of `.ncbirc`, and words that name a
//! format, an algorithm, an entry or a range; the names of an alias file's `MASKLIST`, which
//! blastdbcmd, the searches, makeblastdb's copy and makembindex read; the path of a profile
//! database that rpsblast and rpstblastn name in their errors); and every value of
//! blastdb_aliastool that names a file, a database or a title (`K10d_*`, line RG item K10d: each
//! mode, existing and missing names). NCBI BLAST+ 2.17.0 keeps such a name's bytes (`std::string`) where it opens, lists and writes it; LOSAT wrote U+FFFD, found
//! nothing or failed before. The expected bytes are NCBI's, frozen by `non_utf8_names/make_fixtures.py`
//! (`non_utf8_names/README.md`) from the cases of `non_utf8_names/cases.tsv`, which this test runs
//! the same way: `{D}` in them is the case's directory, and the time of makeblastdb's line
//! `Building a new DB, current time: ` and the seconds of its `added <n> sequences in` are not
//! compared (as `dropin_corpus.py`); the files and
//! directories the run made, changed or removed are, with the contents of the files, the build's
//! dates in them masked (`expected/<case>.files`, `made`, `undated`). `R9` (a `-matrix` name that NCBI's 1,024-byte buffer cuts inside a character) is a case the
//! survey left open.
//!
//! Unix only, and not macOS: a file name there is the system's bytes, and APFS refuses a name that
//! is not UTF-8 (Windows and WASI keep a name as Unicode text; `LOSAT/src/utils/os_bytes.rs`).
#![cfg(all(unix, not(target_os = "macos")))]

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};
use std::process::Command;

mod temp_support;
use temp_support::{expected_losat_bytes, home_and_worker_log, undated_seqid_list, unique_dir};

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
/// `patch`, `symlink`, `cwd`; the directory where the program runs.
fn set_up(setup: &str, d: &Path) -> PathBuf {
    let losat = Path::new(env!("CARGO_MANIFEST_DIR"));
    let db = d.as_os_str().as_bytes();
    let mut cwd = d.to_path_buf();
    for step in setup.split(" ; ") {
        let parts: Vec<&str> = step.split(' ').collect();
        match parts[0] {
            "copy" => {
                std::fs::copy(losat.join(parts[1]), path_in(d, words(parts[2], db))).unwrap();
            }
            "copydir" => {
                let source = losat.join(parts[1]);
                let target = path_in(d, words(parts[2], db));
                std::fs::create_dir_all(&target).unwrap();
                for entry in std::fs::read_dir(&source).unwrap() {
                    let entry = entry.unwrap();
                    if entry.file_type().unwrap().is_file() {
                        std::fs::copy(entry.path(), target.join(entry.file_name())).unwrap();
                    }
                }
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
            "symlink" => std::os::unix::fs::symlink(
                OsString::from_vec(words(parts[1], db)),
                path_in(d, words(parts[2], db)),
            )
            .unwrap(),
            "cwd" => {
                cwd = path_in(d, words(parts[1], db));
                std::fs::create_dir(&cwd).unwrap();
            }
            "mkidx" => {
                let mut cmd = Command::new(BIN);
                cmd.current_dir(&cwd).arg("makembindex");
                for word in &parts[1..] {
                    cmd.arg(OsString::from_vec(words(word, db)));
                }
                let _files = home_and_worker_log(&mut cmd, "non_utf8_mkidx");
                let out = cmd.output().expect("run LOSAT makembindex");
                assert!(
                    out.status.success(),
                    "mkidx failed: {}",
                    String::from_utf8_lossy(&out.stderr)
                );
            }
            other => panic!("unknown step {other}"),
        }
    }
    cwd
}

/// The frozen bytes of `case` (`expected/<case>.<what>`), `{D}` as `d`, in LOSAT's words.
fn expected(case: &str, what: &str, d: &[u8]) -> Vec<u8> {
    let frozen = std::fs::read(fixtures().join(format!("expected/{case}.{what}"))).unwrap();
    without_time(&expected_losat_bytes(replace(&frozen, b"{D}", d)))
}

/// `output` with the time of makeblastdb's `Building a new DB, current time: ` line as `{TIME}`
/// and the seconds of its `added <n> sequences in <seconds> seconds.` as `{SECONDS}`.
fn without_time(output: &[u8]) -> Vec<u8> {
    let output = masked(
        output,
        b"Building a new DB, current time: ",
        b"\n",
        b"{TIME}",
    );
    masked(&output, b" sequences in ", b" seconds.", b"{SECONDS}")
}

/// `output` with the bytes between the first `before` and the next `after` (or the end) as `with`.
fn masked(output: &[u8], before: &[u8], after: &[u8], with: &[u8]) -> Vec<u8> {
    let Some(at) = output.windows(before.len()).position(|w| w == before) else {
        return output.to_vec();
    };
    let start = at + before.len();
    let end = output[start..]
        .windows(after.len())
        .position(|w| w == after)
        .map_or(output.len(), |n| start + n);
    [&output[..start], with, &output[end..]].concat()
}

/// The paths below `d`, `/`-separated, a directory's (or a link to one, which is not followed)
/// with a final `/` (`make_fixtures.py`'s `paths`).
fn paths(d: &Path) -> BTreeSet<Vec<u8>> {
    let mut found = BTreeSet::new();
    let mut stack = vec![(d.to_path_buf(), Vec::new())];
    while let Some((dir, prefix)) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap() {
            let entry = entry.unwrap();
            let mut name = prefix.clone();
            name.extend_from_slice(entry.file_name().as_bytes());
            if entry.path().is_dir() {
                name.push(b'/');
                if !entry.file_type().unwrap().is_symlink() {
                    stack.push((entry.path(), name.clone()));
                }
            }
            found.insert(name);
        }
    }
    found
}

/// The paths below `d` (`paths`), each with the bytes of its file (`None` for a directory, or a
/// link to one) (`make_fixtures.py`'s `snapshot`).
fn snapshot(d: &Path) -> BTreeMap<Vec<u8>, Option<Vec<u8>>> {
    paths(d)
        .into_iter()
        .map(|path| {
            let full = path_in(d, path.clone());
            let data = (path.last() != Some(&b'/') && full.is_file())
                .then(|| std::fs::read(&full).unwrap());
            (path, data)
        })
        .collect()
}

/// A database's index files, whose header holds the build's date in a field of its own length
/// (`make_fixtures.py`'s `INDEX`).
const INDEX: [&str; 2] = [".nin", ".pin"];

/// The span of the date field of an index file's header (its 4-byte length and its bytes with
/// their padding: after the version, the type, in version 5 the volume, the title and in version
/// 5 the LMDB file's name), or `None` for a file without such a header (`make_fixtures.py`'s
/// `index_date`).
fn index_date(data: &[u8]) -> Option<(usize, usize)> {
    let int4 = |at: usize| {
        data.get(at..at + 4)
            .map(|bytes| u32::from_be_bytes(bytes.try_into().unwrap()) as usize)
    };
    let version = int4(0)?;
    let (mut at, strings) = match version {
        4 => (8, 1),
        5 => (12, 2),
        _ => return None,
    };
    for _ in 0..strings {
        at += 4 + int4(at)?;
    }
    let end = at + 4 + int4(at)?;
    (end <= data.len()).then_some((at, end))
}

/// Whether `data` starts with a text of the build's local time, of a fixed length: NCBI's `CTime`
/// form (`MM/DD/YYYY hh:mm:ss`: an alias file's comment, the index file of a column or of a GI
/// mask) or the metadata file's (`YYYY-MM-DDThh:mm:ss`, `last-updated`) (`make_fixtures.py`'s
/// `DATE_TEXTS`).
fn starts_with_date_text(data: &[u8]) -> bool {
    let fits = |form: &[u8]| {
        data.len() >= form.len()
            && form.iter().zip(data).all(|(&f, &b)| match f {
                b'9' => b.is_ascii_digit(),
                _ => f == b,
            })
    };
    fits(b"99/99/9999 99:99:99") || fits(b"9999-99-99T99:99:99")
}

/// `data`, the bytes of the file `path` below `d`, with the build's dates as `<DATE>`: the date
/// field of an index file (`index_date`) and every date text (`starts_with_date_text`); in a
/// metadata file (`.njs`, `.pjs`) `bytes-total` counts the index files it lists without their
/// date fields, whose length changes with the date; in a binary seqid list (`.bsl`) its date
/// (`temp_support::undated_seqid_list`) (`make_fixtures.py`'s `undated`).
fn undated(d: &Path, path: &[u8], data: &[u8]) -> Vec<u8> {
    let mut data = data.to_vec();
    if path.ends_with(b".bsl") {
        if let Some(undated) = undated_seqid_list(&data) {
            data = undated;
        }
    }
    if INDEX.iter().any(|end| path.ends_with(end.as_bytes())) {
        if let Some((start, end)) = index_date(&data) {
            data = [&data[..start], b"<DATE>", &data[end..]].concat();
        }
    } else if path.ends_with(b".njs") || path.ends_with(b".pjs") {
        let directory = path_in(d, path.to_vec()).parent().unwrap().to_path_buf();
        let key = b"\"files\": [";
        let listed = data
            .windows(key.len())
            .position(|w| w == key)
            .map(|at| &data[at + key.len()..])
            .and_then(|rest| rest.split(|&b| b == b']').next())
            .unwrap_or_default();
        let mut dates = 0;
        for name in listed.split(|&b| b == b'"').skip(1).step_by(2) {
            if INDEX.iter().any(|end| name.ends_with(end.as_bytes())) {
                let candidate = directory.join(OsString::from_vec(name.to_vec()));
                let index = if candidate.exists() {
                    std::fs::read(&candidate).unwrap()
                } else {
                    let base = name.rsplit(|&b| b == b'/').next().unwrap_or(name);
                    std::fs::read(directory.join(OsString::from_vec(base.to_vec()))).unwrap()
                };
                let (start, end) = index_date(&index).unwrap();
                dates += end - start;
            }
        }
        let key = b"\"bytes-total\": ";
        if let Some(at) = data.windows(key.len()).position(|w| w == key) {
            let start = at + key.len();
            let end = start
                + data[start..]
                    .iter()
                    .take_while(|b| b.is_ascii_digit())
                    .count();
            let total: usize = std::str::from_utf8(&data[start..end])
                .unwrap()
                .parse()
                .unwrap();
            let total = (total - dates).to_string();
            data = [&data[..start], total.as_bytes(), &data[end..]].concat();
        }
    }
    if path.ends_with(b".log") {
        data = without_time(&data);
    }
    if path.ends_with(b".nal") || path.ends_with(b".pal") {
        data = masked(&data, b"# Alias file created: ", b"\n", b"<DATE>");
    }
    let mut out = Vec::with_capacity(data.len());
    let mut i = 0;
    while i < data.len() {
        if data[i].is_ascii_digit() && starts_with_date_text(&data[i..]) {
            out.extend_from_slice(b"<DATE>");
            i += "MM/DD/YYYY hh:mm:ss".len();
        } else {
            out.push(data[i]);
            i += 1;
        }
    }
    out
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// Each path below `d` that the run made, changed or removed (`before` and `after` are
/// `snapshot`s), one per line in byte order, a byte that is not printable ASCII or is the
/// backslash as `\xHH`; a file's line then has its size and the FNV-1a 64 of its bytes, the
/// build's dates as `<DATE>` (`undated`) and `d` as `{D}`, after tabs; a removed path's line ends
/// with a tab and `gone` (`make_fixtures.py`'s `made`, `expected/<case>.files`).
fn made(
    d: &Path,
    before: &BTreeMap<Vec<u8>, Option<Vec<u8>>>,
    after: &BTreeMap<Vec<u8>, Option<Vec<u8>>>,
) -> Vec<u8> {
    let mut out = Vec::new();
    let all: BTreeSet<&Vec<u8>> = before.keys().chain(after.keys()).collect();
    for path in all {
        if before.get(path).is_some() && before.get(path) == after.get(path) {
            continue;
        }
        for &b in path {
            if b.is_ascii_graphic() && b != b'\\' {
                out.push(b);
            } else {
                out.extend_from_slice(format!("\\x{b:02x}").as_bytes());
            }
        }
        match after.get(path) {
            None => out.extend_from_slice(b"\tgone"),
            Some(Some(data)) => {
                let data = replace(&undated(d, path, data), d.as_os_str().as_bytes(), b"{D}");
                out.extend_from_slice(
                    format!("\t{}\t{:016x}", data.len(), fnv1a64(&data)).as_bytes(),
                );
            }
            Some(None) => {}
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
        let cwd = set_up(setup, &d);
        let mut words_of_argv = argv.split(' ');
        let program = words_of_argv.next().unwrap();
        let mut command = Command::new(BIN);
        command.current_dir(&cwd).arg(program);
        for word in words_of_argv {
            command.arg(OsString::from_vec(words(word, &db)));
        }
        for name in [
            "BLASTDB",
            "NCBI_CONFIG_PATH",
            "NCBI_DONT_USE_LOCAL_CONFIG",
            "NCBI",
            "WINDOW_MASKER_PATH",
        ] {
            command.env_remove(name);
        }
        if env != "-" {
            let (name, value) = env.split_once('=').unwrap();
            command.env(name, OsString::from_vec(words(value, &db)));
        }
        let _files = home_and_worker_log(&mut command, &format!("non_utf8_{case}"));
        let before = snapshot(&d);
        let output = command.output().expect("run LOSAT");
        let made = made(&d, &before, &snapshot(&d));
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
    assert_eq!(ran, 216, "the cases of cases.tsv");
    assert!(failures.is_empty(), "{}", failures.join("\n\n"));
}
