//! `losat blastdbcheck` against NCBI BLAST+ 2.17.0 blastdbcheck: each row of
//! `tests/fixtures/blastdbcheck/cases.tsv` runs in its own fresh run directory (`in/` copied from
//! `tests/fixtures/blastdb_aliastool/in`, an empty `home/`, and the row's case directory built by
//! the operations of `setup.tsv`) and gives NCBI's standard output, standard error (in LOSAT's
//! words, `losat::wording`; the run's directory replaced by `<RUN>`), exit status
//! (`ncbi/<row>.{out,err,rc}`) and, for each file the run made or changed, its size and FNV-1a
//! digest (`ncbi/<row>.files`), frozen from the oracle of row RG-D2b (2026-11-25; README.md of the
//! directory). NCBI is not run here.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1-1560
//! The rows whose output depends on the time-seeded random sample (`seed` not `-`) are skipped
//! here; a crate-internal test covers them.

#![cfg(unix)]

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// Rows where LOSAT must refuse (exit status 1, standard error starting `Error: ` and containing
/// `not supported by LOSAT`) because NCBI's result is a crash or an endless loop that LOSAT does
/// not reproduce: a sequence file of no bytes, a cut string ISAM, a garbage index found by `-dir`
/// (status 139), `-stride` of 0 or less (killed after 5 s); NCBI's frozen result stays beside it.
const LOSAT_REFUSES: [&str; 7] = [
    "d_p4_psq_zero",
    "d_p5_psq_zero",
    "d_p4_psd_cut",
    "r_garbage",
    "r_garbage2",
    "stride0",
    "strideneg",
];

/// Rows where NCBI reads zeros past the end of a short file inside its last mapped page (a
/// truncated index, sequence or header file) and succeeds; LOSAT must refuse them (as above)
/// until its volume readers model that page (row RG-D2b's next PR).
const AWAITING_ZERO_PAGE: [&str; 6] = [
    "d_p4_pin_10",
    "d_p4_pin_cut",
    "d_p4_psq_half",
    "d_p4_phr_half",
    "d_p5_pin_cut",
    "d_n4_nsq_half",
];

/// Rows of an alias file whose count checks overflow a signed integer in NCBI's binary, where
/// the overflow's result (undefined in C) decides the check; LOSAT must refuse them (as above)
/// until the Owner decides the exception (Q-D2B-OV). An overflow that decides nothing gives NCBI's
/// result (a_ov3). `-dir .` of the alias rows (`a_dir`, `a_dir_v1`) now gets past `membx.pal`
/// and `taxmiss.pal` and meets `ov1.pal`.
const AWAITING_Q_D2B_OV: [&str; 4] = ["a_ov1", "a_ov2", "a_dir", "a_dir_v1"];

/// Rows of `-cdd_delta`, which LOSAT refuses until row RG-D2b ports its test of a CDD database's
/// headers (the next PR of the row).
const AWAITING_CDD_DELTA: [&str; 12] = [
    "cdd",
    "c_ok",
    "c_ok_v1",
    "c_nowc",
    "c_noaux",
    "c_nofreq",
    "c_freqmagic",
    "c_wcn",
    "c_wcoff",
    "c_obsr28",
    "c_two",
    "c_two_v4",
];

/// `tests/fixtures/blastdbcheck`.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdbcheck")
}

/// The base databases, shared with the blastdb_aliastool fixtures.
fn base_in() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdb_aliastool/in")
}

/// One row of `cases.tsv`.
struct Row {
    id: String,
    dir: String,
    env: Vec<(String, String)>,
    seed: String,
    args: String,
}

fn rows() -> Vec<Row> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let mut part = line.splitn(5, '\t');
            let mut next = || part.next().unwrap_or("").to_string();
            let (id, dir, env, seed, args) = (next(), next(), next(), next(), next());
            let env = if env == "-" {
                Vec::new()
            } else {
                env.split(' ')
                    .map(|pair| {
                        let (name, value) = pair.split_once('=').unwrap();
                        (name.to_string(), value.to_string())
                    })
                    .collect()
            };
            Row {
                id,
                dir,
                env,
                seed,
                args,
            }
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

/// `text` of a `text` operation with `\n`, `\r`, `\t`, `\\` and `<RUN>` replaced.
fn unescape(text: &str, run: &str) -> String {
    let mut out = String::new();
    let mut rest = text;
    while let Some(c) = rest.chars().next() {
        if let Some(after) = rest.strip_prefix("<RUN>") {
            out.push_str(run);
            rest = after;
        } else if c == '\\' && rest.len() > 1 && "nrt\\".contains(&rest[1..2]) {
            out.push(match &rest[1..2] {
                "n" => '\n',
                "r" => '\r',
                "t" => '\t',
                _ => '\\',
            });
            rest = &rest[2..];
        } else {
            out.push(c);
            rest = &rest[c.len_utf8()..];
        }
    }
    out
}

fn hex_bytes(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn set_mode_644(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(path, fs::Permissions::from_mode(0o644)).unwrap();
}

/// Builds `run/in`, `run/home` and `run/cases/<dir>` of `row` by the operations of `setup.tsv`.
fn build(run: &Path, row: &Row) {
    let run_text = run.to_str().unwrap();
    copy_tree(&base_in(), &run.join("in"));
    fs::create_dir_all(run.join("home")).unwrap();
    let dir = run.join("cases").join(&row.dir);
    fs::create_dir_all(&dir).unwrap();
    let files = fixtures().join("files");
    let table = fs::read_to_string(fixtures().join("setup.tsv")).unwrap();
    for line in table
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let part: Vec<&str> = line.split('\t').collect();
        if part[0] != row.dir {
            continue;
        }
        let (op, a) = (part[1], &part[2..]);
        match op {
            "db" => {
                let sub = a.get(1).copied().unwrap_or(".");
                let prefix = a.get(2).copied().unwrap_or("");
                let target = dir.join(sub);
                fs::create_dir_all(&target).unwrap();
                let stem = format!("{}.", a[0]);
                for entry in fs::read_dir(base_in()).unwrap() {
                    let name = entry.unwrap().file_name().to_string_lossy().into_owned();
                    if name.starts_with(&stem) {
                        let to = target.join(format!("{prefix}{name}"));
                        fs::copy(base_in().join(&name), &to).unwrap();
                        set_mode_644(&to);
                    }
                }
            }
            "tax" => {
                for name in ["taxdb.bti", "taxdb.btd"] {
                    fs::copy(files.join("tax").join(a[0]).join(name), dir.join(name)).unwrap();
                }
            }
            "put" => {
                let to = dir.join(a[0]);
                fs::create_dir_all(to.parent().unwrap()).unwrap();
                fs::copy(files.join(a[1]), &to).unwrap();
                set_mode_644(&to);
            }
            "cp" => {
                fs::copy(dir.join(a[0]), dir.join(a[1])).unwrap();
            }
            "rm" => fs::remove_file(dir.join(a[0])).unwrap(),
            "trunc" => {
                let path = dir.join(a[0]);
                let mut bytes = fs::read(&path).unwrap();
                let size = if a[1] == "half" {
                    bytes.len() / 2
                } else if let Some(drop) = a[1].strip_prefix('-') {
                    bytes.len() - drop.parse::<usize>().unwrap()
                } else {
                    a[1].parse().unwrap()
                };
                bytes.truncate(size);
                fs::write(path, bytes).unwrap();
            }
            "poke" => {
                let path = dir.join(a[0]);
                let mut bytes = fs::read(&path).unwrap();
                let at: usize = a[1].parse().unwrap();
                let value = hex_bytes(a[2]);
                bytes[at..at + value.len()].copy_from_slice(&value);
                fs::write(path, bytes).unwrap();
            }
            "text" => {
                let text = unescape(a.get(1).copied().unwrap_or(""), run_text);
                fs::write(dir.join(a[0]), text).unwrap();
            }
            "mkdir" => fs::create_dir_all(dir.join(a[0])).unwrap(),
            "link" => std::os::unix::fs::symlink(a[1], dir.join(a[0])).unwrap(),
            other => panic!("unknown setup operation {other}"),
        }
    }
}

/// Every file below `dir` with its bytes, by its `/`-separated path relative to `dir`
/// (links to directories are not followed).
fn listing(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut found = BTreeMap::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            let kind = entry.file_type().unwrap();
            if kind.is_dir() {
                stack.push(path);
            } else if kind.is_symlink() && path.is_dir() {
                continue;
            } else if let Ok(bytes) = fs::read(&path) {
                let relative = path.strip_prefix(dir).unwrap();
                let parts: Vec<String> = relative
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                found.insert(parts.join("/"), bytes);
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

/// A command of `binary` with no LOSAT, NCBI or BLAST variable, run in `cwd` in UTC with `home` as
/// its home directory.
fn command(binary: &Path, cwd: &Path, home: &Path, row: &Row) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", home)
        .env("HOME", home)
        .env("TZ", "UTC")
        .env("LC_ALL", "C")
        .env("BLAST_USAGE_REPORT", "false");
    for (name, value) in &row.env {
        command.env(name, value);
    }
    command.current_dir(cwd);
    command
}

/// The run of `losat blastdbcheck` for `row` in the fresh run directory `dir`.
fn run(dir: &Path, row: &Row) -> Run {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    let dir = dir.canonicalize().unwrap();
    build(&dir, row);
    let before = listing(&dir);
    let cwd = dir.join("cases").join(&row.dir);
    let arguments: Vec<String> = words(&row.args)
        .into_iter()
        .map(|word| word.replace("<RUN>", dir.to_str().unwrap()))
        .collect();
    let io = dir.with_extension("io");
    fs::create_dir_all(&io).unwrap();
    let mut child = command(
        Path::new(env!("CARGO_BIN_EXE_losat")),
        &cwd,
        &dir.join("home"),
        row,
    )
    .arg("blastdbcheck")
    .args(&arguments)
    .stdin(Stdio::null())
    .stdout(fs::File::create(io.join("out")).unwrap())
    .stderr(fs::File::create(io.join("err")).unwrap())
    .spawn()
    .unwrap();
    let started = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status.code().unwrap_or(-1);
        }
        if started.elapsed() > Duration::from_secs(120) {
            let _ = child.kill();
            let _ = child.wait();
            break -2;
        }
        std::thread::sleep(Duration::from_millis(10));
    };
    let read = |name: &str| {
        let mut bytes = Vec::new();
        fs::File::open(io.join(name))
            .unwrap()
            .read_to_end(&mut bytes)
            .unwrap();
        bytes
    };
    let (stdout, stderr) = (read("out"), read("err"));
    let _ = fs::remove_dir_all(&io);
    let after = listing(&dir);
    let mut files = String::new();
    for (path, bytes) in &after {
        if before.get(path) == Some(bytes) {
            continue;
        }
        let bytes = with_run(bytes, &dir);
        files.push_str(&format!(
            "{path}\t{}\t{:016x}\n",
            bytes.len(),
            fnv1a64(&bytes)
        ));
    }
    for path in before.keys().filter(|path| !after.contains_key(*path)) {
        files.push_str(&format!("{path}\tREMOVED\t-\n"));
    }
    (
        with_run(&stdout, &dir),
        with_run(&stderr, &dir),
        status,
        files,
    )
}

/// NCBI's run of the row `name`, its standard output and error in LOSAT's words
/// (`losat::wording`). A file the run made whose bytes are frozen (`ncbi/<name>.made/<path>`, the
/// `-logfile` rows: NCBI's exception texts) is given by the size and digest of those bytes in
/// LOSAT's words.
fn expected(name: &str) -> Run {
    let file = |extension: &str| fixtures().join("ncbi").join(format!("{name}.{extension}"));
    let err = losat::wording::losat_bytes(&fs::read(file("err")).unwrap()).into_owned();
    let out = losat::wording::losat_output_bytes(&fs::read(file("out")).unwrap()).into_owned();
    let made = file("made");
    let mut files = String::new();
    for line in fs::read_to_string(file("files")).unwrap().lines() {
        let path = line.split('\t').next().unwrap();
        match fs::read(made.join(path)) {
            Ok(bytes) => {
                let bytes = losat::wording::losat_output_bytes(&bytes);
                files.push_str(&format!(
                    "{path}\t{}\t{:016x}\n",
                    bytes.len(),
                    fnv1a64(&bytes)
                ));
            }
            Err(_) => {
                files.push_str(line);
                files.push('\n');
            }
        }
    }
    (
        out,
        err,
        fs::read_to_string(file("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap(),
        files,
    )
}

fn scratch(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("blastdbcheck-{name}-{}", std::process::id()))
}

/// The text of `bytes` for a report, cut to 800 characters.
fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).chars().take(800).collect()
}

/// Every list of rows that LOSAT must refuse.
const REFUSAL_LISTS: [&[&str]; 4] = [
    &LOSAT_REFUSES,
    &AWAITING_ZERO_PAGE,
    &AWAITING_Q_D2B_OV,
    &AWAITING_CDD_DELTA,
];

fn must_refuse(name: &str) -> bool {
    REFUSAL_LISTS.iter().any(|list| list.contains(&name))
}

#[test]
fn blastdbcheck_gives_ncbis_bytes() {
    let mut differ = Vec::new();
    for row in rows() {
        if row.seed != "-" {
            continue;
        }
        let name = row.id.clone();
        let dir = scratch(&name);
        let got = run(&dir, &row);
        if must_refuse(&name) {
            let refused = got.2 == 1
                && got.1.starts_with(b"Error: ")
                && text(&got.1).contains("not supported by LOSAT");
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

/// Every id in the refusal lists names a row of the table, and no row is in two lists.
#[test]
fn refusal_lists_name_rows() {
    let ids: Vec<String> = rows().into_iter().map(|row| row.id).collect();
    let mut seen = Vec::new();
    for name in REFUSAL_LISTS.iter().flat_map(|list| list.iter()) {
        assert!(ids.iter().any(|id| id == name), "no row {name}");
        assert!(!seen.contains(name), "{name} is in two lists");
        seen.push(name);
    }
}
