//! LOSAT's readers of the GI-indexed masks (`makeblastdb -gi_mask`: the alias file's `MASKLIST`,
//! the `.gni`/`.gno`/`.gnd` files) against NCBI BLAST+ 2.17.0 (the official Linux x64 release):
//! `blastdbcmd` (`-info`, `-outfmt %m`, `-mask_sequence_with`), `blastn`, `blastp`, `blastx`,
//! `tblastn` and `tblastx` with `-db_soft_mask` and `-db_hard_mask` (the protein programs' database
//! masks read through the same `seqdb::mask` entry), a plain search of each program family, and
//! `makembindex` (`-db_mask`, `-show_filters`) on the databases that NCBI's makeblastdb built
//! (`tests/seqdb/gi_mask_read`, `tests/seqdb/make_gi_mask_read_fixtures.py`), among them a data
//! file past 64 KiB (`wide`) and `.nog` files cut or edited (`sv4n*`).
//! The exit status, standard output and standard error of every case of `cases.tsv` equal NCBI's
//! (its messages in LOSAT's words, `losat::wording::losat_output_bytes`); a makembindex row also
//! makes or changes the files NCBI's did, with their bytes (`<case>.files`).
//!
//! Each case runs the `losat` binary in the directory the table names, under a root that holds
//! the databases (`<R>` in the frozen files) and the queries, with the environment cleared of
//! `BLASTDB` and the registry as the fixtures ran NCBI's programs. The searches run again with
//! `-num_threads 4` and give the bytes of NCBI's `-num_threads 1` (Owner decision DB-3,
//! 2026-10-13), also where NCBI's blastn with GI-indexed masks stops or crashes in some runs on
//! more than one thread. A round trip builds `m01` (four volumes) and `sv5` (one volume) with
//! `losat makeblastdb` and runs their cases again: the bytes are NCBI's but for the date of the
//! database (`Date:` of `-info`, `Posted date:` of the pairwise report). The rows `raw_*` copy
//! `m01`, `sv5` and `svp` with `makeblastdb -input_type blastdb` (the GI masks become the copy's
//! mask columns, as NCBI's raw database source registers them) and read the copies back, with
//! the copies' dates and the log's clock masked.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp,
//! c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:117-139 (`LOSAT/src/seqdb/gi_mask.rs`).

// The frozen files are those of the Linux x64 release; the root of a message is the directory as
// the child process sees it.
#![cfg(target_os = "linux")]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::SystemTime;

use crate::makembindex_fixtures::fnv1a64;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/gi_mask_read")
}

struct Case {
    name: String,
    dbs: Vec<String>,
    cwd: String,
    env: Option<(String, String)>,
    words: Vec<String>,
}

fn cases() -> Vec<Case> {
    let table = std::fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .map(|line| {
            let mut fields = line.split('\t').map(str::to_string);
            let name = fields.next().unwrap();
            let dbs = fields
                .next()
                .unwrap()
                .split(',')
                .map(str::to_string)
                .collect();
            let cwd = fields.next().unwrap();
            let env = fields.next().unwrap();
            let env = (env != "-").then(|| {
                let (key, value) = env.split_once('=').unwrap();
                (key.to_string(), value.to_string())
            });
            Case {
                name,
                dbs,
                cwd,
                env,
                words: fields.collect(),
            }
        })
        .collect()
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(label: &str) -> Self {
        #[allow(clippy::disallowed_methods)] // a Linux-only test
        let base = std::env::temp_dir();
        #[allow(clippy::disallowed_methods)]
        let pid = std::process::id();
        let dir = std::fs::canonicalize(base)
            .unwrap()
            .join(format!("losat-gimask-read-{label}-{pid}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The text files of a database (alias and metadata), which hold the root as `<R>`.
fn is_text(name: &str) -> bool {
    [".nal", ".pal", ".njs", ".pjs"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

/// Copies the frozen database `db` to `root/db`, the root put back in its text files.
fn put_database(root: &Path, db: &str) {
    let from = fixtures().join("db").join(db);
    let to = root.join(db);
    if to.exists() {
        return;
    }
    std::fs::create_dir_all(&to).unwrap();
    for entry in std::fs::read_dir(&from).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name().to_string_lossy().into_owned();
        let bytes = std::fs::read(entry.path()).unwrap();
        let bytes = match is_text(&name) {
            true => String::from_utf8(bytes)
                .unwrap()
                .replace("<R>", root.to_str().unwrap())
                .into_bytes(),
            false => bytes,
        };
        std::fs::write(to.join(&name), bytes).unwrap();
    }
}

/// The root of the runs: the queries, and the databases are put there as the cases need them.
fn root(label: &str) -> TempDir {
    let root = TempDir::new(label);
    for query in ["q.fa", "pq.fa"] {
        std::fs::copy(fixtures().join("in").join(query), root.0.join(query)).unwrap();
    }
    root
}

/// Each file below `dir`: its `/`-separated path relative to it, its size and time of change.
fn file_states(dir: &Path) -> BTreeMap<String, (u64, SystemTime)> {
    let mut found = BTreeMap::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in std::fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            let meta = std::fs::metadata(&path).unwrap();
            if meta.is_dir() {
                stack.push(path);
                continue;
            }
            let relative = path.strip_prefix(dir).unwrap();
            let parts: Vec<String> = relative
                .components()
                .map(|part| part.as_os_str().to_string_lossy().into_owned())
                .collect();
            found.insert(parts.join("/"), (meta.len(), meta.modified().unwrap()));
        }
    }
    found
}

/// The files below `root` made or changed since `before`: `path<TAB>size<TAB>FNV-1a 64` lines.
fn made_files(root: &Path, before: &BTreeMap<String, (u64, SystemTime)>) -> String {
    let mut files = String::new();
    for (path, state) in file_states(root) {
        if before.get(&path) == Some(&state) {
            continue;
        }
        let bytes = std::fs::read(root.join(&path)).unwrap();
        files.push_str(&format!(
            "{path}\t{}\t{:016x}\n",
            bytes.len(),
            fnv1a64(&bytes)
        ));
    }
    files
}

/// What a case gave: exit status, standard output, standard error, and for a makembindex row the
/// files it made or changed (`made_files`; empty for the other rows).
type Output = (i32, Vec<u8>, Vec<u8>, String);

fn is_makembindex(case: &Case) -> bool {
    case.words[0] == "makembindex"
}

/// Runs `case` under `root` with `extra` words after its own: exit status, standard output and
/// error, the root as `<R>`, and the files it made (`Output`).
fn run(root: &Path, case: &Case, extra: &[&str]) -> Output {
    let r = root.to_str().unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .current_dir(root.join(&case.cwd))
        // `$HOME/.ncbirc` is read like the one of the current directory.
        .env("HOME", root)
        .args(case.words.iter().map(|word| word.replace("{R}", r)))
        .args(extra);
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy().into_owned();
        // `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
        // (`cargo test --all-features`) writes its worker log there.
        let losat = name.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG";
        if losat || name.starts_with("NCBI_CONFIG") || name == "BLASTDB" {
            command.env_remove(name);
        }
    }
    if let Some((key, value)) = &case.env {
        command.env(key, value.replace("{R}", r));
    }
    let before = is_makembindex(case).then(|| file_states(root));
    let output = command.output().expect("run losat");
    let files = match &before {
        Some(before) => made_files(root, before),
        None => String::new(),
    };
    let unroot = |bytes: Vec<u8>| {
        String::from_utf8_lossy(&bytes)
            .replace(r, "<R>")
            .into_bytes()
    };
    (
        output.status.code().unwrap_or(-1),
        unroot(output.stdout),
        unroot(output.stderr),
        files,
    )
}

/// NCBI's frozen run of `case`, its messages in LOSAT's words.
fn expected(case: &Case) -> Output {
    let name = &case.name;
    let dir = fixtures().join("expected");
    let read = |ext: &str| std::fs::read(dir.join(format!("{name}.{ext}"))).unwrap();
    let rc = String::from_utf8(read("rc"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    let mapped = |bytes: Vec<u8>| losat::wording::losat_output_bytes(&bytes).into_owned();
    let files = match is_makembindex(case) {
        true => String::from_utf8(read("files")).unwrap(),
        false => String::new(),
    };
    (rc, mapped(read("stdout")), mapped(read("stderr")), files)
}

/// The date of the database masked: the text after `Date: ` (up to a tab or the end of the line)
/// and after the pairwise report's `    Posted date:  `.
fn without_dates(bytes: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(bytes);
    text.split_inclusive('\n')
        .map(|line| {
            for start in ["Date: ", "    Posted date:  "] {
                if let Some(rest) = line.strip_prefix(start) {
                    let end = rest.find(['\t', '\n']).unwrap_or(rest.len());
                    return format!("{start}<D>{}", &rest[end..]);
                }
            }
            line.to_string()
        })
        .collect::<String>()
        .into_bytes()
}

/// The first line where `got` and `want` differ: its number, LOSAT's line and NCBI's.
fn first_difference(got: &[u8], want: &[u8]) -> String {
    let got = String::from_utf8_lossy(got);
    let want = String::from_utf8_lossy(want);
    let mut got_lines = got.split_inclusive('\n');
    let mut want_lines = want.split_inclusive('\n');
    for line in 1.. {
        match (got_lines.next(), want_lines.next()) {
            (None, None) => break,
            (a, b) if a != b => return format!("line {line}: {a:?} (NCBI {b:?})"),
            _ => {}
        }
    }
    "equal".to_string()
}

fn describe(got: &Output, want: &Output) -> String {
    format!(
        "rc {} (NCBI {}), stdout {}, stderr {}, files {}",
        got.0,
        want.0,
        first_difference(&got.1, &want.1),
        first_difference(&got.2, &want.2),
        first_difference(got.3.as_bytes(), want.3.as_bytes())
    )
}

fn is_search(case: &Case) -> bool {
    !["blastdbcmd", "makeblastdb", "makembindex"].contains(&case.words[0].as_str())
}

/// The rows `raw_*` copy a database with `makeblastdb -input_type blastdb` and read the copy back:
/// the copy's date (`without_dates`), the log's clock (`current time: `) and the time the copy
/// took (`added N sequences in S seconds.`) are masked.
fn comparable(case: &Case, out: Output) -> Output {
    if !case.name.starts_with("raw_") {
        return out;
    }
    let text = String::from_utf8_lossy(&without_dates(&out.1)).into_owned();
    let stdout = text
        .split_inclusive('\n')
        .map(|line| {
            if let Some(at) = line.find("current time: ") {
                return format!("{}current time: <T>\n", &line[..at]);
            }
            match line.find(" sequences in ") {
                Some(at) if line.ends_with(" seconds.\n") => {
                    format!("{} sequences in <S> seconds.\n", &line[..at])
                }
                _ => line.to_string(),
            }
        })
        .collect::<String>();
    (out.0, stdout.into_bytes(), out.2, out.3)
}

/// Every case of the table, and the searches again with `-num_threads 4`, whose reference is
/// NCBI's `-num_threads 1` output (DB-3).
#[test]
fn readers_give_ncbi_s_output_on_ncbi_s_databases() {
    let root = root("ncbi");
    let mut problems = Vec::new();
    for case in cases() {
        for db in &case.dbs {
            put_database(&root.0, db);
        }
        let want = comparable(&case, expected(&case));
        let got = comparable(&case, run(&root.0, &case, &[]));
        if got != want {
            problems.push(format!("{}: {}", case.name, describe(&got, &want)));
        }
        if is_search(&case) {
            let got = comparable(&case, run(&root.0, &case, &["-num_threads", "4"]));
            if got != want {
                problems.push(format!("{} n4: {}", case.name, describe(&got, &want)));
            }
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

/// The databases `m01` and `sv5` built by `losat makeblastdb` from the same inputs read as NCBI's
/// do, but for their dates.
#[test]
fn readers_give_ncbi_s_output_on_losat_s_databases() {
    let root = root("losat");
    let inputs = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/gi_mask/in");
    let dbs = std::fs::read_to_string(fixtures().join("dbs.tsv")).unwrap();
    for line in dbs.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if !["m01", "sv5"].contains(&fields[0]) {
            continue;
        }
        let dir = root.0.join(fields[0]);
        std::fs::create_dir_all(&dir).unwrap();
        for input in fields[2].split(' ') {
            std::fs::copy(inputs.join(input), dir.join(input)).unwrap();
        }
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .current_dir(&dir)
            .env("HOME", &dir)
            .env_remove("BLASTDB")
            .arg("makeblastdb")
            .args(&fields[3..])
            .output()
            .unwrap();
        assert!(output.status.success(), "{}: {output:?}", fields[0]);
        for input in fields[2].split(' ') {
            std::fs::remove_file(dir.join(input)).unwrap();
        }
    }
    let mut problems = Vec::new();
    let mut ran = 0;
    for case in cases() {
        let own = case.dbs.iter().all(|db| db == "m01" || db == "sv5");
        // `copy2` is the copy of `-in 'db db'`, made by a row skipped here.
        let alias = case
            .words
            .iter()
            .any(|word| word == "al" || word == "db db" || word == "copy2");
        if !own || alias {
            continue;
        }
        ran += 1;
        let want = comparable(&case, expected(&case));
        let want = (want.0, without_dates(&want.1), want.2, want.3);
        let got = comparable(&case, run(&root.0, &case, &[]));
        let got = (got.0, without_dates(&got.1), got.2, got.3);
        if got != want {
            problems.push(format!("{}: {}", case.name, describe(&got, &want)));
        }
    }
    assert!(ran > 40, "{ran}");
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
