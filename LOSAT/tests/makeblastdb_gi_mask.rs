//! `losat makeblastdb -gi_mask` and `-gi_mask_name` against makeblastdb 2.17.0 (the official
//! Linux x64 release): the exit status, the standard output and error, the files of the run
//! directory with their sizes, and the bytes of the files that `-gi_mask` writes or changes (the
//! GI-mask files, the alias file, the metadata file, the mask column files, the log file), for
//! every case of `tests/seqdb/gi_mask/cases.tsv` (`tests/seqdb/make_gi_mask_fixtures.py`).
//!
//! Each case runs the `losat` binary in a directory of its own (a case whose setup starts with
//! `^` in the directory of the case before it), with its inputs copied in first and the names on
//! the command line relative to it, as the fixtures ran makeblastdb. The run directory is `<C>`
//! in the messages and in the alias, metadata and log files; the time of "Building a new DB" is
//! `<T>` and the seconds of "Adding sequences" `<S>`. The bytes that hold the run's date are
//! masked: the date of a `.gmi` or `.gni` file (its length is compared), the creation time of a
//! mask column file (`.naa`, `.paa`), the line "Alias file created" of an alias file, and
//! "last-updated" of a metadata file. Its "bytes-total" is compared as its distance from the sum
//! of the sizes of the files it lists, each time it lists one, in the run directory, the metadata
//! file itself at the size it had before the run (makeblastdb lists the files before it writes the
//! new one, so an old one that the run did not remove counts as it was): 0 for makeblastdb's (the
//! script checks that when it freezes them). A metadata file that the case's run did not write (an
//! earlier case of its directory did) has it masked.
//!
//! The cases of the column `losat` other than `same`: `refuse` (two algorithms that would write
//! one GI-mask data file at the same time, which LOSAT refuses at the record that would open it:
//! exit status 1, and what an error of NCBI's `SetMaskData` leaves, `REFUSED_LEFTOVERS`) and
//! `parser` (an argument rule that LOSAT's parser checks: exit status 2 and no file, where
//! makeblastdb prints its usage; LOSAT's message holds the sentence of makeblastdb's).
//!
//! The expectations of a single-volume database built with `-gi_mask` are NCBI's run with the one
//! change of PD-LOSAT-NCBI-DEFECTS approved exception 22 (`single_volume_exception` of the
//! script): the alias file's `DBLIST` names the volume by its file name, and the run ends as a
//! successful one, with the metadata file of NCBI's `blastdbcmd -metadata` for that database.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp (the files),
//! c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:338-342, 405-438 (their close and the
//! alias file), c++/src/app/blastdb/makeblastdb.cpp:1205, 1231-1241 (`-gi_mask_name`).

// The frozen files are those of the Linux x64 release: the LMDB files follow the host's page
// size, and the run directory of a message is the directory as the child process sees it.
#![cfg(target_os = "linux")]

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/gi_mask")
}

struct Case {
    name: String,
    losat: String,
    setup: String,
    words: Vec<String>,
}

fn cases(dir: &Path) -> Vec<Case> {
    let table = std::fs::read_to_string(dir.join("cases.tsv")).unwrap();
    table
        .lines()
        .map(|line| {
            let mut fields = line.split('\t').map(str::to_string);
            Case {
                name: fields.next().unwrap(),
                losat: fields.next().unwrap(),
                setup: fields.next().unwrap(),
                words: fields.collect(),
            }
        })
        .collect()
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(label: &str) -> Self {
        // The run directories as `getcwd` names them in the messages (`<C>`).
        let dir = std::fs::canonicalize(std::env::temp_dir())
            .unwrap()
            .join(format!("losat-mkdb-gimask-{label}-{}", std::process::id()));
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

/// The setup of a case in `case_dir`: `D/` makes a directory, `F` copies the input `F`, `F>G`
/// copies it to `G`, `F^N>G` writes its first `N` FASTA records to `G`, `F=>G` copies it to `G`
/// with `{C}` replaced by `case_dir` (a file of a database that an earlier run left, not an
/// input: it is listed and compared as the run leaves it). Returns the paths of the inputs put
/// there.
fn setup(case_dir: &Path, tokens: &[&str], inputs: &Path) -> Vec<String> {
    let mut put = Vec::new();
    for token in tokens {
        if token.ends_with('/') {
            std::fs::create_dir_all(case_dir.join(token)).unwrap();
            continue;
        }
        if let Some((source, target)) = token.split_once("=>") {
            let text = std::fs::read_to_string(inputs.join(source)).unwrap();
            let text = text.replace("{C}", case_dir.to_str().unwrap());
            std::fs::write(case_dir.join(target), text).unwrap();
            continue;
        }
        let (source, target) = token.split_once('>').unwrap_or((token, token));
        match source.split_once('^') {
            Some((source, count)) => {
                let text = std::fs::read_to_string(inputs.join(source)).unwrap();
                let records: String = text
                    .split('>')
                    .skip(1)
                    .take(count.parse().unwrap())
                    .map(|record| format!(">{record}"))
                    .collect();
                std::fs::write(case_dir.join(target), records).unwrap();
            }
            None => {
                std::fs::copy(inputs.join(source), case_dir.join(target)).unwrap();
            }
        }
        put.push(target.to_string());
    }
    put
}

/// The run directory, the time and the seconds of a message or a text file masked.
fn masked(text: &str, c: &str) -> String {
    let mut out = text.replace(c, "<C>");
    for (start, end, fill) in [
        ("current time: ", "\n", "<T>"),
        (" sequences in ", " seconds.", "<S>"),
        ("# Alias file created: ", "\n", "<D>"),
        ("\"last-updated\": \"", "\"", "<D>"),
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

/// How the "bytes-total" of a metadata file is compared (`bytes_total_distance`).
#[derive(Clone, Copy)]
enum BytesTotal<'a> {
    /// A file of LOSAT's run under this directory, with the sizes of the run's metadata files
    /// before the run (`metadata_sizes`).
    Run(&'a Path, &'a HashMap<String, u64>),
    /// A file of makeblastdb's run, whose value is the sum (distance 0).
    Frozen,
    /// A file that an earlier case of the directory wrote: masked.
    Masked,
}

/// The "bytes-total" of the metadata file `rel` (`text`) replaced by its distance from the sum of
/// the sizes of the files that it lists (each time it lists one), read in the directory of the
/// metadata file of LOSAT's run (`total`); the metadata file itself counts at its size before the
/// run.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1610-1703 (the files
// that exist when the metadata is made, the old metadata file among them, are listed and summed)
fn bytes_total_distance(rel: &str, text: &str, total: BytesTotal<'_>) -> String {
    let key = "\"bytes-total\": ";
    let Some(at) = text.find(key).map(|at| at + key.len()) else {
        return text.to_string();
    };
    let end = at + text[at..].find(',').unwrap_or(0);
    let distance = match total {
        BytesTotal::Frozen => "0".to_string(),
        BytesTotal::Masked => "N".to_string(),
        BytesTotal::Run(run, before) => {
            let path = run.join(rel);
            let dir = path.parent().unwrap();
            let own = path.file_name().unwrap().to_string_lossy().into_owned();
            let listed = text
                .split_once("\"files\": [")
                .map_or("", |(_, rest)| rest.split(']').next().unwrap_or(""));
            let mut sum = 0i128;
            for name in listed
                .split(',')
                .map(|item| item.trim().trim_matches('"'))
                .filter(|name| !name.is_empty())
            {
                let size = match name == own {
                    true => before.get(rel).copied(),
                    false => std::fs::metadata(dir.join(name))
                        .ok()
                        .map(|meta| meta.len()),
                };
                match size {
                    Some(size) => sum += i128::from(size),
                    None => return format!("{text}\n(listed file {name} is missing)"),
                }
            }
            match text[at..end].trim().parse::<i128>() {
                Ok(total) => (total - sum).to_string(),
                Err(_) => format!("unreadable {:?}", &text[at..end]),
            }
        }
    };
    format!("{}<{distance}>{}", &text[..at], &text[end..])
}

/// `CBlastDbBlob::x_ReadVarInt` (see `LOSAT/src/writedb/gi_mask_tests.rs`).
fn read_var_int(bytes: &[u8], at: &mut usize) -> usize {
    let mut value = 0usize;
    loop {
        let byte = bytes[*at];
        *at += 1;
        if byte & 0x80 != 0 {
            value = (value << 7) | usize::from(byte & 0x7F);
        } else {
            return (value << 6) | usize::from(byte & 0x3F);
        }
    }
}

/// A file of a run with the bytes that hold the run's date or directory masked; `total` says how
/// the "bytes-total" of a metadata file is compared (`BytesTotal::Run` with the run directory).
fn masked_file(rel: &str, bytes: Vec<u8>, c: &str, total: BytesTotal<'_>) -> Vec<u8> {
    if rel.ends_with(".gmi") || rel.ends_with(".gni") {
        let mut bytes = bytes;
        let mut at = 32;
        let desc = read_var_int(&bytes, &mut at);
        at += desc;
        let date = read_var_int(&bytes, &mut at);
        bytes[at..at + date].fill(b'#');
        return bytes;
    }
    if rel.ends_with(".naa") || rel.ends_with(".paa") {
        // The column's creation time (`CTime::AsString`, `MM/DD/YYYY hh:mm:ss`).
        let mut bytes = bytes;
        let pattern = b"00/00/0000 00:00:00";
        let mut at = 0;
        while at + pattern.len() <= bytes.len() {
            let window = &bytes[at..at + pattern.len()];
            let date = window.iter().zip(pattern).all(|(&b, &p)| match p {
                b'0' => b.is_ascii_digit(),
                _ => b == p,
            });
            if date {
                bytes[at..at + pattern.len()].fill(b'#');
                at += pattern.len();
            } else {
                at += 1;
            }
        }
        return bytes;
    }
    let text = [".nal", ".pal", ".njs", ".pjs"]
        .iter()
        .any(|ext| rel.ends_with(ext))
        || rel == "log.txt";
    let metadata = rel.ends_with(".njs") || rel.ends_with(".pjs");
    match (text, metadata) {
        (true, true) => {
            let text = bytes_total_distance(rel, &String::from_utf8_lossy(&bytes), total);
            masked(&text, c).into_bytes()
        }
        (true, false) => masked(&String::from_utf8_lossy(&bytes), c).into_bytes(),
        (false, _) => bytes,
    }
}

/// The files under `dir` (not the inputs `put`), as `path<TAB>size`, `*` for a file whose size
/// follows the date or the directory of the run; sorted.
fn listing(dir: &Path, put: &[String]) -> Vec<String> {
    let mut files = Vec::new();
    let mut stack = vec![PathBuf::new()];
    while let Some(rel) = stack.pop() {
        for entry in std::fs::read_dir(dir.join(&rel)).unwrap() {
            let entry = entry.unwrap();
            let path = rel.join(entry.file_name());
            let name = path.to_string_lossy().into_owned();
            let meta = std::fs::metadata(entry.path()).unwrap();
            if meta.is_dir() {
                stack.push(path);
            } else if !put.contains(&name) {
                let sized = [".nin", ".pin", ".njs", ".pjs", ".nal", ".pal"]
                    .iter()
                    .any(|ext| name.ends_with(ext))
                    || name == "log.txt";
                let size = if sized {
                    "*".to_string()
                } else {
                    meta.len().to_string()
                };
                files.push(format!("{name}\t{size}"));
            }
        }
    }
    files.sort();
    files
}

/// Every file under `dir`, as paths relative to it.
fn files_under(dir: &Path) -> Vec<String> {
    listing(dir, &[])
        .into_iter()
        .map(|line| line.split('\t').next().unwrap().to_string())
        .collect()
}

/// The sizes of the metadata files under `dir` (relative path), taken before a run.
fn metadata_sizes(dir: &Path) -> HashMap<String, u64> {
    files_under(dir)
        .into_iter()
        .filter(|rel| rel.ends_with(".njs") || rel.ends_with(".pjs"))
        .map(|rel| {
            let size = std::fs::metadata(dir.join(&rel)).unwrap().len();
            (rel, size)
        })
        .collect()
}

fn run(case_dir: &Path, words: &[String]) -> (i32, String, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .current_dir(case_dir)
        // `$HOME/.ncbirc` is read like the one of the current directory.
        .env("HOME", case_dir)
        .arg("makeblastdb")
        .args(words);
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy().into_owned();
        if name.starts_with("LOSAT_")
            || name.starts_with("NCBI_CONFIG")
            || name == "BLASTDB"
            || name == "BLASTDB_LMDB_MAP_SIZE"
            || name == "MAX_LMDB_TXN_ENTRY"
            || name.starts_with("LMDB_")
        {
            command.env_remove(name);
        }
    }
    let output = command.output().expect("run losat makeblastdb");
    let c = case_dir.to_str().unwrap();
    (
        output.status.code().unwrap_or(-1),
        masked(&String::from_utf8_lossy(&output.stdout), c),
        masked(&String::from_utf8_lossy(&output.stderr), c),
    )
}

fn first_difference(a: &[u8], b: &[u8]) -> String {
    let at = a.iter().zip(b).position(|(x, y)| x != y);
    format!(
        "{} bytes vs {}, first difference at {:?}",
        a.len(),
        b.len(),
        at.or((a.len() != b.len()).then_some(a.len().min(b.len())))
    )
}

/// The differences of a case of `same` from makeblastdb's run.
fn compare(
    case: &Case,
    case_dir: &Path,
    put: &[String],
    got: (i32, String, String),
    started: std::time::SystemTime,
    before: &HashMap<String, u64>,
) -> Vec<String> {
    let dir = fixtures();
    // makeblastdb's messages in LOSAT's words (`losat::wording`, the one mapping of the freeze).
    let expected = |ext: &str| {
        let frozen =
            std::fs::read(dir.join("expected").join(format!("{}.{ext}", case.name))).unwrap();
        String::from_utf8_lossy(&losat::wording::losat_output_bytes(&frozen)).into_owned()
    };
    let mut out = Vec::new();
    let want_rc: i32 = expected("rc").trim().parse().unwrap();
    for (what, got, want) in [
        ("status", got.0.to_string(), want_rc.to_string()),
        ("stdout", got.1, expected("stdout")),
        ("stderr", got.2, expected("stderr")),
    ] {
        if got != want {
            out.push(format!("{what}:\n    got  {got:?}\n    want {want:?}"));
        }
    }
    let files = listing(case_dir, put);
    let mut want_files: Vec<String> = expected("files").lines().map(str::to_string).collect();
    want_files.sort();
    if files != want_files {
        let missing: Vec<_> = want_files.iter().filter(|f| !files.contains(f)).collect();
        let extra: Vec<_> = files.iter().filter(|f| !want_files.contains(f)).collect();
        out.push(format!(
            "files: missing {missing:?}, not in makeblastdb's {extra:?}"
        ));
    }
    // A case that leaves none of the files kept has no directory there.
    let kept = dir.join("db").join(&case.name);
    let c = case_dir.to_str().unwrap();
    let kept_files = match kept.is_dir() {
        true => files_under(&kept),
        false => Vec::new(),
    };
    for rel in kept_files {
        let mut frozen = std::fs::read(kept.join(&rel)).unwrap();
        if rel == "log.txt" {
            frozen = losat::wording::losat_output_bytes(&frozen).into_owned();
        }
        // A metadata file of an earlier case of the directory lists files that this run changed
        // (`started` is the stamp's time, which the file system gives).
        let fresh = std::fs::metadata(case_dir.join(&rel))
            .and_then(|meta| meta.modified())
            .is_ok_and(|modified| modified >= started);
        let (want_total, got_total) = match fresh {
            true => (BytesTotal::Frozen, BytesTotal::Run(case_dir, before)),
            false => (BytesTotal::Masked, BytesTotal::Masked),
        };
        let want = masked_file(&rel, frozen, c, want_total);
        match std::fs::read(case_dir.join(&rel)) {
            Ok(bytes) => {
                let got = masked_file(&rel, bytes, c, got_total);
                if got != want {
                    let text = [".nal", ".njs", ".pal", ".pjs", ".txt"]
                        .iter()
                        .any(|ext| rel.ends_with(ext));
                    out.push(match text {
                        true => format!(
                            "{rel}:\n    got  {:?}\n    want {:?}",
                            String::from_utf8_lossy(&got),
                            String::from_utf8_lossy(&want)
                        ),
                        false => format!("{rel}: {}", first_difference(&got, &want)),
                    });
                }
            }
            Err(error) => out.push(format!("{rel}: {error}")),
        }
    }
    out
}

#[test]
fn gi_mask_builds_end_as_makeblastdb_s() {
    let dir = fixtures();
    let inputs = dir.join("in");
    let tmp = TempDir::new("cases");
    let mut failures = Vec::new();
    let mut case_dir = PathBuf::new();
    let mut put = Vec::new();
    let mut counts = [0usize; 3];
    for case in cases(&dir) {
        let mut tokens: Vec<&str> = case.setup.split_whitespace().collect();
        if tokens.first() == Some(&"^") {
            tokens.remove(0);
        } else {
            case_dir = tmp.0.join(&case.name);
            std::fs::create_dir_all(&case_dir).unwrap();
            put.clear();
        }
        put.extend(setup(&case_dir, &tokens, &inputs));
        let c = case_dir.to_str().unwrap().to_string();
        let words: Vec<String> = case.words.iter().map(|w| w.replace("{C}", &c)).collect();
        // The time of the run as the file system stamps it (a stamp file written just before):
        // a file of the directory older than that is an earlier case's.
        if case.setup.starts_with('^') {
            std::thread::sleep(std::time::Duration::from_millis(20));
        }
        let stamp = tmp.0.join("stamp");
        std::fs::write(&stamp, case.name.as_bytes()).unwrap();
        let started = std::fs::metadata(&stamp).unwrap().modified().unwrap();
        let before = metadata_sizes(&case_dir);
        let got = run(&case_dir, &words);
        let problems = match case.losat.as_str() {
            "same" => {
                counts[0] += 1;
                compare(&case, &case_dir, &put, got, started, &before)
            }
            "refuse" => {
                counts[1] += 1;
                let mut problems = Vec::new();
                let refusal = "losat makeblastdb: -gi_mask: two masking algorithms would write \
                               one GI-mask data file at the same time";
                if got.0 != 1 || !got.2.starts_with(refusal) || got.2.lines().count() != 1 {
                    problems.push(format!("not refused: {got:?}"));
                }
                let files: Vec<String> = listing(&case_dir, &put)
                    .into_iter()
                    .map(|line| line.split('\t').next().unwrap().to_string())
                    .collect();
                if files != REFUSED_LEFTOVERS {
                    problems.push(format!("files left: {files:?}"));
                }
                problems
            }
            "parser" => {
                counts[2] += 1;
                let files = listing(&case_dir, &put);
                let mut problems = Vec::new();
                if got.0 != 2 || !got.1.is_empty() || !files.is_empty() {
                    problems.push(format!("not a parser error: {got:?}, files {files:?}"));
                }
                // makeblastdb's sentence of a dependency, in LOSAT's message.
                let frozen =
                    std::fs::read(dir.join("expected").join(format!("{}.stderr", case.name)))
                        .unwrap();
                let frozen = String::from_utf8_lossy(&losat::wording::losat_output_bytes(&frozen))
                    .into_owned();
                let sentence = frozen
                    .lines()
                    .find_map(|line| line.split_once("(CArgException::eConstraint) "))
                    .map(|(_, sentence)| sentence);
                if let Some(sentence) = sentence {
                    if !got.2.contains(sentence) {
                        problems.push(format!("not makeblastdb's sentence {sentence:?}: {got:?}"));
                    }
                }
                problems
            }
            other => vec![format!("unknown column {other:?}")],
        };
        if !problems.is_empty() {
            failures.push(format!(
                "{} {:?}:\n  {}",
                case.name,
                case.words,
                problems.join("\n  ")
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} cases differ:\n{}",
        failures.len(),
        counts.iter().sum::<usize>(),
        failures.join("\n")
    );
    // Every kind of case is there (the table was read whole).
    assert!(
        counts[0] >= 75 && counts[1] >= 3 && counts[2] >= 9,
        "{counts:?}"
    );
}

/// What a refused case leaves (`refuse`): the build ends at the refused record as an error of
/// NCBI's `SetMaskData` ends it (case `x03gi`): the files that `ListFiles` names are removed, and
/// the GI-mask files of the algorithm that wrote first (closed with the build), the alias file
/// and the files of the version 5 LMDB step that `ListFiles` does not name stay.
const REFUSED_LEFTOVERS: [&str; 11] = [
    "db.nal", "db.nos", "db.not", "db.ntf", "db.nto", "x.gmd", "x.gmi", "x.gmo", "x.gnd", "x.gni",
    "x.gno",
];

/// `losat <command>` in `dir` (`command`: the program and its words, separated by single spaces),
/// with the environment of `run`: exit status, standard output and error.
fn losat(dir: &Path, command: &str) -> (i32, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(dir)
        .env("HOME", dir)
        .env_remove("BLASTDB")
        .args(command.split(' '))
        .output()
        .expect("run losat");
    (
        output.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// The database of one volume that `losat makeblastdb -gi_mask` writes opens (its alias file
/// names the volume by its file name, PD-LOSAT-NCBI-DEFECTS approved exception 22), and LOSAT's
/// programs that read sequences stop on its `MASKLIST` with LOSAT's refusal, before they print
/// anything: they do not read GI masks yet (line NP item M4c-2).
#[test]
fn readers_refuse_the_gi_masks_of_a_single_volume_database() {
    let inputs = fixtures().join("in");
    let tmp = TempDir::new("readers");
    for name in ["nucl.fa", "dust_nucl.asnb"] {
        std::fs::copy(inputs.join(name), tmp.0.join(name)).unwrap();
    }
    let build = losat(
        &tmp.0,
        "makeblastdb -in nucl.fa -dbtype nucl -parse_seqids -out db -title t \
         -mask_data dust_nucl.asnb -gi_mask -gi_mask_name gm1",
    );
    assert_eq!((build.0, build.2.as_str()), (0, ""), "{build:?}");
    let alias = std::fs::read_to_string(tmp.0.join("db.nal")).unwrap();
    assert!(alias.ends_with("\nDBLIST db\nMASKLIST gm1 \n"), "{alias:?}");
    assert!(tmp.0.join("db.njs").is_file());
    let refusal = "unsupported database: alias file key MASKLIST is not supported";
    for command in [
        "blastdbcmd -db db -info",
        "blastdbcmd -db db -entry all -outfmt %m0",
        "blastn -query nucl.fa -db db -db_soft_mask 0 -outfmt 6",
    ] {
        let got = losat(&tmp.0, command);
        assert!(
            got.0 == 1 && got.1.is_empty() && got.2.contains(refusal),
            "{command}: {got:?}"
        );
    }
}
