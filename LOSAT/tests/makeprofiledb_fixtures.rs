//! `losat makeprofiledb` against NCBI BLAST+ 2.17.0 makeprofiledb (`-dbtype rps`): each row of
//! `tests/fixtures/makeprofiledb/cases.tsv` runs in its own copy of `in/` and gives NCBI's standard
//! output, standard error (in LOSAT's words, `losat::wording`) and exit status
//! (`ncbi/<row>.{out,err,rc}`) and, for each file the run wrote, its size and FNV-1a digest
//! (`ncbi/<row>.files`; the date of a `.pin` or `.pal` file replaced by `<DATE>`, its size `*`),
//! frozen from the oracle of row PS10a (2026-11-21; README.md of the directory). NCBI is not run
//! here. Then `losat rpsblast` on databases that `losat makeprofiledb` built gives NCBI rpsblast's
//! bytes on NCBI makeprofiledb's databases (`rps/`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// `tests/fixtures/makeprofiledb`.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/makeprofiledb")
}

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

/// The bytes of a file that holds the build's date with the date replaced by `<DATE>`: in a
/// `.pin` file the date's 4-byte length and its bytes with their padding (`Nov 21, 2026  5:07
/// AM`), in a `.pal` file the `MM/DD/YYYY hh:mm:ss` of its first line; `None` for another file.
fn masked(name: &str, bytes: &[u8]) -> Option<Vec<u8>> {
    const MONTHS: [&[u8]; 12] = [
        b"Jan ", b"Feb ", b"Mar ", b"Apr ", b"May ", b"Jun ", b"Jul ", b"Aug ", b"Sep ", b"Oct ",
        b"Nov ", b"Dec ",
    ];
    if name.ends_with(".pin") {
        // The date field: its length (with the padding), then a month's name, the time and AM or
        // PM, then NULs.
        let field_length =
            |at: usize| u32::from_be_bytes(bytes[at - 4..at].try_into().unwrap()) as usize;
        let start = (4..bytes.len()).find(|&at| {
            let length = field_length(at);
            MONTHS.iter().any(|month| bytes[at..].starts_with(month))
                && (16..=40).contains(&length)
                && at + length <= bytes.len()
                && {
                    let text = &bytes[at..at + length];
                    let end = text.iter().rposition(|&b| b != 0).map_or(0, |i| i + 1);
                    text[..end].ends_with(b"AM") || text[..end].ends_with(b"PM")
                }
        })?;
        let length = field_length(start);
        let mut out = bytes[..start - 4].to_vec();
        out.extend_from_slice(b"<DATE>");
        out.extend_from_slice(&bytes[start + length..]);
        return Some(out);
    }
    if name.ends_with(".pal") {
        let marker = b"# Alias file created ";
        let at = bytes.windows(marker.len()).position(|w| w == marker)? + marker.len();
        let mut out = bytes[..at].to_vec();
        out.extend_from_slice(b"<DATE>");
        out.extend_from_slice(&bytes[at + "MM/DD/YYYY hh:mm:ss".len()..]);
        return Some(out);
    }
    None
}

/// A command of `binary` with no LOSAT, NCBI or BLAST variable, run in `dir`.
fn command(binary: &Path, dir: &Path) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    command.env("NCBI_CONFIG_PATH", dir).current_dir(dir);
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

/// The files below `dir` outside `in/`, as `/`-separated paths relative to it.
fn files_below(dir: &Path) -> Vec<String> {
    let mut found = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                if path != dir.join("in") {
                    stack.push(path);
                }
            } else {
                let relative = path.strip_prefix(dir).unwrap();
                let parts: Vec<String> = relative
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                found.push(parts.join("/"));
            }
        }
    }
    found.sort();
    found
}

/// What a run gave: stdout, stderr, exit status, and the files it made (`path\tsize\tdigest`).
type Run = (Vec<u8>, Vec<u8>, i32, String);

/// The run of `losat makeprofiledb words` in a fresh copy of `in/`, in the directory `dir`.
fn run(dir: &Path, words: &[String]) -> Run {
    let _ = fs::remove_dir_all(dir);
    copy_tree(&fixtures().join("in"), &dir.join("in"));
    let output = command(Path::new(env!("CARGO_BIN_EXE_losat")), dir)
        .arg("makeprofiledb")
        .args(words)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let mut files = String::new();
    for path in files_below(dir) {
        let bytes = fs::read(dir.join(&path)).unwrap();
        match masked(&path, &bytes) {
            Some(masked) => files.push_str(&format!("{path}\t*\t{:016x}\n", fnv1a64(&masked))),
            None => files.push_str(&format!(
                "{path}\t{}\t{:016x}\n",
                bytes.len(),
                fnv1a64(&bytes)
            )),
        }
    }
    (
        output.stdout,
        output.stderr,
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
        .join(format!("makeprofiledb-{name}-{}", std::process::id()))
}

#[test]
fn makeprofiledb_gives_ncbis_bytes() {
    let mut differ = Vec::new();
    for (name, words) in cases() {
        let dir = scratch(&name);
        let got = run(&dir, &words);
        let want = expected(&name);
        if got != want {
            let text = |bytes: &[u8]| {
                String::from_utf8_lossy(bytes)
                    .chars()
                    .take(800)
                    .collect::<String>()
            };
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

// `in/l_QR*.lst` (row PS10d): queries in a nucleotide coding with fewer bytes than NCBI's
// writer reads for their length (NCBI's two runs of QR01 and QR02 wrote different `.psq` bytes,
// QR04 ended with SIGSEGV; README.md). LOSAT refuses each, whatever the profile: QR05 has
// frequency ratios only, where NCBI also stops with its error of the query's length.
#[test]
fn makeprofiledb_refuses_nucleotide_data_shorter_than_its_length() {
    for (list, what) in [
        (
            "QR01",
            "1 byte of ncbi2na data for a length of 6 (NCBI's makeprofiledb reads 2 bytes",
        ),
        (
            "QR02",
            "2 bytes of ncbi4na data for a length of 6 (NCBI's makeprofiledb reads 3 bytes",
        ),
        (
            "QR03",
            "2 letters of iupacna data for a length of 6 (NCBI's makeprofiledb reads 3 bytes of \
             their ncbi4na form",
        ),
        (
            "QR04",
            "0 bytes of ncbi2na data for a length of 6 (NCBI's makeprofiledb reads 2 bytes",
        ),
        (
            "QR05",
            "0 letters of iupacna data for a length of 6 (NCBI's makeprofiledb reads 3 bytes of \
             their ncbi4na form",
        ),
    ] {
        let dir = scratch(list);
        let words = words(&format!("-in in/l_{list}.lst -out db"));
        let (out, err, status, _) = run(&dir, &words);
        assert_eq!(
            (
                String::from_utf8_lossy(&out),
                String::from_utf8_lossy(&err),
                status
            ),
            (
                "".into(),
                format!(
                    "Error: a profile whose query has {what} there, past the data) is not \
                     supported by LOSAT's MAKEPROFILEDB\n"
                )
                .into(),
                1
            ),
            "{list}"
        );
        let _ = fs::remove_dir_all(&dir);
    }
}

/// `text` with the date of each `Posted date:` line (outfmt 0's database footer: the build date
/// that makeprofiledb wrote into the `.pin` file) replaced by `<DATE>`.
fn posted_date_masked(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        match line.find("Posted date:") {
            Some(at) => {
                out.push_str(&line[..at]);
                out.push_str("Posted date:  <DATE>");
                if line.ends_with('\n') {
                    out.push('\n');
                }
            }
            None => out.push_str(line),
        }
    }
    out
}

// `rps/<row>.<fmt>`: NCBI rpsblast 2.17.0 `-query q3.faa -db db/db -outfmt <fmt>` on the database
// NCBI makeprofiledb built for the row (survey rows R01 and R05, 2026-11-21; QC01, a query of
// lower-case IUPACaa letters, 2026-11-23), its output in LOSAT's words; the database's build date
// (outfmt 0's `Posted date:`) is the day of the run.
#[test]
fn rpsblast_on_losat_databases_gives_ncbis_bytes() {
    for (row, case) in [("R01", "O00"), ("R05", "S_f2"), ("QC01", "QC01")] {
        let (_, words) = cases().into_iter().find(|(name, _)| name == case).unwrap();
        let dir = scratch(&format!("rps-{row}"));
        let (_, err, status, _) = run(&dir, &words);
        assert_eq!(status, 0, "{}", String::from_utf8_lossy(&err));
        fs::copy(fixtures().join("rps/q3.faa"), dir.join("q3.faa")).unwrap();
        for fmt in ["6", "0"] {
            let output = command(Path::new(env!("CARGO_BIN_EXE_losat")), &dir)
                .args(["rpsblast", "-query", "q3.faa", "-db", "db", "-outfmt", fmt])
                .stdin(Stdio::null())
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(0), "{row} outfmt {fmt}");
            let ncbi = fs::read(fixtures().join(format!("rps/{row}.{fmt}"))).unwrap();
            assert_eq!(
                posted_date_masked(&String::from_utf8_lossy(&output.stdout)),
                posted_date_masked(&String::from_utf8_lossy(
                    &losat::wording::losat_output_bytes(&ncbi)
                )),
                "{row} outfmt {fmt}"
            );
        }
        let _ = fs::remove_dir_all(&dir);
    }
}
