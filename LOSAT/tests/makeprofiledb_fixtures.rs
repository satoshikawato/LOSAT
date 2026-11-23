//! `losat makeprofiledb` against NCBI BLAST+ 2.17.0 makeprofiledb (`-dbtype rps`): each row of
//! `tests/fixtures/makeprofiledb/cases.tsv` runs in its own copy of `in/` and gives NCBI's standard
//! output, standard error (in LOSAT's words, `losat::wording`) and exit status
//! (`ncbi/<row>.{out,err,rc}`) and, for each file the run wrote, its size and FNV-1a digest
//! (`ncbi/<row>.files`; the date of a `.pin` or `.pal` file replaced by `<DATE>`, its size `*`;
//! on Windows the files read back to Linux's, `as_on_linux`),
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

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3979-4003 (mdb_env_map)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:438-448
/// The bytes `bytes` of the file `name` that a run wrote, as LOSAT writes them on Linux, the
/// platform of NCBI's frozen files. On Windows LOSAT writes what NCBI's code writes there
/// (`writedb::platform`; not checked against a makeprofiledb.exe): an LMDB file (`.pdb`, `.ptf`)
/// is as long as its map (500000 bytes, grown where the data needs more), which its meta pages
/// record, and the LMDB name of a version 5 index file (`.pin`) is the database's name after its
/// last `\`. The bytes are read back to Linux's after the check that the file is as long as the
/// map of its newer meta page: the file ends after that meta page's last page, and the map size
/// of both meta pages is Linux's (700000000 bytes for `.pdb`, 300000000 for `.ptf`); the LMDB name
/// is the database's name after its last `/`.
fn as_on_linux(name: &str, bytes: Vec<u8>) -> Vec<u8> {
    if !cfg!(windows) {
        return bytes;
    }
    if name.ends_with(".pin") {
        let int4 = |at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
        // Version 5: the version, the molecule type, the volume, the title, then the LMDB name,
        // each string after its length. A refused run leaves an empty file.
        if bytes.len() < 16 || int4(0) != 5 {
            return bytes;
        }
        let lmdb_at = 12 + 4 + int4(12);
        let lmdb_end = lmdb_at + 4 + int4(lmdb_at);
        let lmdb = &bytes[lmdb_at + 4..lmdb_end];
        let base = lmdb.rsplit(|&byte| byte == b'/').next().unwrap();
        let mut out = bytes[..lmdb_at].to_vec();
        out.extend_from_slice(&(base.len() as u32).to_be_bytes());
        out.extend_from_slice(base);
        out.extend_from_slice(&bytes[lmdb_end..]);
        return out;
    }
    let linux_map_size: u64 = if name.ends_with(".pdb") {
        700_000_000
    } else if name.ends_with(".ptf") {
        300_000_000
    } else {
        return bytes;
    };
    // An `MDB_meta` after the page's header of 16 bytes: the magic number and the version (4
    // bytes each), the address (8), the map size (8), two `MDB_db` (48 bytes each), the last page
    // (8), the transaction id (8); LMDB's pages on Windows are of 4096 bytes.
    const PAGE: usize = 4096;
    const MAP_SIZE: usize = 16 + 4 + 4 + 8;
    const LAST_PAGE: usize = MAP_SIZE + 8 + 2 * 48;
    const TXNID: usize = LAST_PAGE + 8;
    let field = |page: usize, at: usize| {
        u64::from_le_bytes(bytes[page * PAGE + at..][..8].try_into().unwrap())
    };
    // A file without its two meta pages (a run refused before the environment was made) stays
    // as it is.
    if bytes.len() < 2 * PAGE {
        return bytes;
    }
    let newer = usize::from(field(1, TXNID) > field(0, TXNID));
    assert_eq!(
        bytes.len() as u64,
        field(newer, MAP_SIZE),
        "{name}: an LMDB file on Windows is as long as its map"
    );
    let end = (field(newer, LAST_PAGE) as usize + 1) * PAGE;
    let mut out = bytes[..end].to_vec();
    for page in 0..2 {
        out[page * PAGE + MAP_SIZE..][..8].copy_from_slice(&linux_map_size.to_le_bytes());
    }
    out
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
        let bytes = as_on_linux(&path, fs::read(dir.join(&path)).unwrap());
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

// `in/l_QR06.lst`, `in/l_QR07.lst` (row PS10d): profiles with frequency ratios only whose query
// NCBI's run does not end in a defined way. QR06 (ncbi8na alone, review O6): NCBI's conversion
// error, then its writer's error posted as it closes, then the taxonomy writer's exception out
// of its destructor (`terminate`, exit status 134). QR07 (a gap): NCBI's conversion reads past
// its table of codings.
#[test]
fn makeprofiledb_refuses_a_query_where_ncbi_has_no_result() {
    for (list, err) in [
        (
            "QR06",
            "Error: LOSAT database creation error: Unable to process sequence for entry \
             [CDD:237976].\nError: an error (ERROR: Unknown conversion.) after which the \
             database's writer stops with another (Error: No tax info for any oid; NCBI's \
             makeprofiledb aborts there) is not supported by LOSAT's MAKEPROFILEDB\n",
        ),
        (
            "QR07",
            "Error: a profile of frequency ratios whose query's Seq-data is a gap (NCBI's \
             sequence conversion reads its coding past the end of its table) is not supported \
             by LOSAT's MAKEPROFILEDB\n",
        ),
    ] {
        let dir = scratch(list);
        let words = words(&format!("-in in/l_{list}.lst -out db"));
        let (out, got, status, _) = run(&dir, &words);
        assert_eq!(
            (
                String::from_utf8_lossy(&out),
                String::from_utf8_lossy(&got),
                status
            ),
            ("".into(), err.into(), 1),
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
