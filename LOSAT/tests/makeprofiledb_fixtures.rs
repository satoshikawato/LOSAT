//! `losat makeprofiledb` against NCBI BLAST+ 2.17.0 makeprofiledb (`-dbtype rps` and
//! `-dbtype delta`): each row of `tests/fixtures/makeprofiledb/cases.tsv` runs in its own copy of
//! `in/` and gives NCBI's standard output, standard error (in LOSAT's words, `losat::wording`)
//! and exit status (`ncbi/<row>.{out,err,rc}`) and, for each file the run wrote, its size and
//! FNV-1a digest (`ncbi/<row>.files`; the date of a `.pin` or `.pal` file replaced by `<DATE>`,
//! its size `*`; on Windows the files read back to Linux's, `as_on_linux`; where LMDB's pages are
//! not Linux's 4096 bytes, as on macOS on Apple silicon, an LMDB file is checked for whole pages
//! of the host's size, `losat_lmdb_line` and `ncbi_lmdb_line`), frozen from the oracle of row
//! PS10a (2026-11-21; README.md of the directory) and subsequent PS10/PS14 rows. NCBI is not run
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

/// The page size of the LMDB files of NCBI's frozen runs (Linux x86-64).
const LINUX_PAGE: usize = 4096;

/// Whether `name` is an LMDB file of a database (`.pdb`, `.ptf`, also of a volume).
fn is_lmdb(name: &str) -> bool {
    name.ends_with(".pdb") || name.ends_with(".ptf")
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:321,395 (GET_PAGESIZE)
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:4400-4402 (the page size of a new file)
/// The page size of the LMDB files that NCBI's makeprofiledb and LOSAT's write on this host
/// (`writedb::platform::host_page_size`): the operating system's page size, `sysconf` of
/// `_SC_PAGE_SIZE` off Windows (16384 bytes on macOS on Apple silicon), up to 32768; on Windows
/// 4096 bytes (`as_on_linux`).
fn host_lmdb_page() -> usize {
    #[cfg(unix)]
    let os = rustix::param::page_size();
    #[cfg(not(unix))]
    let os = LINUX_PAGE;
    os.min(32768)
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1036-1054 (MDB_meta, mm_psize)
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3747 (the meta page records it)
/// The page size that the LMDB file `bytes` records in its first meta page: after the page's
/// header of 16 bytes, the magic number and the version (4 bytes each), the address and the map
/// size (8 each), the `md_pad` of the first `MDB_db`.
fn recorded_page_size(bytes: &[u8]) -> Option<usize> {
    let field = bytes.get(40..44)?;
    Some(u32::from_le_bytes(field.try_into().unwrap()) as usize)
}

/// The line of the LMDB file `name` where LMDB's pages are not Linux's (`host_lmdb_page`): such a
/// file cannot be read back to NCBI's Linux file, whose B-trees fill pages of 4096 bytes, so the
/// line says whether the file is made of whole pages of `page` bytes, in place of its length and
/// digest.
fn pages_line(name: &str, whole: bool, page: usize) -> String {
    let made_of = if whole { "whole" } else { "not whole" };
    format!("{name}\tLMDB file of {made_of} pages of {page} bytes\n")
}

/// The line of LOSAT's LMDB file `name` (`bytes`) on a host of other pages than Linux's: whole
/// pages of the size the file records (`recorded_page_size`), which is the host's in NCBI's
/// (`ncbi_lmdb_line`). `None` for another file, on a host of Linux's pages, and for a file
/// shorter than two pages (no environment: the empty file of a refused run, whose line stays
/// exact).
fn losat_lmdb_line(name: &str, bytes: &[u8]) -> Option<String> {
    let host = host_lmdb_page();
    if host == LINUX_PAGE || !is_lmdb(name) || bytes.len() < 2 * host {
        return None;
    }
    let page = recorded_page_size(bytes)?;
    Some(pages_line(name, bytes.len().is_multiple_of(page), page))
}

/// The line of NCBI's LMDB file `name` of `length` bytes (frozen on Linux) on a host of other
/// pages: NCBI's makeprofiledb writes it here in whole pages of the host's size (LMDB's page
/// size, `host_lmdb_page`). `None` as for `losat_lmdb_line`.
fn ncbi_lmdb_line(name: &str, length: &str) -> Option<String> {
    let host = host_lmdb_page();
    let length: usize = length.parse().ok()?;
    if host == LINUX_PAGE || !is_lmdb(name) || length < 2 * LINUX_PAGE {
        return None;
    }
    Some(pages_line(name, true, host))
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
        if let Some(line) = losat_lmdb_line(&path, &bytes) {
            files.push_str(&line);
            continue;
        }
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
        fs::read_to_string(file("files"))
            .unwrap()
            .lines()
            .map(|line| {
                let mut fields = line.split('\t');
                let (path, length) = (fields.next().unwrap(), fields.next().unwrap_or(""));
                ncbi_lmdb_line(path, length).unwrap_or_else(|| format!("{line}\n"))
            })
            .collect(),
    )
}

fn scratch(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("makeprofiledb-{name}-{}", std::process::id()))
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3747 (the meta page records it)
/// LOSAT's LMDB files record the host's page size (LMDB's `GET_PAGESIZE`), the field that
/// `losat_lmdb_line` reads where the pages are not Linux's.
#[test]
fn makeprofiledb_s_lmdb_files_record_the_host_s_page_size() {
    let (_, words) = cases().into_iter().find(|(name, _)| name == "O00").unwrap();
    let dir = scratch("pages");
    let (_, _, status, _) = run(&dir, &words);
    assert_eq!(status, 0);
    for file in ["db.pdb", "db.ptf"] {
        let bytes = fs::read(dir.join(file)).unwrap();
        assert_eq!(recorded_page_size(&bytes), Some(host_lmdb_page()), "{file}");
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn makeprofiledb_gives_ncbis_bytes() {
    let mut differ = Vec::new();
    for (name, words) in cases() {
        // D36 uses a 248-byte basename so `.freq` (253 B) fits Unix `NAME_MAX` (255) while
        // `.wcounts` (256 B) fails in `x_WrapUpDelta`; on Windows `MAX_PATH` (260) prevents even
        // `.rps` from opening, so `makeprofiledb_delta_defers_wcounts_creation_until_wrapup`
        // covers late `.wcounts` failure on every platform.
        if cfg!(windows) && name == "D36" {
            continue;
        }
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

// `in/l_QR08.lst`-`in/l_QR15.lst` (row PS10b): PSSM and binary defline inputs where NCBI's run
// has no defined result or depends on process state. QR08-QR10: `iupacna` byte >= 0x80 in a pair
// position before any final odd overhang byte, where NCBI's `x_ConvertIupacnaTo4na` indexes its
// table with a signed `char` and reads before the table (odd overhang `B18` and bytes < 0x80 `B19`
// are defined, and `mol na` `B20` stops with NCBI's molecule error first). QR11 and QR14:
// `:Seq-interval` pointer of another type than `PssmParameters` (binary and ASN.1 text; via
// `Names::PROFILE`). QR12: `NULL` `Seq-id` inside `ASN1_BlastDefLine` in version 5 with 0 prior
// OIDs (default, `-taxid 9606`, and empty `-taxid_map`, where NCBI's `~CWriteDB_TaxID` aborts
// after posting 802.5; v4, non-empty `-taxid_map`, `org`, `-index false` on the first slot, and
// v5 after a prior profile are frozen in `RB17`-`RB29`). QR13: 1 GiB string length past EOF
// (host-dependent memory). QR15: `NULL` `Seq-id` profile followed by a missing input file in v5
// and v4, where `~CWriteDB_Impl`'s `Close` throws `WriteDbError::Exception` in the unwinding of
// the input error and NCBI aborts (rc 134); the refusal stays on one line without the trailing
// line feed of `what()`.
#[test]
fn makeprofiledb_refuses_binary_and_iupacna_high_byte_defects() {
    let iupacna_err = |byte: &str, pos: usize| {
        format!(
            "Error: a profile whose query has an iupacna byte {byte} at position {pos} (NCBI's \
             makeprofiledb indexes its iupacna-to-ncbi4na table with a signed char and reads \
             before the table) is not supported by LOSAT's MAKEPROFILEDB\n"
        )
    };
    let other_si_err = |file: &str, at: &str| {
        format!(
            "Error: a profile that LOSAT's PSSM reader does not read ({file}: \
             ':Seq-interval', a value of the type Seq-interval where a PssmParameters belongs \
             ({at}): whether makeprofiledb reads it depends on the types its process has met \
             before) is not supported by LOSAT's MAKEPROFILEDB\n"
        )
    };
    let null_post_802_5 = "Error: (802.5) CObjectOStream: error at byte 15: \
        Blast-def-line-set.E.seqid.E: cannot write unassigned member NULL element while writing \
        container \n";
    let qr12_v5_abort = format!(
        "{null_post_802_5}Error: an error (ERROR: byte 15: cannot write unassigned member NULL \
         element while writing container ) after which the database's writer stops with another \
         (Error: No tax info for any oid; NCBI's makeprofiledb aborts there) is not supported by \
         LOSAT's MAKEPROFILEDB\n"
    );
    let qr15_close_abort = format!(
        "{null_post_802_5}Error: an error (INPUT ERROR: in/nope.smp does not exists) after which \
         the database's writer stops with another (Error: byte 15: cannot write unassigned member \
         NULL element while writing container ; NCBI's makeprofiledb aborts there) is not \
         supported by LOSAT's MAKEPROFILEDB\n"
    );
    let cases = [
        ("QR08", "QR08", "-binary", iupacna_err("0x80", 0)),
        ("QR09", "QR09", "-binary", iupacna_err("0x80", 3)),
        ("QR10", "QR10", "-binary", iupacna_err("0x80", 0)),
        (
            "QR11",
            "QR11",
            "-binary",
            other_si_err("in/b_other_si.smp.bin", "byte 1250"),
        ),
        ("QR12_v5", "QR12", "", qr12_v5_abort.clone()),
        ("QR12_tid", "QR12", "-taxid 9606", qr12_v5_abort.clone()),
        (
            "QR12_tmap_empty",
            "QR12",
            "-taxid_map in/tmap_empty.txt",
            qr12_v5_abort,
        ),
        (
            "QR13",
            "QR13",
            "-binary",
            "Error: a profile that LOSAT's PSSM reader does not read (in/b_str_1g.smp.bin: a \
             string whose length (1073741824 bytes) passes the end of the file: makeprofiledb \
             allocates 1073741825 bytes for it, which succeeds or fails by the host's memory) is \
             not supported by LOSAT's MAKEPROFILEDB\n"
                .to_string(),
        ),
        (
            "QR14",
            "QR14",
            "",
            other_si_err("in/t_other_si.smp", "line 201"),
        ),
        ("QR15_v5", "QR15", "", qr15_close_abort.clone()),
        ("QR15_v4", "QR15", "-blastdb_version 4", qr15_close_abort),
    ];
    for (tag, list, extra, expected_err) in cases {
        let dir = scratch(tag);
        let argstr = if extra.is_empty() {
            format!("-in in/l_{list}.lst -out db")
        } else {
            format!("{extra} -in in/l_{list}.lst -out db")
        };
        let words = words(&argstr);
        let (out, got, status, _) = run(&dir, &words);
        assert_eq!(
            (
                String::from_utf8_lossy(&out),
                String::from_utf8_lossy(&got),
                status
            ),
            ("".into(), expected_err.into(), 1),
            "{tag}"
        );
        let _ = fs::remove_dir_all(&dir);
    }
}

#[cfg(unix)]
#[test]
fn makeprofiledb_binary_handles_non_utf8_paths() {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;

    let dir = scratch("bin-nonutf8");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("in")).unwrap();
    let smp_bytes = fs::read(fixtures().join("in/nv01.smp.bin")).unwrap();
    let bad_name = OsStr::from_bytes(b"bad_\xff_name.smp.bin");
    fs::write(dir.join("in").join(bad_name), &smp_bytes).unwrap();
    fs::write(dir.join("in/list.lst"), b"in/bad_\xff_name.smp.bin\n").unwrap();
    let out_name = OsStr::from_bytes(b"db_\xfe_out");
    let output = command(Path::new(env!("CARGO_BIN_EXE_losat")), &dir)
        .args([
            OsStr::new("makeprofiledb"),
            OsStr::new("-binary"),
            OsStr::new("-in"),
            OsStr::new("in/list.lst"),
            OsStr::new("-out"),
            out_name,
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    for ext in [
        b".rps".as_slice(),
        b".loo",
        b".aux",
        b".freq",
        b".pin",
        b".phr",
        b".psq",
    ] {
        let mut file = b"db_\xfe_out".to_vec();
        file.extend_from_slice(ext);
        assert!(
            dir.join(OsStr::from_bytes(&file)).is_file(),
            "missing output file {:?}",
            String::from_utf8_lossy(&file)
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

// `in/l_QR16.lst`-`in/l_QR22.lst` (row PS14): `-dbtype cobalt` and `-dbtype delta` inputs where
// NCBI BLAST+ 2.17.0 has undefined behaviour: `byRow TRUE` (`QR16`) and `numRows < 28` with
// weighted residue frequencies (`QR17`) copy into an empty `std::list::begin()` (overwriting the
// list header node), empty `numIndeptObsr` (`QR18`) dereferences `end()`, a partial trailing
// column (`QR19`, length not a multiple of 28) increments past `freqs.end()`, and out-of-range or
// NaN floats (`QR20`-`QR22`) convert to `Uint4` with undefined behaviour.
#[test]
fn makeprofiledb_refuses_delta_ub_and_cobalt() {
    let cases = [
        (
            "cobalt",
            "-dbtype cobalt -in in/l_D01.lst -out db",
            "Error: -dbtype cobalt (only -dbtype rps and -dbtype delta are built) is not \
             supported by LOSAT's MAKEPROFILEDB\n",
        ),
        (
            "QR16",
            "-dbtype delta -in in/l_QR16.lst -out db",
            "Error: a profile stored by row (byRow TRUE) for -dbtype delta (NCBI's makeprofiledb \
             copies the transposed weighted residue frequencies into an empty list's begin \
             iterator) is not supported by LOSAT's MAKEPROFILEDB\n",
        ),
        (
            "QR17",
            "-dbtype delta -in in/l_QR17.lst -out db",
            "Error: a profile with numRows 26 (fewer than 28) and weighted residue frequencies \
             for -dbtype delta (NCBI's makeprofiledb copies them into an empty list's begin \
             iterator) is not supported by LOSAT's MAKEPROFILEDB\n",
        ),
        (
            "QR18",
            "-dbtype delta -in in/l_QR18.lst -out db",
            "Error: a profile with an empty numIndeptObsr list for -dbtype delta (NCBI's \
             makeprofiledb dereferences the end iterator of the empty list) is not supported by \
             LOSAT's MAKEPROFILEDB\n",
        ),
        (
            "QR19",
            "-dbtype delta -in in/l_QR19.lst -out db",
            "Error: a profile whose weightedResFreqsPerPos has 29 entries (not a multiple of 28) \
             and every full column sums within 0.0001 of 1.0 for -dbtype delta (NCBI's \
             makeprofiledb increments past the end of the list) is not supported by LOSAT's \
             MAKEPROFILEDB\n",
        ),
        (
            "QR20",
            "-dbtype delta -obsr_threshold -10 -exclude_invalid false -in in/l_QR20.lst -out db",
            "Error: a numIndeptObsr value of -2 for -dbtype delta (NCBI's conversion of (-2 + 1) * \
             1000 to a 32-bit unsigned integer has no defined result) is not supported by LOSAT's \
             MAKEPROFILEDB\n",
        ),
        (
            "QR21",
            "-dbtype delta -exclude_invalid false -in in/l_QR21.lst -out db",
            "Error: a weightedResFreqsPerPos value of -0.5 for -dbtype delta (NCBI's conversion of \
             -0.5 * 1000 to a 32-bit unsigned integer has no defined result) is not supported by \
             LOSAT's MAKEPROFILEDB\n",
        ),
        (
            "QR22",
            "-dbtype delta -in in/l_QR22.lst -out db",
            "Error: a numIndeptObsr value of NaN for -dbtype delta (NCBI's conversion of (NaN + 1) \
             * 1000 to a 32-bit unsigned integer has no defined result) is not supported by LOSAT's \
             MAKEPROFILEDB\n",
        ),
    ];
    for (tag, argstr, expected_err) in cases {
        let dir = scratch(tag);
        let words = words(argstr);
        let (out, got, status, _) = run(&dir, &words);
        assert_eq!(
            (
                String::from_utf8_lossy(&out),
                String::from_utf8_lossy(&got),
                status
            ),
            ("".into(), expected_err.into(), 1),
            "{tag}"
        );
        let _ = fs::remove_dir_all(&dir);
    }
}

#[cfg(unix)]
#[test]
fn makeprofiledb_delta_handles_non_utf8_paths() {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;

    let dir = scratch("delta-nonutf8");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("in")).unwrap();
    let smp_bytes = fs::read(fixtures().join("in/c1.smp")).unwrap();
    let bad_name = OsStr::from_bytes(b"delta_\xff_name.smp");
    fs::write(dir.join("in").join(bad_name), &smp_bytes).unwrap();
    fs::write(dir.join("in/list.lst"), b"in/delta_\xff_name.smp\n").unwrap();
    let out_name = OsStr::from_bytes(b"db_\xfe_delta");
    let output = command(Path::new(env!("CARGO_BIN_EXE_losat")), &dir)
        .args([
            OsStr::new("makeprofiledb"),
            OsStr::new("-dbtype"),
            OsStr::new("delta"),
            OsStr::new("-in"),
            OsStr::new("in/list.lst"),
            OsStr::new("-out"),
            out_name,
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    for ext in [
        b".wcounts".as_slice(),
        b".obsr",
        b".rps",
        b".loo",
        b".aux",
        b".freq",
        b".pin",
        b".phr",
        b".psq",
    ] {
        let mut file = b"db_\xfe_delta".to_vec();
        file.extend_from_slice(ext);
        assert!(
            dir.join(OsStr::from_bytes(&file)).is_file(),
            "missing output file {:?}",
            String::from_utf8_lossy(&file)
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

/// Copies `fixtures()/in` into a fresh `dir` and runs `losat makeprofiledb` with `args` and
/// optional `TMPDIR` override.
fn run_makeprofiledb_in(dir: &Path, tmpdir: Option<&Path>, args: &[&str]) -> std::process::Output {
    let mut cmd = command(Path::new(env!("CARGO_BIN_EXE_losat")), dir);
    if let Some(tmp) = tmpdir {
        cmd.env("TMPDIR", tmp);
    }
    cmd.arg("makeprofiledb")
        .args(args)
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn reset_fixture_dir(tag: &str) -> PathBuf {
    let dir = scratch(tag);
    let _ = fs::remove_dir_all(&dir);
    copy_tree(&fixtures().join("in"), &dir.join("in"));
    dir
}

#[test]
fn makeprofiledb_delta_defers_wcounts_creation_until_wrapup() {
    // Case (a): an input profile named `db.00.wcounts` with `-max_smp_vol 2` and a pre-existing
    // `db.wcounts` (which `s_DeleteMakeprofileDb` removes, skipping numbered `.00.wcounts`).
    // NCBI reads `db.00.wcounts` in passes 2 and 3 and truncates it only in `x_WrapUpDelta`.
    let dir = reset_fixture_dir("delta-wrapup-input-name");
    fs::write(dir.join("db.wcounts"), b"dummy").unwrap();
    fs::copy(dir.join("in/c1.smp"), dir.join("db.00.wcounts")).unwrap();
    fs::write(
        dir.join("in/l_wrapup.lst"),
        b"db.00.wcounts\nin/c2.smp\nin/c3.smp\n",
    )
    .unwrap();
    let wrapup_out = run_makeprofiledb_in(
        &dir,
        None,
        &[
            "-dbtype",
            "delta",
            "-max_smp_vol",
            "2",
            "-title",
            "in/l_D01.lst",
            "-in",
            "in/l_wrapup.lst",
            "-out",
            "db",
        ],
    );
    assert_eq!(
        wrapup_out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&wrapup_out.stderr)
    );
    let want_d05 = expected("D05");
    assert_eq!(wrapup_out.stdout, want_d05.0);
    assert_eq!(wrapup_out.stderr, want_d05.1);
    let mut got_files = String::new();
    for path in files_below(&dir) {
        let bytes = as_on_linux(&path, fs::read(dir.join(&path)).unwrap());
        if let Some(line) = losat_lmdb_line(&path, &bytes) {
            got_files.push_str(&line);
            continue;
        }
        match masked(&path, &bytes) {
            Some(m) => got_files.push_str(&format!("{path}\t*\t{:016x}\n", fnv1a64(&m))),
            None => got_files.push_str(&format!(
                "{path}\t{}\t{:016x}\n",
                bytes.len(),
                fnv1a64(&bytes)
            )),
        }
    }
    assert_eq!(got_files, want_d05.3);
    let _ = fs::remove_dir_all(&dir);

    // Case (b) cross-platform: a non-empty directory at `db.wcounts` blocks `File::create` in
    // `x_WrapUpDelta` after Pass 2 finishes, leaving the v4 protein DB on disk with exit 2.
    let dir = reset_fixture_dir("delta-wrapup-open-fail");
    fs::create_dir_all(dir.join("db.wcounts")).unwrap();
    fs::write(dir.join("db.wcounts/keep"), b"x").unwrap();
    let fail_out = run_makeprofiledb_in(
        &dir,
        None,
        &[
            "-dbtype",
            "delta",
            "-blastdb_version",
            "4",
            "-in",
            "in/l_D36.lst",
            "-out",
            "db",
        ],
    );
    assert_eq!(fail_out.status.code(), Some(2));
    assert!(fail_out.stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&fail_out.stderr),
        "Error: ERROR: Failed to open output .wcounts file\n"
    );
    for ext in [".pin", ".phr", ".psq", ".psd", ".psi", ".pog"] {
        assert!(dir.join(format!("db{ext}")).is_file(), "missing db{ext}");
    }
    for ext in [".rps", ".loo", ".aux", ".freq", ".obsr"] {
        assert!(!dir.join(format!("db{ext}")).exists(), "unexpected db{ext}");
    }
    let _ = fs::remove_dir_all(&dir);
}

#[cfg(unix)]
#[test]
fn makeprofiledb_delta_refuses_missing_tmpdir_after_pass2() {
    let dir = reset_fixture_dir("delta-missing-tmpdir");
    let missing_tmp = dir.join("no_such_tmpdir");

    // When Pass 2 succeeds, missing TMPDIR is refused (NCBI's mkstemp error contains the PID).
    let refused_out = run_makeprofiledb_in(
        &dir,
        Some(&missing_tmp),
        &["-dbtype", "delta", "-in", "in/l_D36.lst", "-out", "db"],
    );
    assert_eq!(refused_out.status.code(), Some(1));
    assert!(refused_out.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&refused_out.stderr);
    assert!(
        stderr.contains("temporary directory")
            && stderr.contains("not supported by LOSAT's MAKEPROFILEDB"),
        "unexpected stderr: {stderr}"
    );

    // When Pass 2 fails (D26: scaling factors do not match), NCBI never enters x_UpdateDelta, so
    // the Pass-2 error takes precedence over missing TMPDIR.
    let pass2_out = run_makeprofiledb_in(
        &dir,
        Some(&missing_tmp),
        &["-dbtype", "delta", "-in", "in/l_D26.lst", "-out", "db2"],
    );
    assert_eq!(pass2_out.status.code(), Some(1));
    assert!(pass2_out.stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&pass2_out.stderr),
        "Error: ERROR: Scaling factors do not match\n"
    );
    let _ = fs::remove_dir_all(&dir);
}

// `rps/<row>.<fmt>`: NCBI rpsblast 2.17.0 `-query <query> -db db/db -outfmt <fmt>` on the database
// NCBI makeprofiledb built for the row (survey rows R01 and R05, 2026-11-21; QC01, a query of
// lower-case IUPACaa letters, 2026-11-23; B01 and B03, ASN.1 binary score and frequency-ratio
// profiles, 2026-11-25; D01 and D04, text and binary `-dbtype delta` profiles, 2026-11-26), its
// output in LOSAT's words; the database's build date (outfmt 0's `Posted date:`) is the day of
// the run.
#[test]
fn rpsblast_on_losat_databases_gives_ncbis_bytes() {
    for (row, case, query) in [
        ("R01", "O00", "q3.faa"),
        ("R05", "S_f2", "q3.faa"),
        ("QC01", "QC01", "q3.faa"),
        ("B01", "B01", "q3.faa"),
        ("B03", "B03", "q3.faa"),
        ("D01", "D01", "qcdd.faa"),
        ("D04", "D04", "qcdd.faa"),
    ] {
        let (_, words) = cases().into_iter().find(|(name, _)| name == case).unwrap();
        let dir = scratch(&format!("rps-{row}"));
        let (_, err, status, _) = run(&dir, &words);
        assert_eq!(status, 0, "{}", String::from_utf8_lossy(&err));
        fs::copy(fixtures().join(format!("rps/{query}")), dir.join(query)).unwrap();
        for fmt in ["6", "0"] {
            let output = command(Path::new(env!("CARGO_BIN_EXE_losat")), &dir)
                .args(["rpsblast", "-query", query, "-db", "db", "-outfmt", fmt])
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
