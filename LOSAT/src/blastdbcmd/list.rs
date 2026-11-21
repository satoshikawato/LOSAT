//! `losat blastdbcmd -list DIR`: the databases of a directory (`FindBlastDBs`) written with
//! `-list_outfmt` (`CBlastDbFormatter`).

use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::Ostream;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

use super::{refusal, Diag};

/// The separator NCBI puts after a directory (`CDirEntry::GetPathSeparator`).
const SEPARATOR: u8 = if cfg!(windows) { b'\\' } else { b'/' };

/// A database the walk found (`SSeqDBInitInfo`): its name in double quotes, as NCBI keeps
/// it, and its molecule type.
#[derive(Clone, Debug, PartialEq, Eq)]
struct FoundDb {
    quoted: Vec<u8>,
    seqtype: SeqType,
}

impl FoundDb {
    /// `NStr::Replace(m_BlastDbName, "\"", kEmptyStr)`: the name without the double quotes.
    fn file_name(&self) -> Vec<u8> {
        self.quoted.iter().copied().filter(|&b| b != b'"').collect()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1366-1427
// (`CBlastDbFinder`, cited at `find_files`): the names found are `std::string`s, opened and
// listed as their bytes.
/// The path of the bytes of a name found in a directory (`utils::os_bytes::path_from_bytes`): a
/// directory or database whose name is not UTF-8 is entered and listed as NCBI does (probe R6b of
/// line LIB-B).
fn path_of(bytes: &[u8]) -> std::path::PathBuf {
    crate::utils::os_bytes::path_from_bytes(bytes)
}

/// The name of the listed database `db` in its quotes, as the database reader takes it, which
/// keeps the names of database files as UTF-8 strings (`seqdb::alias::SearchPath::new`): a
/// name that is not UTF-8 is LOSAT's refusal of the letters that open the database (`%t %n %l
/// %d %U`; NCBI BLAST+ 2.17.0 opens it, probe R6a of line LIB-B).
fn reader_name(db: &FoundDb) -> Result<&str, ListStop> {
    std::str::from_utf8(&db.quoted).map_err(|_| {
        ListStop::Refused(refusal(&format!(
            "{}: the name of the database is not valid UTF-8 (LOSAT keeps the names of database files as UTF-8 strings)",
            String::from_utf8_lossy(&db.file_name())
        )))
    })
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:3019-3076
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1350-1364
// The masks are `*.pin`, `*.pal`, `*.nin`, `*.nal` (with case); a file is not opened, so an
// empty or stray `x.pin` is a database, and so is a link that leads nowhere (`IsDir` is
// false for it). The sub-directories (`*`) are entered with `-recursive`, links to
// directories too. The order of the directory does not matter: the list is sorted.
/// `FindFilesInDir(CDir(dir), masks, {"*"}, finder, fFF_File | recursive)`: the databases of
/// the files of `dir` whose names end with one of `extensions` (`.pin` ...).
fn find_files(dir: &[u8], extensions: &[&str], recursive: bool, found: &mut Vec<FoundDb>) {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3879
    // An empty directory is the current one, its names without a prefix (probes K36, K37).
    let read = if dir.is_empty() {
        std::fs::read_dir(".")
    } else {
        std::fs::read_dir(path_of(dir))
    };
    let Ok(entries) = read else {
        return;
    };
    let mut prefix = dir.to_vec();
    let ends_in_separator = |b: u8| b == SEPARATOR || (cfg!(windows) && b == b'/');
    if !prefix.is_empty() && !prefix.last().is_some_and(|&b| ends_in_separator(b)) {
        prefix.push(SEPARATOR);
    }
    for entry in entries.flatten() {
        let name = entry.file_name();
        let name = name.as_encoded_bytes();
        let mut full = prefix.clone();
        full.extend_from_slice(name);
        let is_dir = || std::fs::metadata(path_of(&full)).is_ok_and(|meta| meta.is_dir());
        let matches = extensions.iter().any(|ext| name.ends_with(ext.as_bytes()));
        let mut known_dir = None;
        if matches {
            let dir = is_dir();
            known_dir = Some(dir);
            if !dir {
                let base = &full[..full.len() - 4];
                let mut quoted = Vec::with_capacity(base.len() + 2);
                quoted.push(b'"');
                quoted.extend_from_slice(base);
                quoted.push(b'"');
                found.push(FoundDb {
                    quoted,
                    seqtype: if full[full.len() - 3] == b'n' {
                        SeqType::Nucleotide
                    } else {
                        SeqType::Protein
                    },
                });
            }
        }
        if recursive && known_dir.unwrap_or_else(is_dir) {
            find_files(&full, extensions, recursive, found);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1366-1427
// A database goes when its quoted name holds a name of a `DBLIST` line of an alias file of
// the list anywhere: `"db/p5"` of blastdb_aliastool takes `"db/p5"` and not `db/p5n` (the
// quotes anchor it), `n5` of a hand-written list takes `n5`, `n55` and `n5x` (probes L03,
// J162, J167 of line BDC). Only spaces part the names; a `\r` stays on the last.
/// `-remove_redundant_dbs`: drops from `dbs` every database an alias file of `dbs` names.
fn remove_alias_components(dbs: &mut Vec<FoundDb>) {
    let mut remove = BTreeSet::new();
    for db in dbs.iter() {
        let name = db.file_name();
        let (alias, index) = match db.seqtype {
            SeqType::Nucleotide => (".nal", ".nin"),
            SeqType::Protein => (".pal", ".pin"),
        };
        let mut path = name.clone();
        path.extend_from_slice(alias.as_bytes());
        if !std::fs::metadata(path_of(&path)).is_ok_and(|meta| meta.is_file()) {
            path = name;
            path.extend_from_slice(index.as_bytes());
        }
        if path.last() != Some(&b'l') {
            continue;
        }
        let Ok(text) = std::fs::read(path_of(&path)) else {
            continue;
        };
        for line in text.split(|&b| b == b'\n') {
            if !line.starts_with(b"DBLIST") {
                continue;
            }
            let tokens: Vec<&[u8]> = line
                .split(|&b| b == b' ')
                .filter(|token| !token.is_empty())
                .collect();
            remove.extend(tokens.iter().skip(1).map(|token| token.to_vec()));
        }
    }
    for token in &remove {
        dbs.retain(|db| !contains(&db.quoted, token));
    }
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.is_empty()
        || haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1428-1462
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdb.hpp:1550-1561
// The names compare as their bytes with the quotes: `"db/a b"` before `"db/a!"` before
// `"db/x"` before `"db/x.00"`; for one name the protein database first (`eProtein` is 0).
/// `FindBlastDBs(dir, dbtype, recursive, true, remove_redundant)`.
fn find_blast_dbs(
    dir: &str,
    dbtype: &str,
    recursive: bool,
    remove_redundant: bool,
) -> Vec<FoundDb> {
    let mut extensions = Vec::new();
    if dbtype != "nucl" {
        extensions.extend([".pin", ".pal"]);
    }
    if dbtype != "prot" {
        extensions.extend([".nin", ".nal"]);
    }
    let mut dbs = Vec::new();
    find_files(dir.as_bytes(), &extensions, recursive, &mut dbs);
    if remove_redundant {
        remove_alias_components(&mut dbs);
    }
    let rank = |seqtype: SeqType| match seqtype {
        SeqType::Protein => 0,
        SeqType::Nucleotide => 1,
    };
    dbs.sort_by(|a, b| {
        a.quoted
            .cmp(&b.quoted)
            .then(rank(a.seqtype).cmp(&rank(b.seqtype)))
    });
    dbs
}

/// `CBlastDbFormatter`: the format with `%%` made `%`, and the place and letter of each
/// replacement.
#[derive(Debug, PartialEq, Eq)]
struct ListFormat {
    spec: Vec<u8>,
    offsets: Vec<usize>,
    types: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_formatter.cpp:42-64
// After the erase the loop steps over the second `%` and the byte after it (`i++` and the
// loop's own); a `%` at the end reads the string's NUL as its letter.
impl ListFormat {
    /// The format `fmt`; `Err` for NCBI's `Invalid format specification`.
    fn parse(fmt: &[u8]) -> Result<Self, ()> {
        let mut spec = fmt.to_vec();
        let (mut offsets, mut types) = (Vec::new(), Vec::new());
        let mut i = 0;
        while i < spec.len() {
            let next = spec.get(i + 1).copied().unwrap_or(0);
            if spec[i] == b'%' && next == b'%' {
                spec.remove(i);
                i += 2;
                continue;
            }
            if spec[i] == b'%' {
                offsets.push(i);
                types.push(next);
            }
            i += 1;
        }
        if offsets.is_empty() {
            return Err(());
        }
        Ok(ListFormat {
            spec,
            offsets,
            types,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_formatter.cpp:175-195
    /// `x_Replacer`: the format with `data[i]` in place of replacement `i`.
    fn replace(&self, data: &[Vec<u8>]) -> Vec<u8> {
        let mut line = Vec::new();
        let mut fmt_idx = 0;
        for (offset, value) in self.offsets.iter().zip(data) {
            line.extend_from_slice(&self.spec[fmt_idx.min(*offset)..*offset]);
            line.extend_from_slice(value);
            fmt_idx = offset + 2;
        }
        if fmt_idx <= self.spec.len() {
            line.extend_from_slice(&self.spec[fmt_idx..]);
        }
        line
    }
}

/// What ends the list before its last database.
enum ListStop {
    /// NCBI's diagnostic line (`BLAST options error: ...`, `BLAST Database error: ...`) and
    /// exit status, after the lines already written.
    Ncbi(Vec<u8>, i32),
    /// LOSAT's refusal: nothing of the list is written.
    Refused(String),
}

impl From<SeqDbError> for ListStop {
    // NCBI's formatter opens the database for the fields that read it and throws to the
    // application's `CATCH_ALL` (blast_app_util.hpp:247-250, cited at `database_error`): an
    // alias file's value that does not convert is `Error: <what()>`, exit status 255.
    fn from(error: SeqDbError) -> Self {
        if let Some(exception) = error.ncbi_exception() {
            return ListStop::Ncbi([b"Error: ".as_slice(), &exception.what()].concat(), 255);
        }
        match error.ncbi_message() {
            Some(text) => ListStop::Ncbi([b"BLAST Database error: ".as_slice(), &text].concat(), 2),
            None => ListStop::Refused(crate::algorithm::common::database::losat_reader_error(
                &error,
                super::PROGRAM,
            )),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2620-2656
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2657-2666
/// `SeqDB_GetFileExtensions`: the extensions (without the molecule letter) of a volume's
/// files.
fn file_extensions(version5: bool) -> Vec<&'static str> {
    let mut extensions = vec!["al", "in", "hr", "sq", "ni", "nd"];
    if !version5 {
        extensions.extend(["si", "sd"]);
    }
    extensions.extend(["pi", "pd", "js"]);
    if version5 {
        extensions.extend(["db", "os", "ot", "tf", "to", "db-lock", "tf-lock"]);
    }
    extensions.extend(["aa", "ab", "ac", "og", "hi", "hd", "ti", "td", "xm"]);
    extensions
}

/// The length of the regular file `path` (`CFile::Exists`, `GetLength`; links followed).
fn file_length(path: &str) -> Option<u64> {
    std::fs::metadata(path)
        .ok()
        .filter(|meta| meta.is_file())
        .map(|meta| meta.len())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1464-1520
// The name list of a listed database is its name in double quotes, so the files of the
// whole database (`.nal`, `.ndb`, ...) are looked for as `<dir>/"<name>".<ext>` and not
// found: `%U` of a database of several volumes is the sum of its volumes' files (probe L02,
// J161 of line BDC). The volumes are those of the whole alias tree, each once.
/// `GetDiskUsage` of the listed database `db` (of version 5 when `version5`).
fn disk_usage(db: &FoundDb, version5: bool, path: &SearchPath) -> Result<u64, ListStop> {
    let name = reader_name(db)?;
    let (volumes, _) = SeqDb::volume_and_alias_paths(name, db.seqtype, path)?;
    let mol = match db.seqtype {
        SeqType::Protein => 'p',
        SeqType::Nucleotide => 'n',
    };
    let extensions = file_extensions(version5);
    let mut total = 0u64;
    let mut dirname = String::new();
    for volume in &volumes {
        for ext in &extensions {
            let file = format!("{volume}.{mol}{ext}");
            if let Some(length) = file_length(&file) {
                total += length;
                dirname = Path::new(&file)
                    .parent()
                    .map(|dir| dir.to_string_lossy().into_owned())
                    .unwrap_or_default();
            }
        }
    }
    if volumes.len() > 1 {
        if name.split(' ').filter(|piece| !piece.is_empty()).count() > 1 {
            return Err(ListStop::Ncbi(
                format!(
                    "BLAST Database error: Cannot compute disk usage for multiple BLASTDBs \
                     (i.e.: '{name}') at once. Please try again using one BLASTDB at a time."
                )
                .into_bytes(),
                2,
            ));
        }
        for ext in &extensions {
            let file = Path::new(&dirname).join(format!("{name}.{mol}{ext}"));
            if let Some(length) = file_length(&file.to_string_lossy()) {
                total += length;
            }
        }
    }
    Ok(total)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_formatter.cpp:68-162
// The database opens at its first letter that reads it (`CSeqDB("\"name\"", type)`, the
// name in its quotes: one name, looked for on the search path as written), so an unknown
// letter after one that opens a damaged database is that database's error.
/// `x_InitBlastDb`: the database `db`, opened at the first call.
fn open<'a>(
    opened: &'a mut Option<SeqDb>,
    db: &FoundDb,
    path: &SearchPath,
) -> Result<&'a SeqDb, ListStop> {
    if opened.is_none() {
        *opened = Some(SeqDb::open_in(reader_name(db)?, db.seqtype, path)?);
    }
    Ok(opened.as_ref().expect("opened above"))
}

/// `Write(db)`: the line of `db` (without its line feed).
fn write_line(format: &ListFormat, db: &FoundDb, path: &SearchPath) -> Result<Vec<u8>, ListStop> {
    let mut opened: Option<SeqDb> = None;
    let mut data = Vec::with_capacity(format.types.len());
    for &kind in &format.types {
        let value = match kind {
            b'f' => db.file_name(),
            b't' => open(&mut opened, db, path)?.title().to_vec(),
            b'n' => open(&mut opened, db, path)?
                .printed_num_seqs()
                .to_string()
                .into_bytes(),
            b'l' => open(&mut opened, db, path)?
                .printed_total_length()
                .to_string()
                .into_bytes(),
            b'p' => match db.seqtype {
                SeqType::Protein => b"Protein".to_vec(),
                SeqType::Nucleotide => b"Nucleotide".to_vec(),
            },
            b'd' => {
                let date = open(&mut opened, db, path)?.date().to_string();
                if date.contains('\u{fffd}') {
                    return Err(ListStop::Refused(refusal(&format!(
                        "{}: the date of the database is not UTF-8 text",
                        String::from_utf8_lossy(&db.file_name())
                    ))));
                }
                date.into_bytes()
            }
            b'U' => {
                let version5 = open(&mut opened, db, path)?.is_version5();
                disk_usage(db, version5, path)?.to_string().into_bytes()
            }
            b'v' => if open(&mut opened, db, path)?.is_version5() {
                b"5"
            } else {
                b"4"
            }
            .to_vec(),
            other => {
                let mut message =
                    b"BLAST options error: Unrecognized format specification: '%".to_vec();
                message.push(other);
                message.push(b'\'');
                return Err(ListStop::Ncbi(message, 1));
            }
        };
        data.push(value);
    }
    Ok(format.replace(&data))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1184-1199
// The format is read after the walk, so a format without `%` is an error also for a
// directory without databases (probe J173). NCBI writes each line as it goes; an error of
// NCBI's comes after the lines before it (`CInvalidDataException`: `BLAST options error:`,
// exit status 1; `CSeqDBException`: `BLAST Database error:`, exit status 2), while LOSAT's
// refusal of a database it does not read (an empty `x.pin` on which NCBI crashes: J168)
// writes no line.
/// `-list`: the databases of the directory, one line each, to `out`; the exit status.
pub(super) fn list(args: &Args, path: &SearchPath, diag: &mut Diag, out: &mut Ostream<'_>) -> i32 {
    let dir = args.value("list").unwrap_or("");
    let dbtype = args.value("dbtype").unwrap_or("guess");
    let dbs = find_blast_dbs(
        dir,
        dbtype,
        args.given("recursive"),
        args.given("remove_redundant_dbs"),
    );
    let fmt = args.value("list_outfmt").unwrap_or("%f %p");
    let Ok(format) = ListFormat::parse(fmt.as_bytes()) else {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
        diag.post("BLAST options error: Invalid format specification");
        return 1;
    };
    let mut text = Vec::new();
    for db in &dbs {
        match write_line(&format, db, path) {
            Ok(line) => {
                text.extend_from_slice(&line);
                text.push(b'\n');
            }
            Err(ListStop::Ncbi(message, status)) => {
                let _ = out.write_all(&text);
                let _ = out.flush();
                // `cerr`'s tie flushed the output first: a closed pipe ends the program there.
                if out.closed() {
                    return 0;
                }
                diag.post_bytes(&message);
                return status;
            }
            Err(ListStop::Refused(reason)) => return diag.refuse(&reason),
        }
    }
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does the run's
    // `ncbi_app::files::Ostream` that `out` is.
    let _ = out.write_all(&text);
    let _ = out.flush();
    0
}

#[cfg(test)]
#[path = "list_tests.rs"]
mod tests;
