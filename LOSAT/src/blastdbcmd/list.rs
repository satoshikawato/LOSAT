//! `losat blastdbcmd -list DIR`: the databases of a directory (`FindBlastDBs`) written with
//! `-list_outfmt` (`CBlastDbFormatter`).

use crate::wording::{DATABASE_ERROR, OPTIONS_ERROR};
use std::io::Write;
use std::path::Path;

use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::Ostream;
use crate::seqdb::find::{find_blast_dbs, FoundDb};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

use super::{refusal, Diag};

/// `CBlastDbFormatter`: the format with `%%` made `%`, and the place and letter of each
/// replacement.
#[derive(Debug, PartialEq, Eq)]
struct ListFormat {
    spec: Vec<u8>,
    offsets: Vec<usize>,
    types: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_formatter.cpp:42-64
// After the erase the loop steps over the second `%` and the byte after it (a step of `i` in its
// body and the loop's own); a `%` at the end reads the string's NUL as its letter.
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
            Some(text) => ListStop::Ncbi([DATABASE_ERROR.as_bytes(), b" ", &text].concat(), 2),
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
fn file_length(path: &Path) -> Option<u64> {
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
    use crate::utils::os_bytes::{path_bytes, path_of_bytes};
    let name = db.quoted.as_slice();
    let (volumes, _) = SeqDb::volume_and_alias_paths(name, db.seqtype, path)?;
    let mol = match db.seqtype {
        SeqType::Protein => 'p',
        SeqType::Nucleotide => 'n',
    };
    let extensions = file_extensions(version5);
    let mut total = 0u64;
    let mut dirname = Vec::new();
    for volume in &volumes {
        for ext in &extensions {
            let file = [volume.as_slice(), format!(".{mol}{ext}").as_bytes()].concat();
            let file = path_of_bytes(&file);
            if let Some(length) = file_length(&file) {
                total += length;
                dirname = file
                    .parent()
                    .map(|dir| path_bytes(dir).into_owned())
                    .unwrap_or_default();
            }
        }
    }
    if volumes.len() > 1 {
        if name
            .split(|&byte| byte == b' ')
            .filter(|piece| !piece.is_empty())
            .count()
            > 1
        {
            return Err(ListStop::Ncbi(
                [
                    DATABASE_ERROR.as_bytes(),
                    b" Cannot compute disk usage for multiple BLASTDBs (i.e.: '",
                    name,
                    b"') at once. Please try again using one BLASTDB at a time.",
                ]
                .concat(),
                2,
            ));
        }
        for ext in &extensions {
            let file_name = [name, format!(".{mol}{ext}").as_bytes()].concat();
            let file = path_of_bytes(&dirname).join(path_of_bytes(&file_name));
            if let Some(length) = file_length(&file) {
                total += length;
            }
        }
    }
    Ok(total)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_formatter.cpp:68-162
// The database opens at its first letter that reads it (`CSeqDB` of the name in its quotes and
// `type`: one name, looked for on the search path as written), so an unknown letter after one
// that opens a damaged database is that database's error.
/// `x_InitBlastDb`: the database `db`, opened at the first call.
fn open<'a>(
    opened: &'a mut Option<SeqDb>,
    db: &FoundDb,
    path: &SearchPath,
) -> Result<&'a SeqDb, ListStop> {
    if opened.is_none() {
        *opened = Some(SeqDb::open_in(&db.quoted, db.seqtype, path)?);
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
                    format!("{OPTIONS_ERROR} Unrecognized format specification: '%").into_bytes();
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
    let dir = args.bytes("list").unwrap_or_default();
    let dbtype = args.value("dbtype").unwrap_or("guess");
    let dbs = find_blast_dbs(
        &dir,
        dbtype,
        args.given("recursive"),
        args.given("remove_redundant_dbs"),
    );
    let fmt = args
        .bytes("list_outfmt")
        .unwrap_or(std::borrow::Cow::Borrowed(b"%f %p"));
    let Ok(format) = ListFormat::parse(&fmt) else {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
        diag.post(&format!("{OPTIONS_ERROR} Invalid format specification"));
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
