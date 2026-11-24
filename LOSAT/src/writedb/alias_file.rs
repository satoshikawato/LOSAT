//! NCBI's writer of alias files (`CWriteDB_CreateAliasFile` and `CWriteDB_CreateOidMaskDB`,
//! c++/src/objtools/blast/seqdb_writer/writedb.cpp), which `blastdb_aliastool` and
//! `makeprofiledb` call, and whose files the reader reads (`seqdb::alias`): the header and its
//! date, `TITLE`, `DBLIST`, the key of a list (`GILIST`, `SEQIDLIST`, `TAXIDLIST`) or of an OID
//! mask (`OIDLIST`, `OID_MASK_TYPE`), and the counts `NSEQ` and `LENGTH` that the reader finds
//! when it opens the new alias file. A write that fails is not reported, as NCBI's `ofstream`
//! does not report it: the reader's open of the new file then fails.

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;

use super::dir_entry;
use super::time::LocalTime;
use crate::seqdb::alias::path_base_name;
use crate::seqdb::oid_mask;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

/// The key of the list that an alias file names (`EAliasFileFilterType`): its text is also the
/// name in the error of a list that keeps no sequence. `None` is makeprofiledb's
/// (`eNoAliasFilterType`). blastdb_aliastool 2.17.0 has no `-process_as_tis` (`NCBI_TI` is not
/// defined), so `TILIST` is never written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FilterType {
    None,
    GiList,
    SeqIdList,
    TaxIdList,
}

impl FilterType {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:342-355
    /// `s_AliasFileFilterTypeToString`: the key of the list.
    fn key(self) -> &'static str {
        match self {
            FilterType::None => "",
            FilterType::GiList => "GILIST",
            FilterType::SeqIdList => "SEQIDLIST",
            FilterType::TaxIdList => "TAXIDLIST",
        }
    }
}

/// An alias file that was not made.
#[derive(Debug)]
pub(crate) enum AliasError {
    /// A `CSeqDBException` of the writer: its text, which the programs write after
    /// `BLAST Database error:`.
    SeqDb(String),
    /// The reader's error where NCBI's reader throws it out of the writer (the open of an OID
    /// mask database), or LOSAT's refusal of a database that NCBI's reader opens: the caller
    /// reports it as its other reader errors.
    Reader(SeqDbError),
    /// LOSAT's error: the clock cannot be read for the header's date.
    Clock,
}

/// The message of `AliasError::Clock`.
pub(crate) const CLOCK_ERROR: &str = "the clock cannot be read for the alias file's date";

/// `Created protein ` or `Created nucleotide `, the start of the log line.
fn molecule(seqtype: SeqType) -> &'static str {
    match seqtype {
        SeqType::Protein => "protein ",
        SeqType::Nucleotide => "nucleotide ",
    }
}

/// The extension of an alias file of `seqtype`.
fn alias_extension(seqtype: SeqType) -> &'static str {
    match seqtype {
        SeqType::Protein => ".pal",
        SeqType::Nucleotide => ".nal",
    }
}

/// The header of an alias file made at `now` (the current `CTime` in its default format).
fn header(now: &LocalTime) -> String {
    format!("#\n# Alias file created {}\n#\n", now.log_time())
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:3917-3922
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-377
// `SplitPath` cuts after the last of `ALL_SEPARATORS`: `/` on Unix, `:`, `/` and `\` on Windows.
/// `CDirEntry::GetName` of `path`: what follows the last path separator of `path`.
pub(crate) fn dir_entry_name(path: &str) -> &str {
    let separator = |c: char| c == '/' || (cfg!(windows) && (c == '\\' || c == ':'));
    path.rfind(separator).map_or(path, |at| &path[at + 1..])
}

/// Whether a reader error is one that NCBI's reader throws too (its own text, another exception of
/// NCBI's, or one whose report LOSAT refuses: `SeqDbError::ncbi_throws`): the writer's catch of any
/// exception then takes it for a database that does not open. LOSAT's other refusals go to the
/// caller.
fn ncbi_error(error: &SeqDbError) -> bool {
    error.ncbi_message().is_some() || error.ncbi_exception().is_some() || error.ncbi_throws()
}

/// `NSEQ` and `LENGTH` of the database `name` as NCBI's reader opens it (`CSeqDB` of the name
/// and type: `GetTotals` of `eFilteredAll`, `GetNumSeqs` and `GetTotalLength`). `Ok(None)` is
/// NCBI's failure to open it.
fn totals(
    name: &str,
    seqtype: SeqType,
    path: &SearchPath,
) -> Result<Option<(i32, u64)>, AliasError> {
    match SeqDb::open_in(name, seqtype, path) {
        Ok(db) => Ok(Some((db.ncbi_num_seqs(), db.total_length()))),
        Err(error) if ncbi_error(&error) => Ok(None),
        Err(error) => Err(AliasError::Reader(error)),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:358-414
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:284-328
// The file is written before the reader opens it by `file_name` (a name with a space is several
// databases to the reader, which then fails), and it is removed when that fails or finds no
// sequence. `CDirEntry::Remove` reports nothing.
/// `CWriteDB_CreateAliasFile` of the databases `databases`: the alias file `<file_name>.pal` (or
/// `.nal`) with `title` (none when empty), the names quoted, the list `list_file` under the key of
/// `filter` (none when empty), then its counts; `Ok` is the log line (`LOG_POST`). The reader
/// opens the new file on the search path `path`.
pub(crate) fn create_alias_file(
    file_name: &str,
    databases: &[String],
    seqtype: SeqType,
    list_file: &str,
    title: &str,
    filter: FilterType,
    path: &SearchPath,
) -> Result<String, AliasError> {
    let now = LocalTime::now().ok_or(AliasError::Clock)?;
    let fname = format!("{file_name}{}", alias_extension(seqtype));
    let mut text = header(&now);
    if !title.is_empty() {
        text.push_str(&format!("TITLE {title}\n"));
    }
    text.push_str("DBLIST ");
    for name in databases {
        text.push_str(&format!("\"{name}\" "));
    }
    text.push('\n');
    if !list_file.is_empty() {
        text.push_str(&format!("{} {list_file}\n", filter.key()));
    }
    let _ = std::fs::write(&fname, text);
    let counts = totals(file_name, seqtype, path).inspect_err(|_| {
        dir_entry::remove(Path::new(&fname));
    })?;
    let Some((num_seqs, length)) = counts else {
        dir_entry::remove(Path::new(&fname));
        return Err(AliasError::SeqDb(
            "BLASTDB alias file creation failed.  Some referenced files may be missing".to_string(),
        ));
    };
    if num_seqs == 0 {
        dir_entry::remove(Path::new(&fname));
        return Err(AliasError::SeqDb(format!(
            "No seqs in {} were found in BLAST database",
            filter.key()
        )));
    }
    if let Ok(mut file) = OpenOptions::new().append(true).open(&fname) {
        let _ = write!(file, "NSEQ {num_seqs}\nLENGTH {length}\n");
    }
    Ok(format!(
        "Created {}BLAST (alias) database {file_name} with {num_seqs} sequences",
        molecule(seqtype)
    ))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:420-429
/// `s_GetNumOfDigits`: the digits of `n`, at least 2 (the volume numbers' width).
fn num_digits(n: usize) -> usize {
    n.to_string().len().max(2)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:432-474
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:257-271
// The vector of names starts with `num_volumes` empty names, to which the volumes are added
// (`push_back`), so `DBLIST` begins with one space for each. The volumes' names are not
// quoted. The counts are those of the names joined by spaces (one after each name), and are 0
// when they do not open; the log names that joined text.
/// `CWriteDB_CreateAliasFile` of `num_volumes` volumes `<file_name>.NN`: each must be found on
/// the search path `path` (NCBI's error names the first that is not), then the alias file
/// `<file_name>.pal` (or `.nal`) of their names, its title and counts; `Ok` is the log line.
pub(crate) fn create_volumes_alias_file(
    file_name: &str,
    num_volumes: u32,
    seqtype: SeqType,
    title: &str,
    path: &SearchPath,
) -> Result<String, AliasError> {
    let width = num_digits(num_volumes as usize);
    let mut volume_names = Vec::with_capacity(num_volumes as usize);
    let mut concatenated = String::new();
    for i in 0..num_volumes {
        let vol_name = format!("{file_name}.{i:0width$}");
        if path.resolve_db(&vol_name, seqtype).is_none() {
            return Err(AliasError::SeqDb(format!(
                "Failed to find {}{vol_name} BLAST database",
                molecule(seqtype)
            )));
        }
        concatenated.push_str(&vol_name);
        concatenated.push(' ');
        volume_names.push(vol_name);
    }
    let (num_seqs, length) = totals(&concatenated, seqtype, path)?.unwrap_or((0, 0));
    let now = LocalTime::now().ok_or(AliasError::Clock)?;
    let mut text = header(&now);
    if !title.is_empty() {
        text.push_str(&format!("TITLE {title}\n"));
    }
    text.push_str("DBLIST ");
    text.push_str(&" ".repeat(num_volumes as usize));
    for name in &volume_names {
        text.push_str(dir_entry_name(name));
        text.push(' ');
    }
    text.push('\n');
    text.push_str(&format!("NSEQ {num_seqs}\nLENGTH {length}\n"));
    let fname = format!("{file_name}{}", alias_extension(seqtype));
    let _ = std::fs::write(fname, text);
    Ok(format!(
        "Created {}BLAST (alias) database {concatenated} with {num_seqs} sequences",
        molecule(seqtype)
    ))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:548-624
// NCBI's errors are `CSeqDBException`s: the input that does not open (any exception), the volume
// without its mask file (after `TITLE` is written, and the earlier volumes' files). The open of the
// new database is not caught: its errors are the reader's. Each line ends with `endl`. The
// `DBLIST` of several volumes names `<output_db>.NN` as given, which the reader resolves.
/// `CWriteDB_CreateOidMaskDB`: the alias file `<output_db>.pal` (or `.nal`) of the volumes of
/// `input_db` with the OID masks of `mask` (`OID_MASK_TYPE`; with its exclude-model bit, each
/// volume's mask file under `OIDLIST`), one alias file `<output_db>.NN.pal` for each volume when
/// there are several, then the counts of the new database.
pub(crate) fn create_oid_mask_db(
    input_db: &str,
    output_db: &str,
    seqtype: SeqType,
    mask: i32,
    title: &str,
    path: &SearchPath,
) -> Result<(), AliasError> {
    match SeqDb::open_in(input_db, seqtype, path) {
        Ok(_) => {}
        Err(error) if ncbi_error(&error) => {
            return Err(AliasError::SeqDb("Invalid input db".to_string()))
        }
        Err(error) => return Err(AliasError::Reader(error)),
    }
    let (vols, _) =
        SeqDb::volume_and_alias_paths(input_db, seqtype, path).map_err(AliasError::Reader)?;
    if vols.is_empty() {
        return Err(AliasError::SeqDb(format!("no vol found for {input_db}")));
    }
    let ext = alias_extension(seqtype);
    let width = num_digits(vols.len());
    let mut ofs = File::create(format!("{output_db}{ext}")).ok();
    let write = |file: &mut Option<File>, text: &str| {
        if let Some(file) = file {
            let _ = file.write_all(text.as_bytes());
        }
    };
    write(&mut ofs, &format!("TITLE {title}\n"));
    for (i, vol) in vols.iter().enumerate() {
        let mut lines = format!("DBLIST {}\nOIDLIST ", path_base_name(vol));
        // NCBI tests the exclude-model bit of the mask as an `int`: a negative value has the bit.
        if (mask as u64) & oid_mask::EXCLUDE_MODEL != 0 {
            let full_path = format!(
                "{vol}.{}",
                oid_mask::exclude_model_extension(seqtype == SeqType::Protein)
            );
            // `CFile::Exists`: a regular file, also through a link.
            if !Path::new(&full_path).is_file() {
                return Err(AliasError::SeqDb(format!(
                    "Exclude oid mask file not found for {vol}"
                )));
            }
            lines.push_str(dir_entry_name(&full_path));
        }
        lines.push_str(&format!("\nOID_MASK_TYPE {mask}\n"));
        if vols.len() > 1 {
            let mut ovs = File::create(format!("{output_db}.{i:0width$}{ext}")).ok();
            write(&mut ovs, &lines);
        } else {
            write(&mut ofs, &lines);
        }
    }
    if vols.len() > 1 {
        let mut dblist = String::from("DBLIST");
        for i in 0..vols.len() {
            dblist.push_str(&format!(" {output_db}.{i:0width$}"));
        }
        dblist.push('\n');
        write(&mut ofs, &dblist);
    }
    let db = SeqDb::open_in(output_db, seqtype, path).map_err(AliasError::Reader)?;
    write(
        &mut ofs,
        &format!(
            "NSEQ {}\nLENGTH {}\n",
            db.ncbi_num_seqs(),
            db.total_length()
        ),
    );
    Ok(())
}

#[cfg(test)]
#[path = "alias_file_tests.rs"]
mod tests;
