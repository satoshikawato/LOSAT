//! A database build (`CBuildDatabase`,
//! `c++/src/objtools/blast/seqdb_writer/build_db.cpp`): the log of
//! makeblastdb, the removal of an existing database of the same name, the
//! loop that hands the sequences of a FASTA input to the writer, and the end
//! of the build, which removes the files of a failed one.

use std::collections::BTreeSet;
use std::io::{Read, Write};
use std::path::Path;

use super::bioseq::{Bioseq, DeltaFix};
use super::deflines;
use super::dir_entry;
use super::fasta_source::{FastaBioseqSource, SourceErrorKind, Warn};
use super::mask_data::range_set::address_order_refusal;
use super::mask_data::{ranges_are_empty, MaskedRangeSet};
use super::platform::Platform;
use super::raw_source::{LeafTaxIds, RawCopyError, RawSeqDbSource};
use super::taxid_set::{seq_id_key, TaxIdSet};
use super::time::LocalTime;
use super::writer::{DbVersion, Indices, WriteDb};
use super::WriteDbError;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_input::Utf8TagOnce;
use crate::seqdb::column::MASK_DATA_TITLE;
use crate::seqdb::{SeqDbError, SeqType};
use crate::seqid::model::ObjectId;
use crate::seqid::rank::{best_rank_score, find_best};
use crate::seqid::{DeflineSet, ParsedSeqId};

/// `CBuildDatabase`: a FASTA input (`add_fasta`) or the sequences of existing databases
/// (`add_raw_sequences`).
pub struct BuildDatabase<'a> {
    is_protein: bool,
    log: &'a mut dyn Write,
    output_db: WriteDb,
    output_dbname: String,
    /// `m_ParseIDs` (`-parse_seqids`).
    parse_ids: bool,
    /// `m_Taxids` (`SetTaxids`).
    taxids: TaxIdSet,
    /// `m_MaskData`: the masked ranges of the sequences, when makeblastdb read mask data
    /// (`SetMaskDataSource`).
    mask_data: Option<MaskedRangeSet>,
    /// `m_FoundMatchingMasks`.
    found_matching_masks: bool,
    /// LOSAT's refusal of the first sequence whose masks NCBI writes in the order of memory
    /// addresses (`mask_data::range_set`). It is raised by `end_build` only when the build would
    /// otherwise succeed, so that every error of NCBI that comes later in the run comes first.
    address_order: Option<WriteDbError>,
    /// `m_Id2Leafs` (`SetLeafTaxIds`, with `m_KeepLeafs` set): the leaf taxids of a copied
    /// database's definition lines by the first Seq-id's `AsFastaString`.
    leaf_tax_ids: LeafTaxIds,
    /// `ERR_POST_X_ONCE` of the run's warning 803.10 of a UTF8String tag (one run is one NCBI
    /// process): the mask files (`mask_data::process`) and the input database's headers
    /// (`raw_source::scan_leaf_tax_ids`) post it once between them.
    utf8_once: Utf8TagOnce,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:852-927
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-824
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:951-960
/// `CDirEntry::NormalizePath` of a Unix path, without following links: `.`
/// and empty components go, `..` takes its parent with it. A relative path
/// stays relative (`.` when nothing is left).
#[cfg(not(windows))]
fn normalize_path(path: &str) -> String {
    if path.is_empty() {
        return String::new();
    }
    let current = path.trim_end_matches('/');
    if current.is_empty() {
        return "/".to_string();
    }
    let mut head: Vec<&str> = Vec::new();
    for next in current.split('/') {
        if let Some(&last) = head.last() {
            if last == "." {
                if !next.is_empty() {
                    head.pop();
                }
            } else if next == "." || next.is_empty() {
                continue;
            } else if next == ".." {
                if last.is_empty() {
                    continue;
                } else if last != ".." {
                    head.pop();
                    continue;
                }
            }
        }
        head.push(next);
    }
    if head.is_empty() || head == [".", ""] {
        return ".".to_string();
    }
    if head.len() == 1 && head[0].is_empty() {
        return "/".to_string();
    }
    head.join("/")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:665-670
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:693-696
// The reference platform is Unix. On Windows the path is made absolute by
// the system's own rule. The WASI targets take Unix's rule: their current
// directory is the C library's (wasi-libc's `getcwd`), which the command sets
// from the `PWD` that its host gives (`cli::main_entry`) and which is `/`
// when the host gives none.
/// `CDirEntry::CreateAbsolutePath` of `dbname`.
fn create_absolute_path(dbname: &str) -> String {
    // `NormalizePath` leaves no separator at the end of a path that is not a
    // root, so `-out sub/` names the database `sub` on every platform.
    #[cfg(windows)]
    {
        let absolute = std::path::absolute(dbname)
            .map(|path| path.to_string_lossy().into_owned())
            .unwrap_or_else(|_| dbname.to_string());
        dir_entry_path(&absolute).to_string()
    }
    #[cfg(not(windows))]
    {
        absolute_path_in(std::env::current_dir().ok().as_deref(), dbname)
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:776-788
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3708-3716
/// `CDirEntry::CreateAbsolutePath` of `dbname` on Unix, where the current
/// directory is `cwd` (`None` when it cannot be read: `GetCwd` gives "").
/// `ConcatPath` cuts the C spaces at the ends of both parts (a tab at the end
/// of `-out` is not in the database's name; makeblastdb refuses a name with a
/// space before); the result is normalized, relative when there is no
/// current directory.
#[cfg(not(windows))]
fn absolute_path_in(cwd: Option<&Path>, dbname: &str) -> String {
    if dbname.starts_with('/') {
        return normalize_path(dbname);
    }
    let cwd = cwd.map(|cwd| cwd.to_string_lossy()).unwrap_or_default();
    let mut path = crate::ncbi_str::truncate_spaces_str(&cwd).to_string();
    if !path.is_empty() && !path.ends_with('/') {
        path.push('/');
    }
    let mut part = crate::ncbi_str::truncate_spaces_str(dbname);
    if !path.is_empty() {
        part = part.strip_prefix('/').unwrap_or(part);
    }
    path.push_str(part);
    normalize_path(&path)
}

/// A separator of directory names (`DIR_SEPARATORS`).
fn is_separator(c: char) -> bool {
    c == '/' || (cfg!(windows) && c == '\\')
}

/// A character that ends the directory part of a path (`ALL_SEPARATORS`).
fn is_any_separator(c: char) -> bool {
    is_separator(c) || (cfg!(windows) && c == ':')
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:298-313
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:465-472
/// The path a `CDirEntry` keeps for `path`: without its trailing separators,
/// unless it is one separator (the root) or, on Windows, a disk name.
fn dir_entry_path(path: &str) -> &str {
    if path.len() == 1 && path.starts_with(is_separator) {
        return path;
    }
    if cfg!(windows) && (path.len() == 2 || path.len() == 3) && path.as_bytes()[1] == b':' {
        return path;
    }
    path.trim_end_matches(is_separator)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:4021-4046,4056-4067
// The mode is the default of a `CDir`, rwxr-xr-x. `[NCBI] FileAPIHonorUmask`
// is false unless the registry sets it (ncbifile.cpp:161; LOSAT refuses a
// registry that does, `crate::blastinput::ncbi_environment`), so the new
// directory gets exactly that mode, whatever the umask.
/// `s_DirCreate` of `path` with `fCreate_Default` and `mode`: makes one directory; `true`
/// also when it is there already.
fn dir_create(path: &str) -> bool {
    match std::fs::create_dir(path) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
            // `GetType()` does not follow a link: a link to a directory
            // is not a directory here.
            return std::fs::symlink_metadata(path).is_ok_and(|meta| meta.is_dir());
        }
        Err(_) => return false,
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).is_ok()
    }
    #[cfg(not(unix))]
    {
        true
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:4106-4139,4165-4186
// `CDirEntry::Exists` is true for an entry of any type, a link that leads
// nowhere included (`GetType()` without following links).
/// `CreatePath` of the `CDir` of `dir`: makes the directory `dir` and its missing
/// parents, after the path was made absolute and normalized without
/// following links (`a/b/../c` is `a/c`: `a/b` is not made).
fn create_path(dir: &str) -> bool {
    if dir.is_empty() {
        return false;
    }
    let mut path = create_absolute_path(dir);
    if path.is_empty() {
        return false;
    }
    if path.ends_with(is_separator) && !(cfg!(windows) && path.len() == 3) {
        path.pop();
    }
    let mut tmp = path.as_str();
    let mut missed_parts: Vec<&str> = Vec::new();
    while !tmp.is_empty() && std::fs::symlink_metadata(tmp).is_err() {
        let Some(pos) = tmp.rfind(is_separator) else {
            break;
        };
        missed_parts.push(&tmp[pos + 1..]);
        tmp = &tmp[..pos];
    }
    if missed_parts.is_empty() {
        return dir_create(&path);
    }
    let mut p = tmp.to_string();
    for part in missed_parts.iter().rev() {
        p.push(std::path::MAIN_SEPARATOR);
        p.push_str(part);
        if !dir_create(&p) {
            return false;
        }
    }
    true
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1051-1071
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-368
// The directory is the part of `dbname` up to its last separator, after the
// trailing separators of `dbname` are dropped: `-out sub/` names the
// database `sub` of the current directory and makes no directory.
// `d.GetName()` is what follows the last separator of the directory's path
// as it was given (`..` for `nodir/..`), and `d.Exists()` is true for a
// directory, also through a link.
//
// `CreatePath` makes the normalized path, but `CheckAccess` then tests the
// path as it was given: with `-out nodir/../db` and no `nodir`, the current
// directory "is created" (it exists), `nodir/..` does not resolve, and
// makeblastdb stops with the write-permission message.
//
// On Unix LOSAT asks the question of `euidaccess` for write access to `path`
// (`check_access_write`): a path that does not resolve and a directory that exists and cannot be
// written to (by its mode, an ACL, a read-only file system) both stop it here, before the log. On
// Windows NCBI reads the effective rights of the directory's ACL
// (`CWinSecurity::GetFilePermissions`), which LOSAT does not: it stops here only for a path that
// does not resolve, and for a directory it cannot write to it goes on, writes the log header and
// fails when it creates the first file of the first volume, with the system's message for that file
// and the same exit status 255. The library is built with `-DNCBI_MODULE=BLASTDB`
// (Makefile.writedb.lib:15), so the report names the function
// `BLASTDB::ncbi::CBuildDatabase::CreateDirectories`; the error code `eOutputFileError` has no name
// of its own (`CMultisourceException::GetErrCodeString` names `eArg` only,
// multisource_util.hpp:66-72) and prints as `eInvalid`.
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/multisource_util.hpp:65-72
const CREATE_PATH_FAILED: ThrowSite = create_directories_site(1063);
/// `s_FixBioseqDeltas`' throw for a protein delta (`CBuildDatabase::delta_fix`).
const PROTEIN_DELTA: ThrowSite = ThrowSite {
    file: "src/objtools/blast/seqdb_writer/build_db.cpp",
    line: 303,
    kind: "CMultisourceException::eArgErr",
    function: Some("BLASTDB::ncbi::s_FixBioseqDeltas"),
};
const NO_WRITE_ACCESS: ThrowSite = create_directories_site(1069);

/// The throws of `AddSequences` of an `IRawSequenceSource` for raw data it cannot add
/// (`CMultisourceException::eArg`).
const fn raw_data_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/objtools/blast/seqdb_writer/build_db.cpp",
        line,
        kind: "CMultisourceException::eArgErr",
        function: Some("BLASTDB::ncbi::CBuildDatabase::AddSequences"),
    }
}

/// The exception of `AddSequences` of an `IRawSequenceSource` for raw data it cannot add, thrown at
/// `line` of build_db.cpp.
fn raw_data_error(line: u32, message: &str) -> WriteDbError {
    WriteDbError::Exception(NcbiException::new(raw_data_site(line), message))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:418-426
/// `x_AddPig`: the PIG of a copied sequence, the first value of its first definition line's
/// `other-info` (0 without one). An `other-info` without values, whose first value NCBI reads
/// all the same, is LOSAT's refusal.
fn raw_pig(headers: &DeflineSet) -> Result<i32, WriteDbError> {
    let Some(first) = headers.deflines.first() else {
        return Ok(0);
    };
    match &first.other_info {
        None => Ok(0),
        Some(values) => values.first().copied().ok_or_else(|| {
            let mut id = Vec::new();
            if let Some(seq_id) = first.ids.first() {
                crate::seqid::label::fasta(seq_id, &mut id);
            }
            WriteDbError::Unsupported(format!(
                "the first definition line of the input sequence {} has an empty other-info, \
                 whose first value NCBI makeblastdb reads as the PIG (no defined result)",
                String::from_utf8_lossy(&id)
            ))
        }),
    }
}

/// A throw of `CMultisourceException` (`eOutputFileError`) at `line` of build_db.cpp.
const fn create_directories_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/objtools/blast/seqdb_writer/build_db.cpp",
        line,
        kind: "CMultisourceException::eInvalid",
        function: Some("BLASTDB::ncbi::CBuildDatabase::CreateDirectories"),
    }
}

/// `CBuildDatabase::CreateDirectories`.
pub(crate) fn create_directories(dbname: &str) -> Result<(), WriteDbError> {
    let dir_entry = dir_entry_path(dbname);
    let Some(at) = dir_entry.rfind(is_any_separator) else {
        return Ok(());
    };
    let d = dir_entry_path(&dir_entry[..=at]);
    let name = d.rfind(is_any_separator).map_or(d, |at| &d[at + 1..]);
    if !Path::new(d).is_dir() && !create_path(d) {
        return Err(WriteDbError::Exception(NcbiException::new(
            CREATE_PATH_FAILED,
            format!("Failed to create directory '{name}'"),
        )));
    }
    if !check_access_write(d) {
        return Err(WriteDbError::Exception(NcbiException::new(
            NO_WRITE_ACCESS,
            format!("You do not have write permissions on '{name}'"),
        )));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:1735-1748
/// `CDirEntry::CheckAccess` for `fWrite` of the path as it was given: `euidaccess` of
/// the path with `F_OK` and `W_OK`. glibc's `euidaccess` calls `access` with the same mode when
/// the process's real and effective user and group ids are equal, and LOSAT makes that call (`F_OK`
/// is 0). For a set-user-ID or set-group-ID process glibc compares the mode bits with the effective
/// ids instead; `access` asks with the real ids there.
#[cfg(unix)]
fn check_access_write(path: &str) -> bool {
    rustix::fs::access(path, rustix::fs::Access::WRITE_OK).is_ok()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:1722-1732
/// `CDirEntry::CheckAccess` for `fWrite` without the rights of the ACL (see
/// `create_directories`): only whether the path resolves.
#[cfg(not(unix))]
fn check_access_write(path: &str) -> bool {
    std::fs::metadata(path).is_ok()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2629-2641
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2647-2656
/// `SeqDB_GetFileExtensions` for version 4 (`eBDB_Version4`) without the
/// letter of the molecule.
pub(crate) const VOLUME_EXTENSIONS: [&str; 20] = [
    "al", "in", "hr", "sq", "ni", "nd", "si", "sd", "pi", "pd", "js", "aa", "ab", "ac", "og", "hi",
    "hd", "ti", "td", "xm",
];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2663-2668
/// `SeqDB_GetLMDBFileExtensions` without the letter of the molecule.
const LMDB_EXTENSIONS: [&str; 7] = ["db", "os", "ot", "tf", "to", "db-lock", "tf-lock"];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1542-1593
// NCBI finds the volumes of the existing database with its reader and
// ignores a database it cannot read. LOSAT's reader refuses some alias files
// that NCBI reads (filters, replaced totals, group alias files); leaving
// such an alias file in place would hide the new database behind it, so an
// existing database of that kind stops the build before anything is removed.
//
// Two departures from `DeleteBlastDb`, both in how the existing database is
// found (`crate::seqdb::alias::find_volume_paths`):
//
// * `dbpath` is one database name. NCBI's reader splits the name at spaces
//   and double quotes, so for an output path with such a character (a
//   directory `my data`) the call above throws: only the LMDB files removed
//   before it are gone, the volume files stay, the new files are written over
//   or beside them (an old alias file then hides the new database), and
//   makeblastdb ends with exit status 2 when it opens its result. LOSAT
//   removes the whole existing database of that path and the build ends with
//   status 0.
// * The volumes that an alias file names are those in the alias file's own
//   directory. NCBI takes a volume of the same name in the current directory
//   first, and so removes files that do not belong to the database being
//   replaced.
/// `DeleteBlastDb`: removes the files of the database `dbpath` (and a directory
/// of one of their names, with its contents); `true` when one was removed.
pub(crate) fn delete_blast_db(dbpath: &str, is_protein: bool) -> Result<bool, WriteDbError> {
    let letter = if is_protein { 'p' } else { 'n' };
    let seqtype = if is_protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let (db_files, alias_files) = match crate::seqdb::alias::find_volume_paths(dbpath, seqtype) {
        Ok(found) => found,
        Err(SeqDbError::Unsupported { path, message }) => {
            return Err(WriteDbError::Other(format!(
                "an existing database is in the way and LOSAT cannot list its files ({}: \
                 {message}); remove it before building {dbpath}",
                path.display()
            )));
        }
        Err(_) => (Vec::new(), Vec::new()),
    };

    // `Remove` of a `CFile` is `CDirEntry::Remove` with `eRecursive` (`super::dir_entry`): a
    // directory of the name is removed with its contents, a symbolic link is removed and its
    // target kept.
    let mut num_files_removed = 0;
    let mut remove = |name: &str| {
        let removed = dir_entry::remove(Path::new(name));
        num_files_removed += u32::from(removed);
        removed
    };
    for lmdb in LMDB_EXTENSIONS {
        if !remove(&format!("{dbpath}.{letter}{lmdb}")) {
            let mut index = 0u32;
            // The removals of the volume files are not counted by NCBI.
            while dir_entry::remove(Path::new(&format!(
                "{dbpath}.{}{}.{letter}{lmdb}",
                index / 10,
                index % 10
            ))) {
                index += 1;
            }
        }
    }
    for file in &db_files {
        for extension in VOLUME_EXTENSIONS {
            remove(&format!("{file}.{letter}{extension}"));
        }
    }
    for file in &alias_files {
        remove(file);
    }
    Ok(num_files_removed != 0)
}

/// The stream output of a `double` with the default format (six significant digits,
/// `%g`): the seconds of the "Adding sequences" line.
fn format_seconds(value: f64) -> String {
    if value == 0.0 || !value.is_finite() {
        return "0".to_string();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let strip = |text: &str| {
        if text.contains('.') {
            text.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            text.to_string()
        }
    };
    if !(-4..6).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        return format!("{}e{sign}{:02}", strip(mantissa), exponent.abs());
    }
    strip(&format!("{value:.*}", (5 - exponent) as usize))
}

impl<'a> BuildDatabase<'a> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1165-1177
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1270-1273
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1094
    /// The constructor, then `SetMembBits` with no bits and the lines' own kept, as makeblastdb
    /// calls it. `now` is the time of the build, `version` the format version
    /// (`dbver`), `platform` the platform whose makeblastdb the files follow,
    /// `indices` `indexing` (`eFullIndex` for `-parse_seqids`, `eAddHash` for
    /// `-hash_index`). `title` is bytes: the title of an input database (`-input_type
    /// blastdb`) is copied as it is stored.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dbname: &str,
        title: &[u8],
        is_protein: bool,
        log: &'a mut dyn Write,
        now: &LocalTime,
        version: DbVersion,
        platform: Platform,
        indices: Indices,
    ) -> Result<Self, WriteDbError> {
        create_directories(dbname)?;
        let output_dbname = create_absolute_path(dbname);
        let mol_type = if is_protein { "Protein" } else { "Nucleotide" };
        let _ = writeln!(
            log,
            "\n\nBuilding a new DB, current time: {}",
            now.log_time()
        );
        let _ = writeln!(log, "New DB name:   {output_dbname}");
        let _ = log.write_all(&[b"New DB title:  ".as_slice(), title, b"\n"].concat());
        let _ = writeln!(log, "Sequence type: {mol_type}");
        if delete_blast_db(&output_dbname, is_protein)? {
            let _ = writeln!(
                log,
                "Deleted existing {mol_type} BLAST database named {output_dbname}"
            );
        }
        let _ = writeln!(log, "Keep MBits: T");
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1119-1132
        // The limit of 1000*1000*1000 that the constructor sets is replaced
        // by makeblastdb before a sequence is added.
        let output_db = WriteDb::new(
            &output_dbname,
            is_protein,
            title,
            now.index_date().as_bytes(),
            1000 * 1000 * 1000,
            version,
            platform,
            indices,
        );
        let mut output_db = output_db;
        output_db.set_column_date_fallback(now.log_time().as_bytes());
        Ok(Self {
            is_protein,
            log,
            output_db,
            output_dbname,
            parse_ids: indices.full,
            taxids: TaxIdSet::default(),
            mask_data: None,
            found_matching_masks: false,
            address_order: None,
            leaf_tax_ids: LeafTaxIds::new(),
            utf8_once: Utf8TagOnce::default(),
        })
    }

    /// `oid_masks` of the `CBuildDatabase` constructor (`-oid_masks`): the OID mask files the
    /// volumes write (`writedb::oid_list`).
    pub fn set_oid_masks(&mut self, oid_masks: u64) {
        self.output_db.set_oid_masks(oid_masks);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1216-1219
    /// `SetTaxids`: the taxids of `-taxid` or `-taxid_map`.
    pub fn set_taxids(&mut self, taxids: TaxIdSet) {
        self.taxids = taxids;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1211-1213
    /// The message of `~CBuildDatabase` about the taxids: `Some` when `-taxid_map` matched no
    /// sequence.
    pub fn taxid_message(&self) -> Option<&'static str> {
        (!self.taxids.has_ever_fixed_id())
            .then_some("Error: [makeblastdb] No sequences matched any of the taxids provided.")
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1583-1594
    /// `RegisterMaskingAlgorithm` of a program, its options and a name: see
    /// [`WriteDb::register_mask_algorithm`].
    pub fn register_masking_algorithm(
        &mut self,
        program: i32,
        options: &[u8],
        name: &[u8],
    ) -> Result<i32, WriteDbError> {
        self.output_db
            .register_mask_algorithm(program, options, name)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1596-1607
    // The first argument of this overload is the id of `-mask_id` (the call
    // in makeblastdb.cpp:966 passes `id_list[i], desc_list[i], opts`).
    /// `RegisterMaskingAlgorithm` of an id, a description and options: see
    /// [`WriteDb::register_mask_algorithm_id`].
    pub fn register_masking_algorithm_id(
        &mut self,
        id: &[u8],
        description: &[u8],
        options: &[u8],
    ) -> Result<i32, WriteDbError> {
        self.output_db
            .register_mask_algorithm_id(id, description, options)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1609-1612
    // The masks are looked up by the Seq-ids of each sequence (`bs->GetId()`): with
    // `-parse_seqids` the parsed ids, else its `local id` (`add_masks_for_seq_id`).
    /// `SetMaskDataSource`: the masked ranges of the sequences. makeblastdb calls it when it
    /// reads the first mask list that has masks, with an empty set, and goes on inserting into
    /// it through [`Self::mask_data_source_mut`] (`process_mask_data` does the same); every
    /// insert comes before the first sequence is added.
    pub fn set_mask_data_source(&mut self, ranges: MaskedRangeSet) {
        self.mask_data = Some(ranges);
    }

    /// The set that `set_mask_data_source` handed over, to insert into.
    pub fn mask_data_source_mut(&mut self) -> Option<&mut MaskedRangeSet> {
        self.mask_data.as_mut()
    }

    /// `m_FoundMatchingMasks`: whether a sequence has been given mask data.
    #[cfg(test)]
    pub fn found_matching_masks(&self) -> bool {
        self.found_matching_masks
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1204-1210
    // The message is posted when the `CBuildDatabase` is destroyed: in makeblastdb, when the
    // application object goes (`m_DB` is a member), so after the log of `EndBuild`, after the
    // exception message of a failed build and after the metadata file. The caller of LOSAT's
    // build calls this at that point.
    /// The destructor's diagnostic: `warn` gets the message (without its final newline)
    /// when mask data were read and no sequence matched.
    pub fn post_destruction_diagnostics(&self, warn: Warn<'_>) {
        if self.mask_data.is_some() && !self.found_matching_masks {
            warn(
                b"Error: [makeblastdb] No sequences matched any of the masks provided.\n\
                  Please ensure that the -parse_seqids option is used in the\n\
                  filtering program as well as makeblastdb.",
            );
        }
    }

    /// Test hook: see [`WriteDb::set_column_date`].
    #[cfg(test)]
    pub(crate) fn set_column_date(&mut self, index: usize, date: &[u8]) {
        self.output_db.set_column_date(index, date);
    }

    /// `GetOutputDbName`: the absolute name of the database.
    pub fn output_dbname(&self) -> &str {
        &self.output_dbname
    }

    /// The log stream (`m_LogFile`).
    pub fn log(&mut self) -> &mut dyn Write {
        self.log
    }

    /// The database being written (`m_OutputDb`).
    pub fn output_db(&self) -> &WriteDb {
        &self.output_db
    }

    /// `SetMaxFileSize`.
    pub fn set_max_file_size(&mut self, max_file_size: u64) {
        self.output_db.set_max_file_size(max_file_size);
    }

    /// The constructor's `use_gi_mask` (`makeblastdb -gi_mask`): see [`WriteDb::use_gi_masks`].
    /// Called before the masking algorithms are registered.
    pub fn use_gi_masks(&mut self) {
        self.output_db.use_gi_masks();
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:859-876
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-487
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:894-901
    // The FASTA source applies `s_FixBioseqDeltas` before it hands a
    // sequence over; its errors arrive as errors of the source. What it makes of a delta
    // that `Seq-inst.length` does not describe meets the sequence here (`delta_fix`).
    /// `AddSequences` of an `IBioseqSource`: `true` when a sequence was added.
    fn add_sequences<R: Read>(
        &mut self,
        source: &mut FastaBioseqSource<R>,
        warn: Warn<'_>,
    ) -> Result<bool, WriteDbError> {
        let mut found = false;
        let start = std::time::Instant::now();
        let mut count = 0u64;
        while let Some(mut bioseq) = source.next_bioseq(warn).map_err(|error| match error.kind {
            SourceErrorKind::Creation => WriteDbError::creation(error.message),
            SourceErrorKind::OutOfMemory => WriteDbError::OutOfMemory,
            SourceErrorKind::Std => {
                WriteDbError::Other(String::from_utf8_lossy(&error.message).into_owned())
            }
        })? {
            let bioseq_id = bioseq.first_id_fasta();
            if let Some(parsed) = bioseq.parsed.as_mut() {
                if let Some(gi) = deflines::large_gi_to_skip(&parsed.ids, &mut *warn) {
                    let _ = writeln!(
                        self.log,
                        "Ignoring gi '{gi}' as it has value larger then 4294967295"
                    );
                    continue;
                }
                deflines::localize_other_molecule(&mut parsed.ids, self.is_protein, &mut *warn);
            }
            if bioseq.is_aa != self.is_protein {
                continue;
            }
            if bioseq.length > 0x7fff_ffff {
                return Err(WriteDbError::Creation(format!(
                    "Sequences longer than 2,147,483,647 bases are not supported.  Offending \
                     sequence is {bioseq_id}"
                )));
            }
            let mut added = false;
            if bioseq.length != 0 {
                if self.parse_ids {
                    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-477,497
                    // `x_EditHeaders` sets the taxids (`FixTaxId`); its other edits do
                    // nothing for makeblastdb's FASTA input (no GIs removed, no
                    // membership bits or leaf taxids given).
                    let mut headers = deflines::extract_bioseq_deflines(&bioseq, &mut *warn);
                    self.output_db.post_extracted_deflines(&headers);
                    self.taxids.fix_tax_id(&mut headers);
                    if let Some(fix) = bioseq.delta_fix.take() {
                        self.delta_fix(fix, &mut bioseq, &mut *warn)?;
                    }
                    if bioseq.seq_data.is_some() {
                        // `bs->GetId()`: the ids after `AddSequences`' edit of the ones of the
                        // other molecule type (`localize_other_molecule`, above).
                        let ids = self
                            .mask_data
                            .is_some()
                            .then(|| bioseq.parsed.as_ref().map(|parsed| parsed.ids.clone()))
                            .flatten();
                        self.output_db.add_sequence(bioseq)?;
                        self.output_db.set_deflines(headers);
                        added = true;
                        if let Some(ids) = ids {
                            self.add_masks_for_seq_id(&ids, &bioseq_id)?;
                        }
                    }
                } else {
                    // The definition line without `-parse_seqids` gets the global taxid
                    // of `-taxid` (or 0); `-taxid_map` requires `-parse_seqids`.
                    self.output_db.extract_bioseq_deflines(&bioseq);
                    if let Some(fix) = bioseq.delta_fix.take() {
                        self.delta_fix(fix, &mut bioseq, &mut *warn)?;
                    }
                    if bioseq.seq_data.is_some() {
                        // `bs->GetId()` is the `local id N` that `CFastaReader` numbers the
                        // sequences by (`CSeqIdGenerator::GenerateID` without a prefix); a
                        // number that is no `int` matches none.
                        let ids: Vec<ParsedSeqId> = self
                            .mask_data
                            .is_some()
                            .then(|| bioseq.local_id.parse::<i32>().ok())
                            .flatten()
                            .map(|number| ParsedSeqId::Local(ObjectId::Id(number)))
                            .into_iter()
                            .collect();
                        self.output_db.add_sequence(bioseq)?;
                        self.output_db.set_taxid(self.taxids.global_taxid());
                        added = true;
                        self.add_masks_for_seq_id(&ids, &bioseq_id)?;
                    }
                }
            }
            // The serializer posts its diagnostics while the sequence is
            // added, before the reader reads (and warns about) the next
            // record; both go to the standard error stream in that order.
            self.post_diagnostics(warn);
            if !added {
                let _ = writeln!(
                    self.log,
                    "Ignoring sequence '{bioseq_id}' as it has no sequence data"
                );
                continue;
            }
            found = true;
            count += 1;
        }
        if count != 0 {
            let _ = writeln!(
                self.log,
                "Adding sequences from FASTA; added {count} sequences in {} seconds.",
                format_seconds(start.elapsed().as_secs_f64())
            );
        }
        Ok(found)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:446-466
    // The Seq-ids are `bs->GetId()`: the parsed ids of a sequence with `-parse_seqids`, its
    // `local id` without it. A sequence with no 31-bit local id has none. The GIs are those of
    // the first definition line (the Bioseq's ids), in its order; `-gi_mask` requires
    // `-parse_seqids`.
    //
    // Where NCBI's order of the ranges is the memory addresses', LOSAT's order is written, so that
    // `SetMaskData`'s checks of all the algorithms and the rest of the build run as in NCBI; the
    // first such sequence is kept (`sequence`: its first Seq-id) for `end_build`. With `-gi_mask`
    // the ranges of a sequence without a GI are written nowhere, and their order does not matter.
    fn add_masks_for_seq_id(
        &mut self,
        ids: &[ParsedSeqId],
        sequence: &str,
    ) -> Result<(), WriteDbError> {
        let Some(mask_data) = self.mask_data.as_mut() else {
            return Ok(());
        };
        let found = mask_data.get_ranges(ids);
        if ranges_are_empty(found.ranges) {
            return Ok(());
        }
        let gi_masks = self.output_db.uses_gi_masks();
        let gis: Vec<i64> = match gi_masks {
            true => ids.iter().filter_map(ParsedSeqId::gi).collect(),
            false => Vec::new(),
        };
        self.output_db.set_mask_data(found.ranges, &gis)?;
        self.found_matching_masks = true;
        let written = !gi_masks || !gis.is_empty();
        if let (None, Some(seq_ids), true) = (&self.address_order, &found.address_ordered, written)
        {
            self.address_order = Some(address_order_refusal(seq_ids, sequence));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1403-1420
    /// `AddFasta`: adds the sequences of one FASTA input. An error ends the
    /// build and removes its files.
    pub fn add_fasta<R: Read>(
        &mut self,
        fasta_file: R,
        warn: Warn<'_>,
    ) -> Result<(), WriteDbError> {
        let mut source =
            FastaBioseqSource::with_parse_ids(fasta_file, self.is_protein, self.parse_ids);
        let result = match self.add_sequences(&mut source, &mut *warn) {
            Ok(true) => Ok(()),
            Ok(false) => Err(WriteDbError::Creation("No sequences added".to_string())),
            Err(error) => Err(error),
        };
        if result.is_err() {
            // The error of the build is the one already in hand.
            let _ = self.end_build(true, warn);
        }
        result
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1278-1286
    /// `SetLeafTaxIds` of `taxids` with the second argument true: the leaf taxids of a copied
    /// database's definition lines (`raw_source::scan_leaf_tax_ids`), which the copy adds to a
    /// line's own.
    pub(crate) fn set_leaf_tax_ids(&mut self, taxids: LeafTaxIds) {
        let _ = writeln!(self.log, "Keep Leaf Taxids: T");
        self.leaf_tax_ids = taxids;
    }

    /// The run's flag of the UTF8String warning, which posts once in a run.
    pub(crate) fn utf8_once(&mut self) -> &mut Utf8TagOnce {
        &mut self.utf8_once
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:428-442
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1478-1512
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1562-1575
    // makeblastdb gives no membership bits and keeps a line's own (`SetMembBits`), and removes
    // no GI: those edits change nothing. A line some of whose keys
    // (`GetDeflineKeys`, `GetSeqIdKey` of each Seq-id) are in the map of leaf taxids gets the
    // ascending union of its own leaf taxids and the map's; any other keeps its own as stored.
    /// `x_EditHeaders` of a copied sequence: `FixTaxId`, then the leaf taxids.
    fn edit_raw_headers(&mut self, headers: &mut DeflineSet) {
        self.taxids.fix_tax_id(headers);
        if self.leaf_tax_ids.is_empty() {
            return;
        }
        for defline in &mut headers.deflines {
            let mut found: Option<BTreeSet<i32>> = None;
            for id in &defline.ids {
                let key = seq_id_key(id);
                if key.is_empty() {
                    continue;
                }
                if let Some(leafs) = self.leaf_tax_ids.get(&key) {
                    found.get_or_insert_with(BTreeSet::new).extend(leafs);
                }
            }
            if let Some(mut taxids) = found {
                taxids.extend(defline.links.iter().flatten().copied());
                defline.links = Some(taxids.into_iter().collect());
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:904-1048
    // The columns of the input other than its mask column are made in the output (`FindColumn`,
    // else `CreateUserColumn`) with all their metadata before the first sequence. Each sequence
    // is checked, its definition lines edited, and it is added in the form of the sequence file
    // with its PIG, its definition lines, the blobs of its columns that are not empty and its
    // masks: the input's ranges, then, for each definition line whose Seq-ids `-mask_data` has
    // ranges of, all the algorithms of that lookup (`GetRanges`), which makes the masks
    // matched. The masks are set when an algorithm has a range.
    //
    // Unlike `add_fasta`, an error does not end the build: makeblastdb leaves it open, and its
    // writer is closed when it is destroyed (`close_at_destruction`).
    /// `AddSequences` of an `IRawSequenceSource`: copies the sequences of `source`, then logs their
    /// count when there were any. `warn` gets the writer's diagnostics as they are posted.
    pub(crate) fn add_raw_sequences(
        &mut self,
        source: &mut RawSeqDbSource<'_>,
        warn: Warn<'_>,
    ) -> Result<(), RawCopyError> {
        let start = std::time::Instant::now();
        let mut output_columns = Vec::with_capacity(source.column_titles().len());
        for index in 0..source.column_titles().len() {
            let title = source.column_titles()[index].clone();
            if title == MASK_DATA_TITLE {
                output_columns.push(None);
                continue;
            }
            let output_id = match self.output_db.find_column(&title) {
                Some(id) => id,
                None => self.output_db.create_user_column(&title)?,
            };
            for (key, value) in source.column_meta_data(index)? {
                self.output_db
                    .add_column_meta_data(output_id, &key, &value)?;
            }
            output_columns.push(Some(output_id));
        }

        let mut count = 0u64;
        while let Some(raw) = source.next()? {
            if raw.sequence.is_empty() {
                return Err(raw_data_error(970, "Error in raw data: no sequence").into());
            }
            if !raw.ambiguities.is_empty() && self.is_protein {
                return Err(raw_data_error(
                    975,
                    "Error in raw data: protein db cannot with ambiguities",
                )
                .into());
            }
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:979-983
            // NCBI's error for a sequence without header is its empty header, which LOSAT's
            // reader gives as a set without lines and the scan of the leaf taxids refused.
            let mut headers = raw.deflines;
            self.edit_raw_headers(&mut headers);
            self.output_db
                .add_raw_sequence(raw.sequence, raw.ambiguities)?;
            self.output_db.set_pig(raw_pig(&headers)?);

            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1010-1031
            // With `-gi_mask` the masks of a copied sequence go under the GIs of all its
            // definition lines, in their order; a sequence without a GI writes none, and the
            // order of its ranges does not matter.
            let gis: Vec<i64> = match self.output_db.uses_gi_masks() {
                true => headers
                    .deflines
                    .iter()
                    .flat_map(|defline| defline.ids.iter().filter_map(ParsedSeqId::gi))
                    .collect(),
                false => Vec::new(),
            };
            let written = !self.output_db.uses_gi_masks() || !gis.is_empty();
            let mut masks = raw.masks;
            if let Some(mask_data) = self.mask_data.as_mut() {
                for defline in &headers.deflines {
                    let found = mask_data.get_ranges(&defline.ids);
                    if ranges_are_empty(found.ranges) {
                        continue;
                    }
                    masks.extend(found.ranges.iter().cloned());
                    self.found_matching_masks = true;
                    if let (None, Some(seq_ids), true) =
                        (&self.address_order, &found.address_ordered, written)
                    {
                        let mut sequence = Vec::new();
                        if let Some(id) = headers.deflines[0].ids.first() {
                            crate::seqid::label::fasta(id, &mut sequence);
                        }
                        self.address_order = Some(address_order_refusal(
                            seq_ids,
                            &String::from_utf8_lossy(&sequence),
                        ));
                    }
                }
            }
            self.output_db.set_deflines(headers);
            for (blob, output_id) in raw.column_blobs.iter().zip(&output_columns) {
                let Some(output_id) = *output_id else {
                    continue;
                };
                if !blob.is_empty() {
                    self.output_db.set_blob_data(output_id, blob)?;
                }
            }
            if !ranges_are_empty(&masks) {
                self.output_db.set_mask_data(&masks, &gis)?;
            }
            self.post_diagnostics(warn);
            count += 1;
        }
        if count != 0 {
            let _ = writeln!(
                self.log,
                "Adding sequences from raw db source; added {count} sequences in {} seconds.",
                format_seconds(start.elapsed().as_secs_f64())
            );
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:100-108
    // An error that ends makeblastdb outside `AddFasta` (which ends its build and removes its
    // files) leaves the build open: the writer is closed when makeblastdb's application is
    // destroyed, after the error's message and `~CBuildDatabase`'s, and writes what was added,
    // the sequence added last included. A `CWriteDBException` of that close is posted; NCBI's
    // destructor lets any other exception end the process, where LOSAT posts the error's
    // message.
    /// The close of the writer by its destructor after an error: the messages it posts (the
    /// diagnostics of the sequence written last, the error of the close), each without its final
    /// newline.
    pub(crate) fn close_at_destruction(&mut self) -> Vec<Vec<u8>> {
        let mut messages = Vec::new();
        let error = self.output_db.close().err();
        self.assert_no_stdout();
        let diagnostics = self.output_db.take_diagnostics();
        if !diagnostics.is_empty() {
            messages.push(
                diagnostics
                    .strip_suffix(b"\n")
                    .unwrap_or(&diagnostics)
                    .to_vec(),
            );
        }
        if let Some(error) = error {
            let _ = self.output_db.destroy_after_failed_close();
            messages.push(match error {
                WriteDbError::Creation(_)
                | WriteDbError::CreationBytes(_)
                | WriteDbError::GiOnlyBetweenKeys { .. } => {
                    [b"Error: [makeblastdb] ".as_slice(), &error.message()].concat()
                }
                other => other.message(),
            });
        }
        messages
    }

    /// `s_FixBioseqDeltas` on a delta that `Seq-inst.length` does not describe
    /// (`Bioseq::delta_fix`), after the definition lines (`x_EditAndAddBioseq`), whose
    /// diagnostics come first where it throws: the length of the raw sequence it made, or its
    /// throw for a protein delta.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:297-304,370,395-397
    #[cold]
    fn delta_fix(
        &mut self,
        fix: DeltaFix,
        bioseq: &mut Bioseq,
        warn: Warn<'_>,
    ) -> Result<(), WriteDbError> {
        let error = match fix {
            DeltaFix::Length(length) => {
                bioseq.length = length;
                return Ok(());
            }
            DeltaFix::ProteinDelta => {
                let best_id = match &bioseq.parsed {
                    Some(parsed) => {
                        let mut out = Vec::new();
                        if let Some(id) = find_best(&parsed.ids, best_rank_score) {
                            crate::seqid::label::fasta(id, &mut out);
                        }
                        out
                    }
                    None => format!("lcl|{}", bioseq.local_id).into_bytes(),
                };
                WriteDbError::Exception(NcbiException::new(
                    PROTEIN_DELTA,
                    [
                        best_id.as_slice(),
                        b": Protein delta sequences are not supported.",
                    ]
                    .concat(),
                ))
            }
        };
        self.post_diagnostics(warn);
        Err(error)
    }

    /// Hands the diagnostics that the writer posted since the last call to
    /// `warn`, the sink of the standard error stream.
    fn post_diagnostics(&mut self, warn: Warn<'_>) {
        let diagnostics = self.output_db.take_diagnostics();
        if !diagnostics.is_empty() {
            warn(diagnostics.strip_suffix(b"\n").unwrap_or(&diagnostics));
        }
        self.assert_no_stdout();
    }

    /// The writer writes to the standard output only for NCBI4na data of another size than its
    /// length needs (`WriteDb::take_stdout`), which makeblastdb's sources never give: the FASTA
    /// source packs exactly that size, and the sequences of `-input_type blastdb` are copied, not
    /// converted.
    fn assert_no_stdout(&mut self) {
        let stdout = self.output_db.take_stdout();
        debug_assert!(
            stdout.is_empty(),
            "makeblastdb's writer wrote {:?} to the standard output",
            String::from_utf8_lossy(&stdout)
        );
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1423-1475
    // `Close` publishes the sequence that was added last, and the
    // serializer posts the diagnostics of its title then: before the lines
    // that `x_EndBuild` writes to the log. `warn` gets them at that point, so
    // that the log and the messages keep makeblastdb's order when both
    // streams have one destination.
    //
    // A `Close` that threw leaves files open: they are removed open (their later writes go to
    // the removed files), and the writer's destructors run after the removal
    // (`WriteDb::destroy_after_failed_close`).
    //
    // LOSAT's refusal of masks whose order NCBI takes from memory addresses (`address_order`)
    // comes here, when the build would otherwise succeed: after a close that did not fail, in an
    // `EndBuild` without `erase`. It ends the build as an error of NCBI's build ends it
    // (`EndBuild` with `erase` true in `AddFasta`, build_db.cpp:1415-1417): the files that
    // `ListFiles` names are removed, which leaves what NCBI leaves after its error at the last
    // sequence (nothing in version 4; the `.nos`, `.not`, `.ntf` and `.nto` files in version 5).
    // The other refusal raised inside the build, `mask_data`'s column that is created too late,
    // comes from `SetMaskData` inside `AddFasta` and so ends the build in the same way. A build
    // that ends with an error (`erase`, a failed close) reports that error and drops the refusal:
    // NCBI's order never reaches a file there.
    /// `EndBuild`: closes the database; with `erase`, or when closing fails,
    /// removes its files. `Ok(false)` when no volume was created. `warn`
    /// takes the diagnostics of the close. A pending refusal of LOSAT (`address_order`) is
    /// raised after a close without `erase` that did not fail, with the files removed.
    pub fn end_build(&mut self, erase: bool, warn: Warn<'_>) -> Result<bool, WriteDbError> {
        let close_error = self.output_db.close().err();
        self.post_diagnostics(warn);
        let refusal = self
            .address_order
            .take()
            .filter(|_| !erase && close_error.is_none());
        let erase = erase || close_error.is_some() || refusal.is_some();
        let vols = self.output_db.list_volumes();
        let files = self.output_db.list_files();
        let _ = writeln!(self.log);
        let mut success = true;
        if vols.is_empty() {
            let _ = writeln!(self.log, "No volumes were created.");
            success = false;
        } else if erase {
            // `Remove` of a `CFile` is `CDirEntry::Remove` with `eRecursive`.
            for file in &files {
                let _ = dir_entry::remove(Path::new(file));
            }
        }
        let _ = writeln!(self.log);
        if let Some(refusal) = refusal {
            return Err(refusal);
        }
        if close_error.is_some() {
            let _ = self.output_db.destroy_after_failed_close();
        }
        // `NCBI_RETHROW` of the close exception as a `CWriteDBException` of `eArgErr` with the
        // message of the close exception: an exception of the close becomes a creation error with
        // its message only (no report).
        match close_error {
            Some(
                WriteDbError::Creation(message)
                | WriteDbError::Options(message)
                | WriteDbError::Database(message)
                | WriteDbError::Other(message),
            ) => Err(WriteDbError::Creation(message)),
            Some(WriteDbError::Exception(e)) => Err(WriteDbError::creation(e.get_msg().to_vec())),
            // A creation error whose text is not UTF-8, an error of a write of LMDB's
            // data file (`Io`) and a failed allocation keep their kind.
            Some(error) => Err(error),
            None => Ok(success),
        }
    }
}

#[cfg(test)]
#[path = "build_db_tests.rs"]
mod tests;
