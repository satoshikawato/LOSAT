//! A database build (`CBuildDatabase`,
//! `c++/src/objtools/blast/seqdb_writer/build_db.cpp`): the log of
//! makeblastdb, the removal of an existing database of the same name, the
//! loop that hands the sequences of a FASTA input to the writer, and the end
//! of the build, which removes the files of a failed one.

use std::io::{Read, Write};
use std::path::Path;

use super::fasta_source::{FastaBioseqSource, Warn};
use super::mask_data::{ranges_are_empty, MaskedRangeSet};
use super::mask_info::MaskSeqId;
use super::time::LocalTime;
use super::writer::WriteDb;
use super::WriteDbError;
use crate::seqdb::{SeqDbError, SeqType};

/// `CBuildDatabase` for a FASTA input without parsed Seq-ids.
pub struct BuildDatabase<'a> {
    is_protein: bool,
    log: &'a mut dyn Write,
    output_db: WriteDb,
    output_dbname: String,
    /// `m_MaskData`: the masked ranges of the sequences, when makeblastdb
    /// read mask data (`SetMaskDataSource`).
    mask_data: Option<MaskedRangeSet>,
    /// `m_FoundMatchingMasks`.
    found_matching_masks: bool,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:852-927
// ```c++
//     while ( !current.empty()  ||  !tail.empty() ) {
//         std::list<string> pretail;
//         if ( !current.empty() ) {
//             NStr::Split(current, kSep, pretail);
//             current.erase();
//             if (pretail.front().empty()
//     ...
//                 ) {
//                 head.clear();
//     ...
//         if ( !head.empty() ) { // empty heads should accept anything
//             string& last = head.back();
//             if (last == DIR_CURRENT) {
//                 if (!next.empty()) {
//                     head.pop_back();
//                 }
//     ...
//             } else if (next == DIR_CURRENT) {
//                 continue;
//             } else if (next.empty()) {
//                 continue; // leave out empty components in most cases
//             } else if (next == DIR_PARENT) {
//     ...
//                 if (last.empty()) {
//                     continue;
//     ...
//                 } else if (last != DIR_PARENT) {
//                     head.pop_back();
//                     continue;
//                 }
//             }
//         }
// ```
/// `CDirEntry::NormalizePath` of an absolute Unix path, without following
/// links: `.` and empty components go, `..` takes its parent with it.
#[cfg(not(windows))]
fn normalize_absolute_path(path: &str) -> String {
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
    if head.len() == 1 && head[0].is_empty() {
        return "/".to_string();
    }
    head.join("/")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:665-670
// ```c++
// string CDirEntry::CreateAbsolutePath(const string& path, ERelativeToWhat rtw)
// {
//     if ( IsAbsolutePath(path) ) {
//         return NormalizePath(path);
//     }
//     string result;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:693-696
// ```c++
//     switch (rtw) {
//         case eRelativeToCwd:
//             result = ConcatPath(CDir::GetCwd(), path);
//             break;
// ```
// The reference platform is Unix. On Windows the path is made absolute by
// the system's own rule; under WASI, where a program sees only the
// directories it was given, the name is used as it is.
/// `CDirEntry::CreateAbsolutePath(dbname)`.
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
    #[cfg(target_os = "wasi")]
    {
        dir_entry_path(dbname).to_string()
    }
    #[cfg(not(any(windows, target_os = "wasi")))]
    {
        if dbname.starts_with('/') {
            return normalize_absolute_path(dbname);
        }
        match std::env::current_dir() {
            Ok(cwd) => {
                let cwd = cwd.to_string_lossy();
                normalize_absolute_path(&format!("{}/{dbname}", cwd.trim_end_matches('/')))
            }
            Err(_) => dbname.to_string(),
        }
    }
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
// ```c++
// void CDirEntry::Reset(const string& path)
// {
//     m_Path = path;
//     size_t len = path.length();
//     // Root dir
//     if ((len == 1)  &&  IsPathSeparator(path[0])) {
//         return;
//     }
//     // Disk name
// #  if defined(DISK_SEPARATOR)
//     if ( (len == 2 || len == 3) && (path[1] == DISK_SEPARATOR) ) {
//         return;
//     }
// #  endif
//     m_Path = DeleteTrailingPathSeparator(path);
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:465-472
// ```c++
// string CDirEntry::DeleteTrailingPathSeparator(const string& path)
// {
//     size_t pos = path.find_last_not_of(DIR_SEPARATORS);
//     if (pos + 1 < path.length()) {
//         return path.substr(0, pos + 1);
//     }
//     return path;
// }
// ```
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
// ```c++
// inline bool s_DirCreate(const string&path, CDir::TCreateFlags flags, mode_t mode)
// {
//     errno = 0;
// #if defined(NCBI_OS_MSWIN)
//     int res = NcbiSys_mkdir(_T_XCSTRING(path));
// #elif defined(NCBI_OS_UNIX)
//     int res = NcbiSys_mkdir(_T_XCSTRING(path), mode);
// #endif
//     if (res != 0) {
//         if (errno != EEXIST) {
//             LOG_ERROR_ERRNO(52, "s_DirCreate(): Cannot create directory: " + path);
//             return false;
//         }
//         // Entry with such name already exists, check its type
//         CDirEntry::EType type = CDirEntry(path).GetType();
//         if (type != CDirEntry::eDir) {
//             LOG_ERROR_NCBI(53, "s_DirCreate(): Path already exist and is not a directory: " + path, CNcbiError::eNotADirectory);
//             return false;
//         }
//         if (F_ISSET(flags, CDir::fCreate_ErrorIfExists)) {
//             LOG_ERROR_NCBI(54, "s_DirCreate(): Directory already exist: " + path, CNcbiError::eFileExists);
//             return false;
//         }
//         if (!F_ISSET(flags, CDir::fCreate_UpdateIfExists)) {
//             return true;
//         }
// ...
//     if ( F_ISSET(flags, CDir::fCreate_PermByUmask)  ||
//           (!F_ISSET(flags, CDir::fCreate_PermByUmask)  &&  !F_ISSET(flags, CDir::fCreate_PermAsParent)  &&
//            NCBI_PARAM_TYPE(NCBI, FileAPIHonorUmask)::GetDefault()) ) {
//         // nothing to do if (umask) or (default mode with "honor umask" global flag)
//         return true;
//     }
//     // Change directory permissions
//     if (NcbiSys_chmod(_T_XCSTRING(path), mode) != 0) {
//         LOG_ERROR_ERRNO(55, "CDir::Create(): Cannot set mode for directory: " + path);
//         return false;
//     }
//     return true;
// ```
// The mode is the default of a `CDir`, rwxr-xr-x. `[NCBI] FileAPIHonorUmask`
// is false unless the registry sets it (ncbifile.cpp:161; LOSAT refuses a
// registry that does, `crate::blastinput::ncbi_environment`), so the new
// directory gets exactly that mode, whatever the umask.
/// `s_DirCreate(path, fCreate_Default, mode)`: makes one directory; `true`
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
// ```c++
// bool CDir::CreatePath(TCreateFlags flags) const
// {
//     if (GetPath().empty()) {
//         LOG_ERROR_NCBI(59, "CDir::CreatePath(): Path is empty", CNcbiError::eInvalidArgument);
//         return false;
//     }
//     string path(CreateAbsolutePath(GetPath()));
//     if (path.empty()) {
//         LOG_ERROR_NCBI(60, "CDir::CreatePath(): Cannot create absolute path from: " + path, CNcbiError::eInvalidArgument);
//         return false;
//     }
//     if (path[path.length()-1] == GetPathSeparator()
// #if defined(NCBI_OS_MSWIN)
//         &&  path.length() != 3
//         // Special case -- for path like "C:\" dont remove a last separator, it represent a root directory
// #endif
//         ) {
//         path.erase(path.length() - 1);
//     }
//
//     // Find all missed parts of the path
//
//     CTempString tmp(path); // existent part of a path
//     std::list<CTempString> missed_parts;
//
//     while (!tmp.empty()  &&  !CDirEntry(tmp).Exists()) {
//         size_t pos = tmp.find_last_of(DIR_SEPARATORS);
//         if (pos == NPOS) {
//             break;
//         }
//         CTempString part(tmp.substr(pos+1));
//         missed_parts.push_front(part);
//         tmp.erase(pos);
//     }
// ...
//     // Path exists?
//     if (missed_parts.empty()) {
//         // check existence and behave depends on flags
//         if (!s_DirCreate(path, flags, mode)) {
//             LOG_ERROR(96, "CDir::CreatePath(): Cannot create path: " + GetPath());
//             return false;
//         }
//         return true;
//     }
//
//     // Create missed subdirectories
//     string p = tmp;
//     for (auto i : missed_parts) {
//         p += GetPathSeparator();
//         p += i;
//         if (!s_DirCreate(p, flags, mode)) {
//             LOG_ERROR(97, "CDir::CreatePath(): Cannot create path: " + GetPath());
//             return false;
//         }
//     }
//     return true;
// }
// ```
// `CDirEntry::Exists` is true for an entry of any type, a link that leads
// nowhere included (`GetType()` without following links).
/// `CDir(dir).CreatePath()`: makes the directory `dir` and its missing
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
    // NCBI's path is absolute here. LOSAT's is relative under WASI
    // (`create_absolute_path`): its first component has no separator
    // before it and can be missing too.
    let mut relative_start = false;
    while !tmp.is_empty() && std::fs::symlink_metadata(tmp).is_err() {
        let Some(pos) = tmp.rfind(is_separator) else {
            missed_parts.push(tmp);
            tmp = "";
            relative_start = true;
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
        if !std::mem::take(&mut relative_start) {
            p.push(std::path::MAIN_SEPARATOR);
        }
        p.push_str(part);
        if !dir_create(&p) {
            return false;
        }
    }
    true
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1051-1071
// ```c++
// void CBuildDatabase::CreateDirectories(const string& dbname)
// {
//     CDirEntry dir_entry(dbname);
//     string dir_name = dir_entry.GetDir(CDirEntry::eIfEmptyPath_Empty);
//     if (dir_name.empty()) {
//         return;
//     }
//
//     CDir d(dir_name);
//     if ( !d.Exists() ) {
//         if ( !d.CreatePath() ) {
//             string msg("Failed to create directory '" + d.GetName() + "'");
//             NCBI_THROW(CMultisourceException, eOutputFileError, msg);
//         }
//     }
//     if (!d.CheckAccess(CDirEntry::fWrite)) {
//         string msg("You do not have write permissions on '" +
//                    d.GetName() + "'");
//         NCBI_THROW(CMultisourceException, eOutputFileError, msg);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-368
// ```c++
// void CDirEntry::SplitPath(const string& path, string* dir,
//                           string* base, string* ext)
// {
//     // Get file name
//     size_t pos = path.find_last_of(ALL_SEPARATORS);
//     string filename = (pos == NPOS) ? path : path.substr(pos+1);
//
//     // Get dir
//     if ( dir ) {
//         *dir = (pos == NPOS) ? kEmptyStr : path.substr(0, pos+1);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:1735-1748
// ```c++
//     const char* path = GetPath().c_str();
//     int         mode = F_OK;
//
//     if ( access_mode & fRead)    mode |= R_OK;
//     if ( access_mode & fWrite)   mode |= W_OK;
//     if ( access_mode & fExecute) mode |= X_OK;
//
//     // Use euidaccess() where possible
// #  if defined(HAVE_EUIDACCESS)
//     if (euidaccess(path, mode) != 0) {
//         CNcbiError::SetFromErrno(path);
//         return false;
//     }
//     return true;
// ```
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
// makeblastdb stops with the write-permission message. LOSAT reproduces that
// part of `CheckAccess`: a path that does not resolve.
//
// Not reproduced: the test of the write permission itself (`euidaccess`
// with `W_OK`). The standard library has no call that asks for it, and a
// trial file would touch the directory. A directory that exists and cannot
// be written to therefore does not stop LOSAT here, before the log, with
// "You do not have write permissions on ..."; LOSAT goes on, writes the log
// header, removes what it can of an existing database (nothing, in such a
// directory), and fails when it creates the first file of the first volume,
// with the system's message for that file and the same exit status 255.
/// `CBuildDatabase::CreateDirectories`.
fn create_directories(dbname: &str) -> Result<(), WriteDbError> {
    let dir_entry = dir_entry_path(dbname);
    let Some(at) = dir_entry.rfind(is_any_separator) else {
        return Ok(());
    };
    let d = dir_entry_path(&dir_entry[..=at]);
    let name = d.rfind(is_any_separator).map_or(d, |at| &d[at + 1..]);
    if !Path::new(d).is_dir() && !create_path(d) {
        return Err(WriteDbError::Other(format!(
            "Failed to create directory '{name}'"
        )));
    }
    if std::fs::metadata(d).is_err() {
        return Err(WriteDbError::Other(format!(
            "You do not have write permissions on '{name}'"
        )));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2629-2641
// ```c++
//     extn.push_back(kExtnMol + "al");   // alias file
//     extn.push_back(kExtnMol + "in");   // index file
//     extn.push_back(kExtnMol + "hr");   // header file
//     extn.push_back(kExtnMol + "sq");   // sequence file
//     extn.push_back(kExtnMol + "ni");   // ISAM numeric index file
//     extn.push_back(kExtnMol + "nd");   // ISAM numeric data file
//     if (dbver == eBDB_Version4) {
//     	extn.push_back(kExtnMol + "si");   // ISAM string index file
//     	extn.push_back(kExtnMol + "sd");   // ISAM string data file
//     }
//     extn.push_back(kExtnMol + "pi");   // ISAM PIG index file
//     extn.push_back(kExtnMol + "pd");   // ISAM PIG data file
//     extn.push_back(kExtnMol + "js");   // BLASTDB metadata
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2647-2656
// ```c++
//     // Contain masking information
//     extn.push_back(kExtnMol + "aa");   // ISAM mask index file
//     extn.push_back(kExtnMol + "ab");   // ISAM mask data file (big-endian)
//     extn.push_back(kExtnMol + "ac");   // ISAM mask data file (little-endian)
//     extn.push_back(kExtnMol + "og");   // OID to GI file
//     extn.push_back(kExtnMol + "hi");   // ISAM sequence hash index file
//     extn.push_back(kExtnMol + "hd");   // ISAM sequence hash data file
//     extn.push_back(kExtnMol + "ti");   // ISAM trace id index file
//     extn.push_back(kExtnMol + "td");   // ISAM trace id data file
//     extn.push_back(kExtnMol + "xm");   // OID mask file to exclude models
// ```
/// `SeqDB_GetFileExtensions(is_protein, extn, eBDB_Version4)` without the
/// letter of the molecule.
const VOLUME_EXTENSIONS: [&str; 20] = [
    "al", "in", "hr", "sq", "ni", "nd", "si", "sd", "pi", "pd", "js", "aa", "ab", "ac", "og", "hi",
    "hd", "ti", "td", "xm",
];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2663-2668
// ```c++
// 	static const char * ext[]={"db", "os", "ot", "tf", "to", "db-lock", "tf-lock", NULL};
// 	extn.clear();
// 	const string kExtnMol(1, db_is_protein ? 'p' : 'n');
// 	for(const char ** p=ext; *p != NULL; p++) {
// 		extn.push_back(kExtnMol + (*p));
// 	}
// ```
/// `SeqDB_GetLMDBFileExtensions` without the letter of the molecule.
const LMDB_EXTENSIONS: [&str; 7] = ["db", "os", "ot", "tf", "to", "db-lock", "tf-lock"];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1542-1593
// ```c++
// bool DeleteBlastDb(const string& dbpath, CSeqDB::ESeqType seq_type)
// {
//     int num_files_removed = 0;
//     vector<string> db_files, alias_files;
//     bool is_protein = (seq_type == CSeqDB::eProtein);
//
//     vector<string> extn;
//     SeqDB_GetFileExtensions( is_protein, extn, eBDB_Version4);
//     vector<string> lmdb_extn;
//     SeqDB_GetLMDBFileExtensions(is_protein, lmdb_extn);
//     ITERATE(vector<string>, lmdb, lmdb_extn) {
//     	CNcbiOstrstream oss;
//     	oss << dbpath << "." << *lmdb;
//         const string fname = CNcbiOstrstreamToString(oss);
//         if (CFile(fname).Remove()) {
//         	LOG_POST(Info << "Deleted " << fname);
//             num_files_removed++;
//         }
//         else {
//         	unsigned int index = 0;
//         	string vfname = dbpath + "." + NStr::IntToString(index/10) +
//         			        NStr::IntToString(index%10) + "." + *lmdb;
//         	while (CFile(vfname).Remove()) {
//         		index++;
//         		vfname = dbpath + "." + NStr::IntToString(index/10) +
//         	    	     NStr::IntToString(index%10) + "." + *lmdb;
//
//         	}
//         }
//     }
//
//     try { CSeqDB::FindVolumePaths(dbpath, seq_type, db_files, &alias_files); }
//     catch (...) {}    // ignore any errors from the invocation above
//     ITERATE(vector<string>, f, db_files) {
//         ITERATE(vector<string>, e, extn) {
//             CNcbiOstrstream oss;
//             oss << *f << "." << *e;
//             const string fname = CNcbiOstrstreamToString(oss);
//             if (CFile(fname).Remove()) {
//                 LOG_POST(Info << "Deleted " << fname);
//                 num_files_removed++;
//             }
//         }
//     }
//     ITERATE(vector<string>, f, alias_files) {
//         if (CFile(*f).Remove()) {
//             LOG_POST(Info << "Deleted " << *f);
//             num_files_removed++;
//         }
//     }
//     return static_cast<bool>(num_files_removed != 0);
// }
// ```
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
//   directory `my data`) the call above throws, nothing is removed, the new
//   files are written beside the old ones (an old alias file then hides the
//   new database), and makeblastdb ends with exit status 2 when it opens its
//   result. LOSAT removes the existing database of that path and the build
//   ends with status 0.
// * The volumes that an alias file names are those in the alias file's own
//   directory. NCBI takes a volume of the same name in the current directory
//   first, and so removes files that do not belong to the database being
//   replaced.
/// `DeleteBlastDb`: removes the files of the database `dbpath`; `true` when
/// a file was removed.
fn delete_blast_db(dbpath: &str, is_protein: bool) -> Result<bool, WriteDbError> {
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

    let mut num_files_removed = 0;
    let mut remove = |name: &str| {
        let removed = std::fs::remove_file(name).is_ok();
        num_files_removed += u32::from(removed);
        removed
    };
    for lmdb in LMDB_EXTENSIONS {
        if !remove(&format!("{dbpath}.{letter}{lmdb}")) {
            let mut index = 0u32;
            // The removals of the volume files are not counted by NCBI.
            while std::fs::remove_file(format!(
                "{dbpath}.{}{}.{letter}{lmdb}",
                index / 10,
                index % 10
            ))
            .is_ok()
            {
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

/// `ostream << double` with the default format (six significant digits,
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
    // ```c++
    //     CreateDirectories(dbname);
    //     const string output_dbname = CDirEntry::CreateAbsolutePath(dbname);
    //     m_LogFile << "\n\nBuilding a new DB, current time: "
    //               << CTime(CTime::eCurrent).AsString() << endl;
    //
    //     m_LogFile << "New DB name:   " << output_dbname << endl;
    //     m_LogFile << "New DB title:  " << title << endl;
    //     const string mol_type(is_protein ? "Protein" : "Nucleotide");
    //     m_LogFile << "Sequence type: " << mol_type << endl;
    //     if (DeleteBlastDb(output_dbname, ParseMoleculeTypeString(mol_type))) {
    //         m_LogFile << "Deleted existing " << mol_type
    //             << " BLAST database named " << output_dbname << endl;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1270-1273
    // ```c++
    // void CBuildDatabase::SetMembBits(const TLinkoutMap & membbits,
    //                                  bool                keep_mbits)
    // {
    //     m_LogFile << "Keep MBits: " << (keep_mbits ? "T" : "F") << endl;
    // ```
    /// The constructor, then `SetMembBits(no_bits, true)` as makeblastdb
    /// calls it. `now` is the time of the build.
    pub fn new(
        dbname: &str,
        title: &str,
        is_protein: bool,
        log: &'a mut dyn Write,
        now: &LocalTime,
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
        let _ = writeln!(log, "New DB title:  {title}");
        let _ = writeln!(log, "Sequence type: {mol_type}");
        if delete_blast_db(&output_dbname, is_protein)? {
            let _ = writeln!(
                log,
                "Deleted existing {mol_type} BLAST database named {output_dbname}"
            );
        }
        let _ = writeln!(log, "Keep MBits: T");
        // The limit of 1000*1000*1000 that the constructor sets
        // (build_db.cpp:1200) is replaced by makeblastdb before a sequence
        // is added.
        let output_db = WriteDb::new(
            &output_dbname,
            is_protein,
            title.as_bytes(),
            now.index_date().as_bytes(),
            1000 * 1000 * 1000,
        );
        let mut output_db = output_db;
        output_db.set_column_date_fallback(now.log_time().as_bytes());
        Ok(Self {
            is_protein,
            log,
            output_db,
            output_dbname,
            mask_data: None,
            found_matching_masks: false,
        })
    }

    /// `GetOutputDbName`: the absolute name of the database.
    pub fn output_dbname(&self) -> &str {
        &self.output_dbname
    }

    /// The log stream (`m_LogFile`).
    pub fn log(&mut self) -> &mut dyn Write {
        self.log
    }

    /// `SetMaxFileSize`.
    pub fn set_max_file_size(&mut self, max_file_size: u64) {
        self.output_db.set_max_file_size(max_file_size);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1583-1594
    // ```c++
    // int
    // CBuildDatabase::RegisterMaskingAlgorithm(EBlast_filter_program program,
    //                                          const string        & options,
    //                                          const string        & name)
    // {
    //     return m_OutputDb->RegisterMaskAlgorithm(program, options, name);
    // }
    // ```
    /// `RegisterMaskingAlgorithm(program, options, name)`: see
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
    // ```c++
    // int
    // CBuildDatabase::RegisterMaskingAlgorithm(const string        & program,
    //                                          const string        & description,
    //                                          const string        & options)
    // {
    //     return m_OutputDb->RegisterMaskAlgorithm(program, description, options);
    // }
    // ```
    // The first argument of this overload is the id of `-mask_id` (the call
    // in makeblastdb.cpp:966 passes `id_list[i], desc_list[i], opts`).
    /// `RegisterMaskingAlgorithm(id, description, options)`: see
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
    // ```c++
    // void CBuildDatabase::SetMaskDataSource(IMaskDataSource & ranges)
    // {
    //     m_MaskData.Reset(& ranges);
    // }
    // ```
    /// `SetMaskDataSource`: the masked ranges of the sequences. makeblastdb
    /// calls it when the first mask file with masks has been read, and
    /// keeps adding to the set after that, so the caller hands the finished
    /// set over before the first sequence is added.
    pub fn set_mask_data_source(&mut self, ranges: MaskedRangeSet) {
        self.mask_data = Some(ranges);
    }

    /// The set that `set_mask_data_source` handed over, to add to. makeblastdb
    /// creates the set when it reads the first mask list that is not empty,
    /// hands it over at once and goes on inserting; a caller that follows it
    /// calls `set_mask_data_source(MaskedRangeSet::new())` at that point and
    /// inserts through this.
    pub fn mask_data_source_mut(&mut self) -> Option<&mut MaskedRangeSet> {
        self.mask_data.as_mut()
    }

    /// `m_FoundMatchingMasks`: whether a sequence has been given mask data.
    pub fn found_matching_masks(&self) -> bool {
        self.found_matching_masks
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1204-1210
    // ```c++
    // CBuildDatabase::~CBuildDatabase()
    // {
    //     if (m_MaskData.NotEmpty() && !m_FoundMatchingMasks) {
    //         ERR_POST(Error << "No sequences matched any of the masks provided.\n"
    //                        << "Please ensure that the -parse_seqids option is used "
    //                        << "in the\nfiltering program as well as makeblastdb.");
    //     }
    // ```
    // The message is posted when the `CBuildDatabase` is destroyed: in
    // makeblastdb, when the application object goes (`m_DB` is a member), so
    // after the log of `EndBuild`, after the exception message of a failed
    // build and after the metadata file. The caller of LOSAT's build calls
    // this at that point.
    /// The destructor's diagnostic: `warn` gets the message (without its
    /// final newline) when mask data were read and no sequence matched.
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:859-876
    // ```c++
    //         if(bs->IsAa() != m_IsProtein ){
    //                     bs = src.GetNext();
    //                     continue;
    //        }
    //
    //        if (bs->GetLength() > 0x7fffffff)
    //        {
    //                    NCBI_THROW(CWriteDBException,
    //                    eArgErr,
    //                    "Sequences longer than 2,147,483,647 bases are not supported.  Offending sequence is " + bioseq_id);
    //        }
    //
    //         if ((bs->GetLength() == 0) || (!x_EditAndAddBioseq(bs, NULL, add_pig))){
    //             m_LogFile << "Ignoring sequence '" << bioseq_id
    //                       << "' as it has no sequence data" << endl;
    //             bs = src.GetNext();
    //             continue;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-487
    // ```c++
    //     CRef<CBlast_def_line_set> headers =
    //         CWriteDB::ExtractBioseqDeflines(*bs, m_ParseIDs, m_LongIDs,
    //                 m_ScanBioseq4CFastaReaderUsrObjct);
    //
    //     x_EditHeaders(headers);
    //
    //     // Add the sequence
    //     if (sv) {
    //         m_OutputDb->AddSequence(*bs, *sv);
    //     } else {
    //         bs = s_FixBioseqDeltas(bs);
    //         if(bs->GetInst().CanGetSeq_data())
    //         	m_OutputDb->AddSequence(*bs);
    //         else
    //         	return false;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:894-901
    // ```c++
    //     if (count) {
    //         double t = sw.Elapsed();
    //
    //         m_LogFile << "Adding sequences from FASTA; added "
    //                   << count << " sequences in " << t << " seconds." << endl;
    //     }
    //
    //     return found;
    // ```
    // The FASTA source applies `s_FixBioseqDeltas` before it hands a
    // sequence over; its errors arrive as errors of the source.
    /// `AddSequences(IBioseqSource &)`: `true` when a sequence was added.
    fn add_sequences<R: Read>(
        &mut self,
        source: &mut FastaBioseqSource<R>,
        warn: Warn<'_>,
    ) -> Result<bool, WriteDbError> {
        let mut found = false;
        let start = std::time::Instant::now();
        let mut count = 0u64;
        while let Some(bioseq) = source.next_bioseq(warn).map_err(|error| {
            if error.out_of_memory {
                WriteDbError::OutOfMemory
            } else {
                WriteDbError::Creation(error.message)
            }
        })? {
            let bioseq_id = format!("lcl|{}", bioseq.local_id);
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
                self.output_db.extract_bioseq_deflines(&bioseq);
                if bioseq.seq_data.is_some() {
                    // The masks are looked up by the number of the
                    // sequence; one that is no `local id` of 31 bits
                    // matches none.
                    let mask_id = self
                        .mask_data
                        .is_some()
                        .then(|| bioseq.local_id.parse::<i32>().ok())
                        .flatten();
                    self.output_db.add_sequence(bioseq)?;
                    added = true;
                    self.add_masks_for_seq_id(mask_id)?;
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
    // ```c++
    // void
    // CBuildDatabase::x_AddMasksForSeqId(const list< CRef<CSeq_id> >& ids)
    // {
    //     if (m_MaskData.Empty()) {
    //         return;
    //     }
    //
    //     const CMaskedRangesVector& rng = m_MaskData->GetRanges(ids);
    //     if (rng.empty()) {
    //         return;
    //     }
    //
    //     vector <TGi> gis;
    //     ITERATE(list< CRef<CSeq_id> >, id, ids) {
    //         if ((*id)->IsGi()) {
    //             gis.push_back((*id)->GetGi());
    //         }
    //     }
    //     m_OutputDb->SetMaskData(rng, gis);
    //     m_FoundMatchingMasks = true;
    // }
    // ```
    // The Seq-ids of a sequence without `-parse_seqids` are its `local id`
    // (`mask_id`, `None` when it has no 31-bit one: it matches nothing, but
    // the other Seq-ids of masks are compared the same way).
    fn add_masks_for_seq_id(&mut self, mask_id: Option<i32>) -> Result<(), WriteDbError> {
        let Some(mask_data) = self.mask_data.as_mut() else {
            return Ok(());
        };
        let ids: &[MaskSeqId] = match mask_id {
            Some(n) => &[MaskSeqId::LocalId(n)],
            None => &[],
        };
        let ranges = mask_data.get_ranges(ids)?;
        if ranges_are_empty(ranges) {
            return Ok(());
        }
        self.output_db.set_mask_data(ranges)?;
        self.found_matching_masks = true;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1403-1420
    // ```c++
    //     if (fasta_file) {
    //         CFastaBioseqSource fbs(fasta_file,
    //                                m_IsProtein,
    //                                m_ParseIDs,
    //                                m_LongIDs);
    //
    //         try {
    //             success = AddSequences(fbs);
    // 	    if (success == false)
    //             	NCBI_THROW(CWriteDBException, eFileErr, "No sequences added");
    //
    //         }
    //         catch (...) {
    //             EndBuild(true);
    //             throw;
    //         }
    //     }
    //     return success;
    // ```
    /// `AddFasta`: adds the sequences of one FASTA input. An error ends the
    /// build and removes its files.
    pub fn add_fasta<R: Read>(
        &mut self,
        fasta_file: R,
        warn: Warn<'_>,
    ) -> Result<(), WriteDbError> {
        let mut source = FastaBioseqSource::new(fasta_file, self.is_protein);
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

    /// Hands the diagnostics that the writer posted since the last call to
    /// `warn`, the sink of the standard error stream.
    fn post_diagnostics(&mut self, warn: Warn<'_>) {
        let diagnostics = self.output_db.take_diagnostics();
        if !diagnostics.is_empty() {
            warn(diagnostics.strip_suffix(b"\n").unwrap_or(&diagnostics));
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1423-1475
    // ```c++
    // bool CBuildDatabase::EndBuild(bool erase)
    // {
    //     try {
    //         m_OutputDb->Close();
    //         return x_EndBuild(erase, NULL);
    //     } catch (const CException& e) {
    //         return x_EndBuild(true, erase ? NULL : &e);
    //     ...
    // bool CBuildDatabase::x_EndBuild(bool erase, const CException * close_exception)
    // {
    //     bool success = true;
    //
    //     vector<string> vols;
    //     vector<string> files;
    //
    //     m_OutputDb->ListVolumes(vols);
    //     m_OutputDb->ListFiles(files);
    //
    //     m_LogFile << endl;
    //
    //     _ASSERT(vols.empty() == files.empty());
    //
    //     if (vols.empty()) {
    //         m_LogFile << "No volumes were created."
    //                   << endl;
    //
    //         success = false;
    //     } else {
    //         ITERATE(vector<string>, iterf, files) {
    //             if (erase) {
    //                 CFile(*iterf).Remove();
    //                 _TRACE("Deleting " << *iterf);
    //             }
    //         }
    //     }
    //
    //     m_LogFile << endl;
    //
    //     if (close_exception) {
    //         NCBI_RETHROW(*close_exception, CWriteDBException, eArgErr,
    //                      close_exception->GetMsg());
    //     }
    //
    //     return success;
    // }
    // ```
    // `Close` publishes the sequence that was added last, and the
    // serializer posts the diagnostics of its title then: before the lines
    // that `x_EndBuild` writes to the log. `warn` gets them at that point, so
    // that the log and the messages keep makeblastdb's order when both
    // streams have one destination.
    /// `EndBuild`: closes the database; with `erase`, or when closing fails,
    /// removes its files. `Ok(false)` when no volume was created. `warn`
    /// takes the diagnostics of the close.
    pub fn end_build(&mut self, erase: bool, warn: Warn<'_>) -> Result<bool, WriteDbError> {
        let close_error = self.output_db.close().err();
        self.post_diagnostics(warn);
        let erase = erase || close_error.is_some();
        if erase {
            self.output_db.abandon();
        }
        let vols = self.output_db.list_volumes();
        let files = self.output_db.list_files();
        let _ = writeln!(self.log);
        let mut success = true;
        if vols.is_empty() {
            let _ = writeln!(self.log, "No volumes were created.");
            success = false;
        } else if erase {
            for file in &files {
                let _ = std::fs::remove_file(file);
            }
        }
        let _ = writeln!(self.log);
        match close_error {
            Some(
                WriteDbError::Creation(message)
                | WriteDbError::Options(message)
                | WriteDbError::Other(message),
            ) => Err(WriteDbError::Creation(message)),
            // LOSAT's own errors (a file that could not be written or
            // renamed) and a failed allocation keep their kind.
            Some(error) => Err(error),
            None => Ok(success),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(not(windows))]
    #[test]
    fn absolute_paths_are_normalized_as_ncbi_does() {
        assert_eq!(normalize_absolute_path("/a//b/./c/../d/"), "/a/b/d");
        assert_eq!(normalize_absolute_path("/a/b/.."), "/a");
        assert_eq!(normalize_absolute_path("/../a"), "/a");
        assert_eq!(normalize_absolute_path("/a/./"), "/a");
        assert_eq!(normalize_absolute_path("/"), "/");
        assert_eq!(normalize_absolute_path("/a/.."), "/");
    }

    // Values makeblastdb 2.17.0 printed in "added N sequences in T seconds."
    #[test]
    fn seconds_print_as_six_significant_digits() {
        assert_eq!(format_seconds(0.000254154), "0.000254154");
        assert_eq!(format_seconds(0.00107908), "0.00107908");
        assert_eq!(format_seconds(0.0001719), "0.0001719");
        assert_eq!(format_seconds(0.000169992), "0.000169992");
        assert_eq!(format_seconds(5.5), "5.5");
        assert_eq!(format_seconds(12.3456789), "12.3457");
        assert_eq!(format_seconds(0.00001234567), "1.23457e-05");
        assert_eq!(format_seconds(1234567.0), "1.23457e+06");
        assert_eq!(format_seconds(100.0), "100");
        assert_eq!(format_seconds(0.0), "0");
    }
}
