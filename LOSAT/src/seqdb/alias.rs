//! Alias files (`.pal` / `.nal`) and the list of volumes behind a database name.
//!
//! LOSAT reads the keys a multi-volume database of makeblastdb uses, `TITLE`
//! and `DBLIST`, and the keys `NSEQ` and `LENGTH` of blastdb_aliastool when
//! their values equal the sums over the databases they cover (`check_sums`).
//! The keys that filter sequences or replace the statistics are explicit
//! errors (not supported yet); unknown keys are ignored, as in NCBI.
//!
//! Two ways lead from a name to its volumes. `resolve_in` is NCBI's: the
//! name of a search program's `-db` is looked for in the directories of a
//! [`SearchPath`] (the current directory, `BLASTDB`, the registry value
//! `[BLAST] BLASTDB`), names inside alias files are looked for there too when
//! they are found nowhere else, and volumes are recorded under normalised
//! directory names. `resolve_one_name` takes the path of one database as it
//! is and looks nowhere else (`LOSAT makedb`).

use std::collections::VecDeque;
use std::ffi::OsStr;
use std::path::{Path, MAIN_SEPARATOR};

use super::file::DbFile;
use super::volume::with_extension;
use super::{SeqDbError, SeqType};

/// One node of the alias tree: the database name given by the caller (the
/// top node) or one alias file.
#[derive(Debug, Default)]
pub(super) struct AliasNode {
    /// `TITLE` of the alias file, when it has that key.
    pub title: Option<Vec<u8>>,
    /// Raw value of the `NSEQ` key of the alias file, when it has that key.
    pub nseq: Option<Vec<u8>>,
    /// Raw value of the `LENGTH` key of the alias file, when it has that key.
    pub length: Option<Vec<u8>>,
    /// Alias files named in `DBLIST`, in order.
    pub subnodes: Vec<AliasNode>,
    /// Volumes named in `DBLIST` (paths without extension), in order.
    pub vol_names: Vec<String>,
    /// Path of the alias file of this node (`m_ThisName`); `None` for the
    /// top node, which NCBI names "-".
    pub this_name: Option<String>,
}

/// The resolved database: the alias tree and its volumes in NCBI's order.
#[derive(Debug)]
pub(super) struct AliasTree {
    pub root: AliasNode,
    /// Volume paths without extension, sorted as NCBI opens them.
    pub volumes: Vec<String>,
}

/// Alias keys that LOSAT recognises and does not support: they select
/// sequences (alias_files.txt "Type: filter"), replace the values used for
/// statistics ("meta-data"), or name data LOSAT does not read. `NSEQ` and
/// `LENGTH` are read (`check_sums`) and are not among them. They are the
/// keys NCBI's reader looks up (seqdbalias.cpp), except `MAX_SEQ_LENGTH`,
/// which it never applies (seqdbimpl.cpp:972-976 takes the maximum from the
/// volumes) and which is ignored here too, plus `VDBLIST` of
/// blastdb_aliastool.
const UNSUPPORTED_KEYS: [&str; 14] = [
    "FIRST_OID",
    "LAST_OID",
    "GILIST",
    "TILIST",
    "SEQIDLIST",
    "TAXIDLIST",
    "OIDLIST",
    "MEMB_BIT",
    "MASKLIST",
    "VDBLIST",
    "STATS_NSEQ",
    "STATS_TOTLEN",
    "MIN_SEQ_LENGTH",
    "OID_MASK_TYPE",
];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:160-182
// ```c++
// bool CSeqDBAtlas::GetFileSizeL(const string & fname, TIndx &length)
// {
//     ...
//     pair<bool, TIndx> val;
//     CFile whole(fname);
//     Int8 file_length = whole.GetLength();
//
//     if (file_length >= 0) {
//         val.first = true;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3204-3216
// ```c++
// Int8 CFile::GetLength(void) const
// {
//     TNcbiSys_stat st;
//     if (NcbiSys_stat(_T_XCSTRING(GetPath()), &st) != 0) {
//         LOG_ERROR_ERRNO(40, "CFile:GetLength(): stat() failed for: " + GetPath());
//         return -1L;
//     }
//     if ( GetType(st) != eFile ) {
//         LOG_ERROR_NCBI(40, "CFile:GetLength(): Not a file: " + GetPath(), CNcbiError::eOperationNotPermitted);
//         return -1L;
//     }
//     return st.st_size;
// }
// ```
/// `CSeqDBAtlas::DoesFileExist`: `path` is a regular file, or a symbolic link
/// to one.
pub(super) fn exists(path: &Path) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_file())
}

fn dir_name(path: &str) -> &str {
    path.rfind(MAIN_SEPARATOR).map_or("", |at| &path[..at])
}

fn base_name(path: &str) -> &str {
    path.rfind(MAIN_SEPARATOR)
        .map_or(path, |at| &path[at + 1..])
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:88-107
// ```c++
//     int slen = s.Size();
//
//     if (slen > 4) {
//         string extn(s.GetEnd()-4, s.GetEnd());
//         string extn2(extn, 2, 4);
//         ...
//         if (extn[0] == '.'     &&
//             (extn[1] == 'n' || extn[1] == 'p') &&
//             (extn2 == "al"  || extn2 == "in" || extn2 == "db")) {
//             ...
//             s.Resize(slen - 4);
//         }
//     }
// ```
// `CSeqDB_Path::FindBaseName` is `SeqDB_RemoveExtn(SeqDB_RemoveDirName(..))`
// (seqdbgeneral.hpp:975-979): the file name without a final `.pal`, `.nal`,
// `.pin`, `.nin`, `.pdb` or `.ndb`. `SeqDB_CompareVolume` compares these.
fn path_base_name(path: &str) -> &str {
    remove_extn(base_name(path))
}

/// `SeqDB_RemoveExtn` (quoted above): `s` without a final `.pal`, `.nal`,
/// `.pin`, `.nin`, `.pdb` or `.ndb`, when something precedes it.
fn remove_extn(s: &str) -> &str {
    let bytes = s.as_bytes();
    let len = bytes.len();
    if len > 4
        && bytes[len - 4] == b'.'
        && (bytes[len - 3] == b'n' || bytes[len - 3] == b'p')
        && matches!(&bytes[len - 2..], b"al" | b"in" | b"db")
    {
        &s[..len - 4]
    } else {
        s
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:284-295
// ```c++
// void SeqDB_ConvertOSPath(string & dbs)
// {
//     // See also CDirEntry::ConvertToOSPath()
//
//     char delim = CDirEntry::GetPathSeparator();
//
//     for(size_t i = 0; i<dbs.size(); i++) {
//         if (dbs[i] == '/' || dbs[i] == '\\') {
//             dbs[i] = delim;
//         }
//     }
// }
// ```
fn fix_delimiters(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c == '/' || c == '\\' {
                MAIN_SEPARATOR
            } else {
                c
            }
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:140-181
// ```c++
//     if (two.Empty()) {
//         // We only use the extension if there is a filename.
//         one.GetString(outp);
//         return;
//     }
//
//     bool only_two = false;
//
//     if (one.Empty() || two[0] == delim) {
//         only_two = true;
//     }
//
//     // Drive letter test for CP/M derived systems
//     if (delim == '\\'   &&
//         two.Size() > 3  &&
//         isalpha(two[0]) &&
//         two[1] == ':'   &&
//         two[2] == '\\') {
//
//         only_two = true;
//     }
//     ...
//     one.GetString(outp);
//
//     if (outp[outp.size() - 1] != delim) {
//         outp += delim;
//     }
//
//     outp.append(two.GetBegin(), two.GetEnd());
// ```
fn combine_path(one: &str, two: &str) -> String {
    if two.is_empty() {
        return one.to_string();
    }
    let bytes = two.as_bytes();
    let mut only_two = one.is_empty() || two.starts_with(MAIN_SEPARATOR);
    if MAIN_SEPARATOR == '\\'
        && bytes.len() > 3
        && bytes[0].is_ascii_alphabetic()
        && bytes[1] == b':'
        && bytes[2] == b'\\'
    {
        only_two = true;
    }
    if only_two {
        return two.to_string();
    }
    let mut out = one.to_string();
    if !out.ends_with(MAIN_SEPARATOR) {
        out.push(MAIN_SEPARATOR);
    }
    out.push_str(two);
    out
}

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
// `DIR_SEPARATORS` is "/" on Unix and "/\\" on Windows (ncbifile.cpp:110, 114).
fn delete_trailing_path_separator(path: &str) -> &str {
    if cfg!(windows) {
        path.trim_end_matches(['/', '\\'])
    } else {
        path.trim_end_matches('/')
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-972
// ```c++
// string CDirEntry::NormalizePath(const string& path, EFollowLinks follow_links)
// {
//     if ( path.empty() ) {
//         return path;
//     }
//     static const char kSep[] = { DIR_SEPARATOR, '\0' };
//
//     std::list<string> head;       // already resolved to our satisfaction
//     std::list<string> tail;       // to resolve afterwards
//     string  current;              // to resolve next
//     int     link_depth = 0;
//
//     // Delete trailing slash for all paths except similar to 'd:\'
// #  ifdef DISK_SEPARATOR
//     if ( path.find(DISK_SEPARATOR) == NPOS ) {
//         current = DeleteTrailingPathSeparator(path);
//     } else {
//         current = path;
//     }
// #  else
//     current = DeleteTrailingPathSeparator(path);
// #  endif
//
//     if ( current.empty() ) {
//         // root dir
//         return string(1, DIR_SEPARATOR);
//     }
//
// #ifdef NCBI_OS_MSWIN
//     NStr::ReplaceInPlace(current, "/", "\\");
// #endif
//
//     while ( !current.empty()  ||  !tail.empty() ) {
//         std::list<string> pretail;
//         if ( !current.empty() ) {
//             NStr::Split(current, kSep, pretail);
//             current.erase();
//             if (pretail.front().empty()
// #ifdef DISK_SEPARATOR
//                 ||  pretail.front().find(DISK_SEPARATOR) != NPOS
// #endif
//                 ) {
//                 // Absolute path
//                 head.clear();
// #ifdef NCBI_OS_MSWIN
//                 // Remove leading "\\?\". Replace leading "\\?\UNC\" with "\\".
//                 static const char* const kUNC[] = { "", "", "?", "UNC" };
//                 std::list<string>::iterator it = pretail.begin();
//                 unsigned int matched = 0;
//                 while (matched < 4  &&  it != pretail.end()
//                        &&  !NStr::CompareNocase(*it, kUNC[matched])) {
//                     ++it;
//                     ++matched;
//                 }
//                 pretail.erase(pretail.begin(), it);
//                 switch (matched) {
//                 case 2: case 4: // got a UNC path (\\... or \\?\UNC\...)
//                     head.push_back(kEmptyStr);
//                     // fall through
//                 case 1:         // normal volume-less absolute path
//                     head.push_back(kEmptyStr);
//                     break;
//                 case 3/*?*/:    // normal path, absolute or relative
//                     break;
//                 }
// #endif
//             }
//             tail.splice(tail.begin(), pretail);
//         }
//
//         string next;
//         if (!tail.empty()) {
//             next = tail.front();
//             tail.pop_front();
//         }
//         if ( !head.empty() ) { // empty heads should accept anything
//             string& last = head.back();
//             if (last == DIR_CURRENT) {
//                 if (!next.empty()) {
//                     head.pop_back();
//                 }
// #ifdef DISK_SEPARATOR
//             } else if (!last.empty() && last[last.size()-1] == DISK_SEPARATOR) {
//                 // Allow almost anything right after a volume specification
// #endif
//             } else if (next == DIR_CURRENT) {
//                 // Leave out, since we already have content
//                 continue;
//             } else if (next.empty()) {
//                 continue; // leave out empty components in most cases
//             } else if (next == DIR_PARENT) {
// #ifdef DISK_SEPARATOR
//                 SIZE_TYPE pos;
// #endif
//                 // Back up if possible, assuming existing path to be "physical"
//                 if (last.empty()) {
//                     // Already at the root; .. is a no-op
//                     continue;
// #ifdef DISK_SEPARATOR
//                 } else if ((pos = last.find(DISK_SEPARATOR) != NPOS)) {
//                     last.erase(pos + 1);
// #endif
//                 } else if (last != DIR_PARENT) {
//                     head.pop_back();
//                     continue;
//                 }
//             }
//         }
// #ifdef NCBI_OS_UNIX
//         // Is there a Windows equivalent for readlink?
//         if ( follow_links ) {
//             string s(head.empty() ? next : NStr::Join(head, string(1, DIR_SEPARATOR)) + DIR_SEPARATOR + next);
//             char buf[PATH_MAX];
//             int  length = (int)readlink(s.c_str(), buf, sizeof(buf));
//             if (length > 0) {
//                 current.assign(buf, length);
//                 if (++link_depth >= 1024) {
//                     ERR_POST_X(1, Warning << "CDirEntry::NormalizePath(): "
//                                "Reached symlink depth limit " <<
//                                link_depth << " when resolving " << path);
//                     CNcbiError::Set(CNcbiError::eTooManySymbolicLinkLevels);
//                     follow_links = eIgnoreLinks;
//                 }
//                 continue;
//             }
//         }
// #endif
//         // Normal case: just append the next element to head
//         head.push_back(next);
//     }
//
//     // Special cases
//     if ( (head.size() == 0)  ||
//          (head.size() == 2  &&  head.front() == DIR_CURRENT  &&  head.back().empty()) ) {
//         // current dir
//         return DIR_CURRENT;
//     }
//     if (head.size() == 1  &&  head.front().empty()) {
//         // root dir
//         return string(1, DIR_SEPARATOR);
//     }
// #ifdef DISK_SEPARATOR
//     if (head.front().find(DISK_SEPARATOR) != NPOS) {
//         if ((head.size() == 2  &&  head.back() == DIR_CURRENT) ||
//             (head.size() == 3  &&  *(++head.begin()) == DIR_CURRENT  &&  head.back().empty()) ) {
//             // root dir on drive X:
//             return head.front() + DIR_SEPARATOR;
//         }
//     }
// #endif
//     // Compose path
//     return NStr::Join(head, string(1, DIR_SEPARATOR));
// }
// ```
// `DIR_SEPARATOR` is '/' on Unix and '\\' on Windows, where `DISK_SEPARATOR`
// is ':' (ncbifile.cpp:102-116). `NStr::Split` without flags keeps empty
// pieces. NCBI applies the function to a whole value such as `BLASTDB`, the
// separators of a list included, so `/a/x:/../b` becomes `/a/b`.
//
// The symbolic links are read as NCBI reads them, one path component at a
// time and only outside Windows; a relative path is read relative to the
// current directory of the process. The warning NCBI posts at the depth
// limit is not printed. A link whose target is not valid UTF-8 is an
// explicit error (LOSAT keeps path names as UTF-8 strings).
/// `CDirEntry::NormalizePath(path, eFollowLinks)`.
fn normalize_path(path: &str) -> Result<String, SeqDbError> {
    if path.is_empty() {
        return Ok(String::new());
    }
    let windows = cfg!(windows);
    let sep = MAIN_SEPARATOR;
    let disk_separator = ':';

    let mut head: Vec<String> = Vec::new();
    let mut tail: VecDeque<String> = VecDeque::new();
    let mut link_depth = 0;
    let mut follow_links = true;

    let mut current = if windows && path.contains(disk_separator) {
        path.to_string()
    } else {
        delete_trailing_path_separator(path).to_string()
    };
    if current.is_empty() {
        return Ok(sep.to_string());
    }
    if windows {
        current = current.replace('/', "\\");
    }

    while !current.is_empty() || !tail.is_empty() {
        if !current.is_empty() {
            let mut pretail: VecDeque<String> = current.split(sep).map(str::to_string).collect();
            current.clear();
            let absolute = pretail.front().is_some_and(|front| {
                front.is_empty() || (windows && front.contains(disk_separator))
            });
            if absolute {
                head.clear();
                if windows {
                    const UNC: [&str; 4] = ["", "", "?", "UNC"];
                    let mut matched = 0;
                    while matched < 4
                        && pretail
                            .front()
                            .is_some_and(|front| front.eq_ignore_ascii_case(UNC[matched]))
                    {
                        pretail.pop_front();
                        matched += 1;
                    }
                    match matched {
                        2 | 4 => {
                            head.push(String::new());
                            head.push(String::new());
                        }
                        1 => head.push(String::new()),
                        _ => {}
                    }
                }
            }
            while let Some(piece) = pretail.pop_back() {
                tail.push_front(piece);
            }
        }

        let next = tail.pop_front().unwrap_or_default();
        if let Some(last) = head.last_mut() {
            if last == "." {
                if !next.is_empty() {
                    head.pop();
                }
            } else if windows && last.ends_with(disk_separator) {
                // Allow almost anything right after a volume specification
            } else if next == "." || next.is_empty() {
                continue;
            } else if next == ".." {
                if last.is_empty() {
                    continue;
                } else if windows && last.contains(disk_separator) {
                    // `pos` receives the result of the comparison, 1, so
                    // `last.erase(pos + 1)` keeps the first two characters.
                    if last.is_char_boundary(2) {
                        last.truncate(2);
                    }
                } else if last != ".." {
                    head.pop();
                    continue;
                }
            }
        }
        if follow_links && !windows {
            let link = if head.is_empty() {
                next.clone()
            } else {
                format!("{}{sep}{next}", head.join(&sep.to_string()))
            };
            if let Ok(target) = std::fs::read_link(&link) {
                let target =
                    target
                        .into_os_string()
                        .into_string()
                        .map_err(|_| SeqDbError::Unsupported {
                            path: Path::new(&link).to_path_buf(),
                            message: "the target of the symbolic link is not valid UTF-8"
                                .to_string(),
                        })?;
                if !target.is_empty() {
                    current = target;
                    link_depth += 1;
                    if link_depth >= 1024 {
                        follow_links = false;
                    }
                    continue;
                }
            }
        }
        head.push(next);
    }

    if head.is_empty() || (head.len() == 2 && head[0] == "." && head[1].is_empty()) {
        return Ok(".".to_string());
    }
    if head.len() == 1 && head[0].is_empty() {
        return Ok(sep.to_string());
    }
    if windows
        && head[0].contains(disk_separator)
        && ((head.len() == 2 && head[1] == ".")
            || (head.len() == 3 && head[1] == "." && head[2].is_empty()))
    {
        return Ok(format!("{}{sep}", head[0]));
    }
    Ok(head.join(&sep.to_string()))
}

/// The separator of the directories of a search path: `s_GetPathSplitter`
/// (seqdbcommon.cpp:270-281), ":" on Unix and ";" elsewhere.
fn path_splitter() -> char {
    if cfg!(windows) {
        ';'
    } else {
        ':'
    }
}

/// The directories where NCBI BLAST+ looks for a database name: the current
/// directory, the environment variable `BLASTDB`, and the registry value
/// `[BLAST] BLASTDB` (`.ncbirc`, the program's `.ini` file, or the variable
/// `NCBI_CONFIG__BLAST__BLASTDB`), in this order.
///
/// A relative directory of `BLASTDB` or of the registry value is used as it
/// is written, so the operating system resolves it against the current
/// directory of the process, as for NCBI.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SearchPath {
    /// `CDir::GetCwd()`: where a volume named in an alias file is looked for
    /// first ("local copy").
    cwd: String,
    /// The string of `CSeqDBAtlas::GenerateSearchPath`.
    path: String,
}

impl SearchPath {
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
    // ```c++
    //     static const string GenerateSearchPath() {
    //         string splitter;
    //         string path;
    // #if defined(NCBI_OS_UNIX)
    //         splitter = ":";
    // #else
    //         splitter = ";";
    // #endif
    //         // Local directory first;
    //         path  = CDirEntry::NormalizePath(CDir::GetCwd(),eFollowLinks);
    //         path += splitter;
    //         // Then, BLASTDB;
    //         CNcbiEnvironment env;
    //         path += CDirEntry::NormalizePath(env.Get("BLASTDB"),eFollowLinks);
    //         path += splitter;
    //         // Finally, the config file.
    //         CNcbiApplicationAPI* app = CNcbiApplicationAPI::Instance();
    //         if (app) {
    //             const CNcbiRegistry& registry = app->GetConfig();
    //             if (registry.HasEntry("BLAST", "BLASTDB")) {
    //                 path += CDirEntry::NormalizePath(registry.Get("BLAST", "BLASTDB"),eFollowLinks);
    //                 path += splitter;
    //             }
    //         }
    //         return path;
    //     }
    // ```
    // `env.Get` gives the empty string for a variable that is not set
    // (ncbienv.cpp:109-130), and `NormalizePath` returns an empty string as
    // it is, so the string is `CWD::` without `BLASTDB` and without a
    // registry entry.
    /// The search path of a process whose current directory is `cwd`, whose
    /// environment variable `BLASTDB` is `blastdb` (`None`: not set) and
    /// whose registry has the entry `[BLAST] BLASTDB` with the value
    /// `registry` (`None`: no such entry;
    /// `blastinput::ncbi_environment::ApplicationSettings::blastdb`).
    ///
    /// A value that is not valid UTF-8 is an explicit error: LOSAT keeps the
    /// names of database files as UTF-8 strings.
    pub fn new(
        cwd: &Path,
        blastdb: Option<&OsStr>,
        registry: Option<&str>,
    ) -> Result<Self, SeqDbError> {
        let utf8 = |value: &OsStr, what: &str| {
            value
                .to_str()
                .map(str::to_string)
                .ok_or_else(|| SeqDbError::Unsupported {
                    path: Path::new(value).to_path_buf(),
                    message: format!("{what} is not valid UTF-8"),
                })
        };
        let cwd = utf8(cwd.as_os_str(), "the current directory")?;
        let blastdb = match blastdb {
            Some(value) => utf8(value, "the environment variable BLASTDB")?,
            None => String::new(),
        };
        let splitter = path_splitter();
        let mut path = normalize_path(&cwd)?;
        path.push(splitter);
        path.push_str(&normalize_path(&blastdb)?);
        path.push(splitter);
        if let Some(registry) = registry {
            path.push_str(&normalize_path(registry)?);
            path.push(splitter);
        }
        Ok(Self { cwd, path })
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3708-3716
    // ```c++
    // string CDir::GetCwd(void)
    // {
    //     TXChar buf[4096];
    //     if ( NcbiSys_getcwd(buf, sizeof(buf)/sizeof(TXChar) - 1) ) {
    //         return _T_CSTRING(buf);
    //     }
    //     LOG_ERROR_ERRNO(90, "CDir::GetCwd(): Cannot get current directory");
    //     return kEmptyCStr;
    // }
    // ```
    /// The search path of this process: its current directory (the empty
    /// string when it cannot be read, as `CDir::GetCwd` gives), its
    /// environment variable `BLASTDB`, and `registry` as in [`Self::new`].
    pub fn from_environment(registry: Option<&str>) -> Result<Self, SeqDbError> {
        let cwd = std::env::current_dir().unwrap_or_default();
        let blastdb = std::env::var_os("BLASTDB");
        Self::new(&cwd, blastdb.as_deref(), registry)
    }

    /// The string NCBI names in its "No alias or index file found" message:
    /// the directories, each followed by the separator.
    pub fn as_str(&self) -> &str {
        &self.path
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:330-374
    // ```c++
    // static string s_SeqDB_TryPaths(const string         & blast_paths,
    //                                const string         & dbname,
    //                                char                   dbtype,
    //                                bool                   exact,
    //                                CSeqDB_FileExistence & access,
    //                                bool                   linkoutdb_search = false)
    // {
    //     ...
    //     vector<string> roads;
    //     NStr::Split(blast_paths, s_GetPathSplitter(), roads, NStr::fSplit_Tokenize);
    //
    //     string result;
    //     string attempt;
    //
    //     ITERATE(vector<string>, road, roads) {
    //         attempt.erase();
    //
    //         SeqDB_CombinePath(CSeqDB_Substring(SeqDB_MakeOSPath(*road)),
    //                           CSeqDB_Substring(dbname),
    //                           0,
    //                           attempt);
    //
    //         if (exact) {
    //             if (access.DoesFileExist(attempt)) {
    //                 result = attempt;
    //                 break;
    //             }
    //         } else {
    //             if (s_SeqDB_DBExists(attempt, dbtype, access, linkoutdb_search)) {
    //                 result = attempt;
    //                 break;
    //             }
    //         }
    //     }
    //
    //     return result;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:232-262
    // ```c++
    // static bool s_SeqDB_DBExists(const string         & dbname,
    //                              char                   dbtype,
    //                              CSeqDB_FileExistence & access,
    //                              bool linkoutdb_search)
    // {
    //     ...
    //         path.append(".-al");
    //
    //         path[path.size()-3] = dbtype;
    //
    //         if (access.DoesFileExist(path)) {
    //             return true;
    //         }
    //
    //         path[path.size()-2] = 'i';
    //         path[path.size()-1] = 'n';
    //
    //         if (access.DoesFileExist(path)) {
    //             return true;
    // ```
    // `fSplit_Tokenize` is `fSplit_MergeDelimiters | fSplit_Truncate`
    // (ncbistr.hpp:2500-2510): empty directories are left out. With
    // `seqtype` the name is a database name and `s_SeqDB_DBExists` looks for
    // its alias file, then for its index file, inside one directory before
    // the next directory is tried; without it the name is a file name
    // ("exact"). A file exists when it is a regular file (`exists`).
    /// The first directory of the path that holds `dbname`, joined with it.
    fn find(&self, dbname: &str, seqtype: Option<SeqType>) -> Option<String> {
        self.path
            .split(path_splitter())
            .filter(|road| !road.is_empty())
            .map(|road| combine_path(&fix_delimiters(road), dbname))
            .find(|attempt| match seqtype {
                Some(seqtype) => {
                    let base = Path::new(attempt);
                    exists(&with_extension(base, seqtype, "al"))
                        || exists(&with_extension(base, seqtype, "in"))
                }
                None => exists(Path::new(attempt)),
            })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:453-462
    // ```c++
    // string SeqDB_ResolveDbPath(const string & filename)
    // {
    //     CSeqDB_SimpleAccessor access;
    //
    //     return s_SeqDB_FindBlastDBPath(filename,
    //                                    '-',
    //                                    0,
    //                                    true,
    //                                    access);
    // }
    // ```
    // `exact` is true, so `s_SeqDB_TryPaths` takes the file name as it is
    // (the `dbtype` `'-'` is not used) and returns the first directory of the
    // search path in which `CFile::GetLength` finds a regular file of that
    // name; the search path is `CSeqDBAtlas::GenerateSearchPath()` (the
    // string of this `SearchPath`).
    /// `SeqDB_ResolveDbPath`: the first directory of the path that holds a
    /// regular file `filename`, joined with it (`None`: NCBI's empty string).
    pub fn resolve_file(&self, filename: &str) -> Option<String> {
        self.find(filename, None)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1762-1807
// ```c++
// void SeqDB_SplitQuoted(const string             & dbname,
//                        vector<CSeqDB_Substring> & dbs,
//                        bool keep_quote)
// {
//     // split names
//
//     const char * sp = dbname.data();
//
//     bool quoted = false;
//     unsigned begin = 0;
//
//     for(unsigned i = 0; i < dbname.size(); i++) {
//         char ch = dbname[i];
//
//         if (quoted) {
//             // Quoted mode sees '"' as the only actionable token.
//             if (ch == '"') {
//                 if (begin < i) {
//                 	if(keep_quote) i++;
//                     dbs.push_back(CSeqDB_Substring(sp + begin, sp + i));
//                 }
//                 begin = i + 1;
//                 quoted = false;
//             }
//         } else {
//             // Non-quote mode: Space or quote starts the next string.
//
//             if (ch == ' ') {
//                 if (begin < i) {
//                     dbs.push_back(CSeqDB_Substring(sp + begin, sp + i));
//                 }
//                 begin = i + 1;
//             } else if (ch == '"') {
//                 if (begin < i) {
//                     dbs.push_back(CSeqDB_Substring(sp + begin, sp + i));
//                 }
//                 begin = keep_quote ? i : i + 1;
//                 quoted = true;
//             }
//         }
//     }
//
//     if (begin < dbname.size()) {
//         dbs.push_back(CSeqDB_Substring(sp + begin, sp + dbname.size()));
//     }
// }
// ```
// `keep_quote` is false for a `DBLIST` (seqdbalias.cpp:85-88).
fn split_quoted(dbname: &[u8]) -> Vec<&[u8]> {
    let mut dbs = Vec::new();
    let mut quoted = false;
    let mut begin = 0;
    for (i, &ch) in dbname.iter().enumerate() {
        if quoted {
            if ch == b'"' {
                if begin < i {
                    dbs.push(&dbname[begin..i]);
                }
                begin = i + 1;
                quoted = false;
            }
        } else if ch == b' ' {
            if begin < i {
                dbs.push(&dbname[begin..i]);
            }
            begin = i + 1;
        } else if ch == b'"' {
            if begin < i {
                dbs.push(&dbname[begin..i]);
            }
            begin = i + 1;
            quoted = true;
        }
    }
    if begin < dbname.len() {
        dbs.push(&dbname[begin..]);
    }
    dbs
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:689-710
// ```c++
//     while(p < ep) {
//         // Skip spaces
//         while((p < ep) && (*p == ' ')) {
//             p++;
//         }
//
//         const char * eolp = p;
//
//         while((eolp < ep) && !SEQDB_ISEOL(*eolp)) {
//             eolp++;
//         }
//
//         // Non-empty line, so read it.
//         if (eolp != p) {
//             x_ReadLine(p, eolp, name_s, value_s);
//         }
//
//         p = eolp + 1;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:406-436
// ```c++
//     const char * p = bp;
//
//     // If first nonspace char is '#', line is a comment, so skip.
//     if (*p == '#') {
//         return;
//     }
//
//     // Find name
//     const char * spacep = p;
//
//     while((spacep < ep) && ((*spacep != ' ') && (*spacep != '\t')))
//         spacep ++;
//
//     s_SeqDB_QuickAssign(name, p, spacep);
//
//     // Skip spaces, tabs, to find value
//     while((spacep < ep) && ((*spacep == ' ') || (*spacep == '\t')))
//         spacep ++;
//
//     // Strip spaces, tabs from end
//     while((spacep < ep) && ((ep[-1] == ' ') || (ep[-1] == '\t')))
//         ep --;
//
//     s_SeqDB_QuickAssign(value, spacep, ep);
//
//     for(size_t i = 0; i<value.size(); i++) {
//         if (value[i] == '\t') {
//             value[i] = ' ';
//         }
//     }
// ```
// `SEQDB_ISEOL(x)` is `(x) == '\n' || (x) == '\r'` (seqdbgeneral.hpp:192). A
// key that appears again replaces the earlier value (seqdbalias.cpp:446-449).
/// The key / value lines of an alias file, in file order.
fn read_values(bytes: &[u8]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut values: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let ep = bytes.len();
    let mut p = 0;
    while p < ep {
        while p < ep && bytes[p] == b' ' {
            p += 1;
        }
        let mut eolp = p;
        while eolp < ep && bytes[eolp] != b'\n' && bytes[eolp] != b'\r' {
            eolp += 1;
        }
        if eolp != p && bytes[p] != b'#' {
            let line = &bytes[p..eolp];
            let mut spacep = 0;
            while spacep < line.len() && line[spacep] != b' ' && line[spacep] != b'\t' {
                spacep += 1;
            }
            let name = &line[..spacep];
            while spacep < line.len() && (line[spacep] == b' ' || line[spacep] == b'\t') {
                spacep += 1;
            }
            let mut end = line.len();
            while spacep < end && (line[end - 1] == b' ' || line[end - 1] == b'\t') {
                end -= 1;
            }
            let value: Vec<u8> = line[spacep..end]
                .iter()
                .map(|&c| if c == b'\t' { b' ' } else { c })
                .collect();
            if !name.is_empty() {
                match values.iter_mut().find(|(key, _)| key == name) {
                    Some(entry) => entry.1 = value,
                    None => values.push((name.to_vec(), value)),
                }
            }
        }
        p = eolp + 1;
    }
    values
}

struct Resolver<'a> {
    seqtype: SeqType,
    /// Alias files being expanded, outermost first.
    recurse: Vec<String>,
    /// The directory where a volume named in an alias file is looked for
    /// before the alias file's own directory: NCBI's "local copy", in the
    /// current directory. `None`: the alias file's directory only.
    local_dir: Option<String>,
    /// NCBI's search path, when the name is resolved as a search program
    /// resolves it (`resolve_in`). `None`: the name is the path of one
    /// database (`resolve_one_name`), and names are looked for where they
    /// point only.
    search: Option<&'a SearchPath>,
}

/// Where the volumes that an alias file names are looked for.
#[derive(Clone, Copy)]
pub(super) enum LocalCopy<'a> {
    /// In the current directory, then in the alias file's (NCBI's rule).
    CurrentDir,
    /// In this directory, then in the alias file's: the rule above with the
    /// directory given by the caller.
    #[cfg_attr(not(test), allow(dead_code))]
    Dir(&'a str),
    /// In the alias file's directory only.
    No,
}

impl Resolver<'_> {
    fn alias_path(&self, base: &str) -> std::path::PathBuf {
        with_extension(Path::new(base), self.seqtype, "al")
    }

    fn index_path(&self, base: &str) -> std::path::PathBuf {
        with_extension(Path::new(base), self.seqtype, "in")
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:851-862
    // ```c++
    //         if (found) {
    //             string normal_name;
    //             if (m_ExpandLinks) {
    //                 // normalize this_name
    //                 string dir_name, base_name;
    //                 bp.FindDirName().GetString(dir_name);
    //                 bp.FindBaseName().GetString(base_name);
    //                 normal_name = CDirEntry::NormalizePath(dir_name, eFollowLinks) +
    //                   CDirEntry::GetPathSeparator() + base_name;
    //             } else {
    //                 normal_name = bp.GetBasePathS();
    //             }
    // ```
    // `m_ExpandLinks` is true for the alias file of a database that is opened
    // (seqdbalias.hpp:1005-1008, the default of `CSeqDBAliasFile`;
    // seqdbimpl.cpp:56 passes no value). The same lines stand at the three
    // other places where NCBI records a volume (seqdbalias.cpp:777-786,
    // 894-903) or a resolved name (seqdbalias.cpp:325-335).
    /// The name under which NCBI records the volume `base`: its directory
    /// normalised, with symbolic links followed. `resolve_one_name` keeps
    /// the name as it is.
    fn normal_name(&self, base: &str) -> Result<String, SeqDbError> {
        if self.search.is_none() {
            return Ok(base.to_string());
        }
        Ok(format!(
            "{}{MAIN_SEPARATOR}{}",
            normalize_path(dir_name(base))?,
            base_name(base)
        ))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:203-222
    // ```c++
    // bool
    // CSeqDBAliasSets::FindAliasPath(const CSeqDB_Path & dbpath,
    //                                CSeqDB_Path       * resolved,
    //                                CSeqDBLockHold    & locked)
    // {
    //     CSeqDB_Path     aset_path;
    //     CSeqDB_FileName alias_fname;
    //
    //     x_DbToIndexName(dbpath, aset_path, alias_fname);
    //
    //     CSeqDB_Path resolved_aset;
    //
    //     if (! FindBlastDBPath(aset_path, resolved_aset)) {
    //         return false;
    //     }
    //
    //     CSeqDB_Path afpath(resolved_aset.FindDirName(),
    //                        alias_fname.GetFileNameSub());
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:453-460
    // ```c++
    // void CSeqDBAliasSets::x_DbToIndexName(const CSeqDB_Path & dbpath,
    //                                       CSeqDB_Path       & index_path,
    //                                       CSeqDB_FileName   & alias_fname)
    // {
    //     index_path.ReplaceFilename(dbpath,
    //                                CSeqDB_Substring(kSeqDBGroupAliasFileName));
    //     alias_fname.Assign(dbpath.FindFileName());
    // }
    // ```
    // `kSeqDBGroupAliasFileName` is "index.alx" (seqdbcommon.cpp:48), and
    // `FindBlastDBPath` of a `CSeqDB_Path` looks for that file name in the
    // directories of the search path (seqdbalias.hpp:285-300, "exact").
    // NCBI looks a name up in the group alias file before the alias file or
    // the volume of that name, and goes on when the group file has no entry
    // for it. LOSAT does not read group alias files, so the group file NCBI
    // would read is an explicit error whatever it holds.
    /// Rejects the group alias file that NCBI would consult for the alias
    /// file or volume `base`.
    fn reject_group_alias(&self, base: &str) -> Result<(), SeqDbError> {
        let aset_path = combine_path(dir_name(base), "index.alx");
        let found = match self.search {
            Some(search) => search.find(&aset_path, None),
            None => exists(Path::new(&aset_path)).then_some(aset_path),
        };
        match found {
            Some(group) => Err(SeqDbError::Unsupported {
                path: Path::new(&group).to_path_buf(),
                message: "the directory has a group alias file; LOSAT does not read group \
                          alias files"
                    .to_string(),
            }),
            None => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:169-176
    // ```c++
    //     recurse.Push(m_ThisName);
    //
    //     x_ReadValues(m_ThisName, locked);
    //     x_Tokenize(m_Values["DBLIST"]);
    //
    //     CSeqDB_BasePath basepath(m_ThisName.FindBasePath());
    //
    //     x_ExpandAliases(basepath, prot_nucl, recurse, locked);
    // ```
    /// Reads the alias file of `base` and expands its `DBLIST`.
    fn alias_node(&mut self, base: &str) -> Result<AliasNode, SeqDbError> {
        let path = self.alias_path(base);
        let unsupported = |message: String| SeqDbError::Unsupported {
            path: path.clone(),
            message,
        };
        // NCBI maps the alias file (`CSeqDBFileMemMap lease(m_Atlas,path.GetPathS())`,
        // seqdbalias.cpp:662-667); a file that cannot be opened gives the
        // errors of `DbFile::open`.
        let bytes = DbFile::open(&path)?.read_all()?;
        let values = read_values(&bytes);
        let mut title = None;
        let mut nseq = None;
        let mut length = None;
        let mut dblist: &[u8] = &[];
        for (key, value) in &values {
            match key.as_slice() {
                b"TITLE" => title = Some(value.clone()),
                b"NSEQ" => nseq = Some(value.clone()),
                b"LENGTH" => length = Some(value.clone()),
                b"DBLIST" => dblist = value,
                other => {
                    if let Some(key) = UNSUPPORTED_KEYS.iter().find(|k| k.as_bytes() == other) {
                        return Err(unsupported(format!(
                            "alias file key {key} is not supported; LOSAT reads TITLE, DBLIST, NSEQ and \
                             LENGTH"
                        )));
                    }
                }
            }
        }
        let mut entries = Vec::new();
        for entry in split_quoted(dblist) {
            let entry = std::str::from_utf8(entry)
                .map_err(|_| unsupported("a DBLIST name is not valid UTF-8".to_string()))?;
            entries.push(fix_delimiters(entry));
        }

        self.recurse.push(path.to_string_lossy().into_owned());
        let mut node = AliasNode {
            this_name: Some(path.to_string_lossy().into_owned()),
            title,
            nseq,
            length,
            ..AliasNode::default()
        };
        let result = self.expand(&mut node, dir_name(base), base, &entries, false);
        self.recurse.pop();
        result.map(|()| node)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:740-926
    // ```c++
    //     if (m_DBList.empty()) {
    //         ...
    //         NCBI_THROW(CSeqDBException,
    //                    eArgErr,
    //                    string("No database names were ") + situation);
    //     }
    //
    //     for(size_t i = 0; i < m_DBList.size(); i++) {
    //         ...
    //         if (m_DBList[i].FindDirName().Empty()) {
    //             if (m_DBList[i].FindBaseName() == this_name.FindBaseName()) {
    //                 ...
    //                 if (!found) {
    //                     m_VolNames.push_back(CSeqDB_BasePath(normal_name));
    //                 }
    //                 continue;
    //             }
    //         }
    //         ...
    //         CSeqDB_BasePath base(m_DBPath, m_DBList[i]);
    //         CSeqDB_Path new_db_path( base, prot_nucl, 'a', 'l' );
    //
    //         if ( recurse.Exists(new_db_path) ) {
    //             NCBI_THROW(CSeqDBException,
    //                        eFileErr,
    //                        "Illegal configuration: DB alias files are mutually recursive.");
    //         }
    //         ...
    //         if ( m_AliasSets.FindAliasPath(new_db_path, 0, locked) ||
    //              m_Atlas.DoesFileExist(new_db_path) ) {
    //
    //             x_AppendSubNode(base, prot_nucl, recurse, locked);
    //             continue;
    //         }
    //         ...
    //         // Always check local copy first unless full path is given
    //         if (!m_SkipLocal[i]) {
    //             // GCC 3.0.4 requires the extra parentheses below.
    //             CSeqDB_BasePath local_base((CSeqDB_DirName(CDir::GetCwd())),
    //                                    CSeqDB_BaseName(m_DBList[i].FindBaseName()));
    //             CSeqDB_Path new_local_vol_path(local_base, prot_nucl, 'i', 'n' );
    //             if (m_Atlas.DoesFileExist(new_local_vol_path)) {
    //                 bp = CSeqDB_BasePath(new_local_vol_path.FindBasePath());
    //                 found = true;
    //             }
    //         }
    //
    //         if (!found) {
    //             CSeqDB_Path new_vol_path( base, prot_nucl, 'i', 'n' );
    //             if (m_Atlas.DoesFileExist(new_vol_path)) {
    //                 bp = CSeqDB_BasePath(new_db_path.FindBasePath() );
    //                 found = true;
    //             }
    //         }
    //         ...
    //         ostringstream oss;
    //         oss << "Could not find volume or alias file ("
    //             << m_DBList[i].GetBasePathS() << ") referenced in alias file ("
    //             << this_name.GetBasePathS() << ").";
    //
    //         NCBI_THROW(CSeqDBException, eFileErr, oss.str());
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:877-917
    // ```c++
    //         // If all that failed, "restart" the search using the blast DB
    //         // path configuration.  This ensures always finding the local
    //         // copy of a database even if the alias file resides on remote site.
    //
    //         CSeqDB_BasePath result;
    //         CSeqDB_BasePath restart(m_DBList[i]);
    //
    //         if (m_AliasSets.FindBlastDBPath(restart, prot_nucl, result)) {
    //             // either alias or volume file exists for this entry
    //             CSeqDB_Path new_alias( result, prot_nucl, 'a', 'l' );
    //             CSeqDB_Path new_volume( result, prot_nucl, 'i', 'n' );
    //
    //             if (m_Atlas.DoesFileExist(new_alias)) {
    //                 x_AppendSubNode( result, prot_nucl, recurse, locked );
    //             } else if (m_Atlas.DoesFileExist(new_volume)) {
    //                 string normal_name;
    //                 if (m_ExpandLinks) {
    //                     ...
    //                 }
    //                 ...
    //                 if (!found) {
    //                     m_VolNames.push_back(CSeqDB_BasePath(normal_name));
    //                 }
    //             }
    //             continue;
    //         }
    // ```
    // `resolve_one_name` (no search path) does not take this last step and
    // keeps the names as they are written; `resolve_in` follows NCBI: the
    // name as the alias file writes it is looked up in the directories of the
    // search path, and a volume is recorded under its normalised name
    // (`normal_name`).
    //
    // NCBI does not test an alias file found by this last step against the
    // files being expanded, so alias files that name each other through the
    // search path make it recurse without end. LOSAT reports them as the
    // mutual recursion they are, as an unsupported database.
    fn expand(
        &mut self,
        node: &mut AliasNode,
        db_path: &str,
        this_name: &str,
        db_list: &[String],
        skip_local: bool,
    ) -> Result<(), SeqDbError> {
        if db_list.is_empty() {
            return Err(SeqDbError::Format {
                path: self.alias_path(this_name),
                message: format!("No database names were found in alias file [{this_name}]."),
            });
        }
        for entry in db_list {
            if dir_name(entry).is_empty() && base_name(entry) == base_name(this_name) {
                push_unique(&mut node.vol_names, self.normal_name(this_name)?);
                continue;
            }

            let base = combine_path(db_path, entry);
            let new_db_path = self.alias_path(&base);
            let new_db_name = new_db_path.to_string_lossy().into_owned();
            if self.recurse.contains(&new_db_name) {
                return Err(SeqDbError::Format {
                    path: new_db_path,
                    message: "Illegal configuration: DB alias files are mutually recursive."
                        .to_string(),
                });
            }
            // NCBI looks a name up in the group alias file of its directory
            // before the alias file itself (`FindAliasPath` above,
            // seqdbalias.cpp:206-247); LOSAT does not read group alias files.
            self.reject_group_alias(&base)?;
            if exists(&new_db_path) {
                let subnode = self.alias_node(&base)?;
                node.subnodes.push(subnode);
                continue;
            }

            let mut found = None;
            if !skip_local {
                if let Some(local_dir) = &self.local_dir {
                    let local_base = combine_path(local_dir, base_name(entry));
                    if exists(&self.index_path(&local_base)) {
                        found = Some(local_base);
                    }
                }
            }
            if found.is_none() && exists(&self.index_path(&base)) {
                found = Some(base);
            }
            if let Some(volume) = found {
                push_unique(&mut node.vol_names, self.normal_name(&volume)?);
                continue;
            }

            if let Some(result) = self.search.and_then(|s| s.find(entry, Some(self.seqtype))) {
                let new_alias = self.alias_path(&result);
                if exists(&new_alias) {
                    if self
                        .recurse
                        .contains(&new_alias.to_string_lossy().into_owned())
                    {
                        return Err(SeqDbError::Unsupported {
                            path: new_alias,
                            message: "the alias files name each other through the database \
                                      search path (NCBI BLAST+ does not end on them)"
                                .to_string(),
                        });
                    }
                    let subnode = self.alias_node(&result)?;
                    node.subnodes.push(subnode);
                } else if exists(&self.index_path(&result)) {
                    push_unique(&mut node.vol_names, self.normal_name(&result)?);
                }
                continue;
            }

            return Err(SeqDbError::NotFound(format!(
                "Could not find volume or alias file ({entry}) referenced in alias file \
                 ({this_name})."
            )));
        }
        Ok(())
    }
}

fn push_unique(names: &mut Vec<String>, name: String) {
    if !names.contains(&name) {
        names.push(name);
    }
}

fn collect_volumes(node: &AliasNode, out: &mut Vec<String>) {
    out.extend(node.vol_names.iter().cloned());
    for subnode in &node.subnodes {
        collect_volumes(subnode, out);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:190-197
// ```c++
// bool SeqDB_CompareVolume(const string & s1, const string & s2)
// {
//     string x1, x2;
//     CSeqDB_Path(s1).FindBaseName().GetString(x1);
//     CSeqDB_Path(s2).FindBaseName().GetString(x2);
//     if (x1 != x2) return (x1 < x2);
//     else return (s1 < s2);
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:932-958
// ```c++
//     set<string> volset;
//     set<string> aliset;
//
//     if (recursive) {
//         x_FindVolumePaths(volset, aliset);
//     } else {
//     ...
//     vols.clear();
//     ITERATE(set<string>, iter, volset) {
//         vols.push_back(*iter);
//     }
//
//     // Sort to insure deterministic order.
//     sort(vols.begin(), vols.end(), SeqDB_CompareVolume);
// ```
// The volumes are opened, and ordinal ids are given, in the sorted order of
// the volume names, not in `DBLIST` order.
//
// Two cases are rejected because LOSAT does not reproduce what NCBI derives
// from them: a volume reached more than once (NCBI's alias walk counts it
// once per mention, seqdbalias.cpp:1588-1592), and two volumes of the same
// base name in different directories.
/// The volumes of the alias tree `root` in NCBI's order.
fn volumes_in_order(root: AliasNode) -> Result<AliasTree, SeqDbError> {
    let mut mentions = Vec::new();
    collect_volumes(&root, &mut mentions);
    let mut volumes = mentions.clone();
    volumes.sort_by(|a, b| {
        path_base_name(a)
            .cmp(path_base_name(b))
            .then_with(|| a.cmp(b))
    });
    volumes.dedup();
    let unsupported = |volume: &str, message: &str| SeqDbError::Unsupported {
        path: Path::new(volume).to_path_buf(),
        message: message.to_string(),
    };
    if volumes.len() != mentions.len() {
        let repeated = volumes
            .iter()
            .find(|volume| mentions.iter().filter(|m| m == volume).count() > 1)
            .map_or("", String::as_str);
        return Err(unsupported(
            repeated,
            "the volume is listed more than once by the alias files",
        ));
    }
    if let Some(pair) = volumes.windows(2).find(|pair| {
        path_base_name(&pair[0]) == path_base_name(&pair[1])
            && dir_name(&pair[0]) != dir_name(&pair[1])
    }) {
        return Err(unsupported(
            &pair[1],
            "two volumes of the database have the same name in different directories",
        ));
    }
    Ok(AliasTree { root, volumes })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:155-165
// ```c++
// CSeqDB::CSeqDB(const string & dbname,
//                ESeqType       seqtype,
//                CSeqDBGiList * gi_list,
//                bool           use_atlas_lock)
//
// {
//     if (dbname.size() == 0) {
//         NCBI_THROW(CSeqDBException,
//                    eArgErr,
//                    "Database name is required.");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:249-337
// ```c++
// void CSeqDBAliasNode::x_ResolveNames(char prot_nucl, CSeqDBLockHold & locked)
// {
//     m_DBPath = CSeqDB_DirName(".");
//
//     size_t i = 0;
//
//     for(i = 0; i < m_DBList.size(); i++) {
//         // skip local DB search only if absolute path is given
//         if(m_DBList[i].GetBasePathS().find(CDirEntry::GetPathSeparator()) != string::npos) {
//             m_SkipLocal[i] = true;
//         }
//
//         const CSeqDB_Path db_path( (CSeqDB_BasePath(m_DBList[i])), prot_nucl, 'a', 'l' );
//
//         CSeqDB_Path resolved_path;
//
//         // search for X/kSeqDBGroupAliasFileName
//         if (! m_AliasSets.FindAliasPath(db_path, & resolved_path, locked)) {
//             CSeqDB_BasePath base(db_path.FindBasePath());
//             CSeqDB_BasePath resolved_bp;
//
//             // search for X/base.nal/nin
//             if (m_AliasSets.FindBlastDBPath(base, prot_nucl, resolved_bp)) {
//                 resolved_path = CSeqDB_Path(resolved_bp, prot_nucl, 'a', 'l');
//             }
//         }
//
//         if (! resolved_path.Valid()) {
//             string p_or_n;
//
//             switch(prot_nucl) {
//             case 'p':
//                 p_or_n = "protein";
//                 break;
//
//             case 'n':
//                 p_or_n = "nucleotide";
//                 break;
//             ...
//             }
//
//             // Do over (to get the search path).  This doesnt use the
//             // resolution map, since speed is not of the essence.
//
//             string search_path;
//             string input_path;
//
//             db_path.FindBasePath().GetString(input_path);
//
//             SeqDB_FindBlastDBPath(input_path,
//                                   prot_nucl,
//                                   & search_path,
//                                   false,
//                                   m_Atlas);
//
//             ostringstream oss;
//             oss << "No alias or index file found for " << p_or_n
//                 << " database [" << m_DBList[i].GetBasePathS()
//                 << "] in search path [" << search_path << "]. "
//                 << "Please verify the spelling of the BLAST "
//                 << "database and its molecule type.";
//
//             string msg(oss.str());
//
//             NCBI_THROW(CSeqDBException,
//                        eFileErr,
//                        msg);
//         } else {
//             // full dereferenced name but without suffix /X/base
//             if (m_ExpandLinks) {
//                 string dir_name, base_name;
//                 resolved_path.FindDirName().GetString(dir_name);
//                 resolved_path.FindBaseName().GetString(base_name);
//                 m_DBList[i].Assign(CSeqDB_Substring(
//                     CDirEntry::NormalizePath(dir_name, eFollowLinks) +
//                     CDirEntry::GetPathSeparator() +
//                     base_name  ));
//             } else {
//                 m_DBList[i].Assign(resolved_path.FindBasePath());
//             }
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:343-383
// ```c++
//     size_t common = m_DBList[0].GetBasePathS().size();
//
//     // Reduce common length to length of min db path.
//     for(i = 1; common && (i < m_DBList.size()); i++) {
//         ...
//     }
//
//     if (common) {
//         --common;
//     }
//     ...
//     // Adjust back to whole path component.
//     while(common && (first[common-1] != CFile::GetPathSeparator())) {
//         --common;
//     }
//
//     if (common > 1) {
//         // Factor out common path components.
//         m_DBPath.Assign( CSeqDB_Substring(first.data(), first.data() + common) );
//
//         for(i = 0; i < m_DBList.size(); i++) {
//             CSeqDB_Substring sub(m_DBList[i].GetBasePathS());
//             sub.EraseFront((int) common);
//
//             m_DBList[i].Assign(sub);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:112-127
// ```c++
//     m_Values["DBLIST"] = dbname_list;
//
//     x_Tokenize(dbname_list);
//     ...
//     x_ResolveNames(prot_nucl, locked);
//
//     CSeqDBAliasStack recurse;
//
//     x_ExpandAliases(CSeqDB_BasePath("-"), prot_nucl, recurse, locked);
// ```
// The name is a list (`x_Tokenize`, quoted at `resolve_one_name`): spaces
// separate databases and double quotes keep a name with spaces together.
// LOSAT searches one database, so a list of several names is an explicit
// error. A name carries no extension: `db.pin` is looked for as `db.pin.pal`
// and `db.pin.pin`. A name with a path separator is joined to every
// directory of the search path like any other (an absolute name stays as it
// is, `combine_path`); it only turns off the look for a volume of the same
// base name in the current directory (`m_SkipLocal`).
/// Resolves `dbname`, the value of a search program's `-db`, as NCBI BLAST+
/// resolves it: the name is looked for in the directories of `search`, an
/// alias file before an index file inside each directory.
pub(super) fn resolve_in(
    dbname: &str,
    seqtype: SeqType,
    search: &SearchPath,
) -> Result<AliasTree, SeqDbError> {
    if dbname.is_empty() {
        return Err(SeqDbError::NotFound(
            "Database name is required.".to_string(),
        ));
    }
    let names = split_quoted(dbname.as_bytes());
    let name = match names.as_slice() {
        // `x_ExpandAliases` on an empty list (quoted at `Resolver::expand`).
        [] => {
            return Err(SeqDbError::NotFound(
                "No database names were passed to CSeqDB::CSeqDB().".to_string(),
            ));
        }
        [name] => fix_delimiters(&String::from_utf8_lossy(name)),
        _ => {
            return Err(SeqDbError::Unsupported {
                path: Path::new(dbname).to_path_buf(),
                message: "the name lists several databases; LOSAT searches one database \
                          (a name with spaces goes between double quotes)"
                    .to_string(),
            });
        }
    };
    let mut resolver = Resolver {
        seqtype,
        recurse: Vec::new(),
        local_dir: Some(search.cwd.clone()),
        search: Some(search),
    };
    let skip_local = name.contains(MAIN_SEPARATOR);

    // `db_path` is the name with the extension of an alias file, and
    // `db_path.FindBasePath()` takes that extension off again.
    let db_path = resolver.alias_path(&name);
    let db_path = db_path.to_string_lossy();
    resolver.reject_group_alias(&db_path)?;
    let Some(resolved) = search.find(remove_extn(&db_path), Some(seqtype)) else {
        return Err(SeqDbError::NotFound(format!(
            "No alias or index file found for {} database [{name}] in search path [{}]. \
             Please verify the spelling of the BLAST database and its molecule type.",
            seqtype.name(),
            search.as_str()
        )));
    };
    // `resolved_path` is the resolved name with the extension of an alias
    // file; its directory is normalised and its base name loses the extension.
    let resolved_path = resolver.alias_path(&resolved);
    let resolved_path = resolved_path.to_string_lossy();
    let first = format!(
        "{}{MAIN_SEPARATOR}{}",
        normalize_path(dir_name(&resolved_path))?,
        remove_extn(base_name(&resolved_path))
    );

    // One name: the common prefix is its directory, separator included.
    let bytes = first.as_bytes();
    let mut common = bytes.len().saturating_sub(1);
    while common > 0 && char::from(bytes[common - 1]) != MAIN_SEPARATOR {
        common -= 1;
    }
    let (db_dir, entry) = if common > 1 {
        first.split_at(common)
    } else {
        (".", first.as_str())
    };

    let mut root = AliasNode::default();
    // The top node is named "-" (seqdbalias.cpp:106, 127).
    resolver.expand(&mut root, db_dir, "-", &[entry.to_string()], skip_local)?;
    volumes_in_order(root)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:85-97
// ```c++
// void CSeqDBAliasNode::x_Tokenize(const string & dbnames)
// {
//     vector<CSeqDB_Substring> dbs;
//     SeqDB_SplitQuoted(dbnames, dbs);
//
//     m_DBList.resize(dbs.size());
//     m_SkipLocal.resize(dbs.size(),false);
//
//     for(size_t i = 0; i<dbs.size(); i++) {
//         m_DBList[i].Assign(dbs[i]);
//         m_DBList[i].FixDelimiters();
//     }
// }
// ```
// This entry point is for a caller that holds the path of one database, not
// a list in the syntax above (`LOSAT makedb`, whose output name is one name
// by makeblastdb's own check): the name is not split, so a space or a double
// quote in it is a character of the path. NCBI has no such entry point: its
// reader splits every name it is given.
/// Resolves the one database `db_base` (a path without extension, taken as
/// it is): an alias file when `db_base.{p|n}al` exists, else the single
/// volume `db_base.{p|n}in`. `local_copy` says where the volumes named in an
/// alias file are looked for.
pub(super) fn resolve_one_name(
    db_base: &str,
    seqtype: SeqType,
    local_copy: LocalCopy<'_>,
) -> Result<AliasTree, SeqDbError> {
    let mut resolver = Resolver {
        seqtype,
        recurse: Vec::new(),
        local_dir: match local_copy {
            LocalCopy::CurrentDir => std::env::current_dir()
                .ok()
                .map(|cwd| cwd.to_string_lossy().into_owned()),
            LocalCopy::Dir(dir) => Some(dir.to_string()),
            LocalCopy::No => None,
        },
        search: None,
    };
    let db_base = fix_delimiters(db_base);
    if !exists(&resolver.alias_path(&db_base)) && !exists(&resolver.index_path(&db_base)) {
        return Err(SeqDbError::NotFound(format!(
            "No alias or index file found for {} database [{db_base}].",
            seqtype.name()
        )));
    }
    let mut root = AliasNode::default();
    // The top node is named "-" (seqdbalias.cpp:106, 127).
    resolver.expand(&mut root, "", "-", std::slice::from_ref(&db_base), true)?;
    volumes_in_order(root)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:970-984
// ```c++
// void CSeqDBAliasNode::x_FindVolumePaths(set<string> & vols, set<string> & alias) const
// {
//     ITERATE(TVolNames, iter, m_VolNames) {
//         vols.insert(iter->GetBasePathS());
//     }
//
//     string alias_base = m_ThisName.GetPathS();
//     if (alias_base != "-") {
//         alias.insert(m_ThisName.GetPathS());
//     }
//
//     ITERATE(TSubNodeList, iter, m_SubNodes) {
//         (*iter)->x_FindVolumePaths(vols, alias);
//     }
// }
// ```
fn collect_alias_files(node: &AliasNode, out: &mut std::collections::BTreeSet<String>) {
    if let Some(name) = &node.this_name {
        out.insert(name.clone());
    }
    for subnode in &node.subnodes {
        collect_alias_files(subnode, out);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1517-1532
// ```c++
// void
// CSeqDBImpl::FindVolumePaths(const string   & dbname,
//                             char             prot_nucl,
//                             vector<string> & paths,
//                             vector<string> * alias_paths,
//                             bool             recursive,
//                             bool             expand_links)
// {
//     bool use_atlas_lock = true;
//     CSeqDBAtlasHolder AH(NULL, use_atlas_lock);
//     CSeqDBAtlas & atlas(AH.Get());
//
//     // This constructor handles its own locking.
//     CSeqDBAliasFile aliases(atlas, dbname, prot_nucl, expand_links);
//     aliases.FindVolumePaths(paths, alias_paths, recursive);
// }
// ```
// The only caller is the removal of an existing database before a build
// (`DeleteBlastDb`), and for that use LOSAT departs from NCBI on purpose in
// two points:
//
// * `db_base` is one name (`resolve_one_name`). NCBI splits it at spaces and
//   double quotes, fails to resolve the pieces of a path that has one, and
//   so leaves the existing database in place.
// * A volume named in an alias file is looked for in the alias file's own
//   directory only (an absolute name as it is). NCBI looks in the current
//   directory first ("local copy", seqdbalias.cpp:831-840, quoted at
//   `Resolver::expand`), so `makeblastdb -out sub/db` over an existing
//   `sub/db.pal` removes the files `./db.00.*` of the current directory when
//   they exist, although they are not the files of `sub/db`. The files to
//   remove are those the alias file sits beside.
/// The files of the one database `db_base` for its removal: its volumes
/// (paths without extension, in NCBI's order) and the paths of the alias
/// files that name them (`CSeqDB::FindVolumePaths(dbname, seqtype, paths,
/// &alias_paths)`, with the two differences above).
pub(crate) fn find_volume_paths(
    db_base: &str,
    seqtype: SeqType,
) -> Result<(Vec<String>, Vec<String>), SeqDbError> {
    let tree = resolve_one_name(db_base, seqtype, LocalCopy::No)?;
    let mut alias_files = std::collections::BTreeSet::new();
    collect_alias_files(&tree.root, &mut alias_files);
    Ok((tree.volumes, alias_files.into_iter().collect()))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1578-1592
// ```c++
//     TVarList::const_iterator value =
//         m_Values.find(walker->GetFileKey());
//
//     if (value != m_Values.end()) {
//         walker->AddString( (*value).second );
//         return;
//     }
//
//     ITERATE(TSubNodeList, node, m_SubNodes) {
//         (*node)->WalkNodes( walker, volset );
//     }
//
//     ITERATE(TVolNames, volname, m_VolNames) {
//         if (const CSeqDBVol * vptr = volset.GetVol(volname->GetBasePathS())) {
//             walker->Accumulate( *vptr );
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:480-506
// ```c++
// void SeqDB_JoinDelim(string & a, const string & b, const string & delim)
// {
//     if (b.empty()) {
//         return;
//     }
//
//     if (a.empty()) {
//         // a has no size - but might have capacity
//         s_SeqDB_QuickAssign(a, b);
//         return;
//     }
//     ...
//     a += delim;
//     a += b;
// }
// ```
// `CSeqDB_TitleWalker` (seqdbalias.cpp:1011-1026) passes a volume's title and
// an alias file's `TITLE` to `SeqDB_JoinDelim(m_Value, value, "; ")`.
/// Appends the title of `node` to `title`: its `TITLE` when it has one, else
/// the titles below it joined with "; ".
pub(super) fn walk_title(
    node: &AliasNode,
    volume_title: &dyn Fn(&str) -> Option<Vec<u8>>,
    title: &mut Vec<u8>,
) {
    if let Some(value) = &node.title {
        join_delim(title, value);
        return;
    }
    for subnode in &node.subnodes {
        walk_title(subnode, volume_title, title);
    }
    for volume in &node.vol_names {
        if let Some(value) = volume_title(volume) {
            join_delim(title, &value);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1575-1596
// ```c++
// CSeqDBAliasNode::WalkNodes(CSeqDB_AliasWalker * walker,
//                            const CSeqDBVolSet & volset) const
// {
//     TVarList::const_iterator value =
//         m_Values.find(walker->GetFileKey());
//
//     if (value != m_Values.end()) {
//         walker->AddString( (*value).second );
//         return;
//     }
//
//     ITERATE(TSubNodeList, node, m_SubNodes) {
//         (*node)->WalkNodes( walker, volset );
//     }
//
//     ITERATE(TVolNames, volname, m_VolNames) {
//         if (const CSeqDBVol * vptr = volset.GetVol(volname->GetBasePathS())) {
//             walker->Accumulate( *vptr );
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1189-1192
// ```c++
//     virtual void Accumulate(const CSeqDBVol & vol)
//     {
//         m_Value += vol.GetNumOIDs();
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1201-1204
// ```c++
//     virtual void AddString(const string & value)
//     {
//         m_Value += NStr::StringToInt8(value);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1290-1293
// ```c++
//     virtual void Accumulate(const CSeqDBVol & vol)
//     {
//         m_Value += vol.GetVolumeLength();
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1302-1305
// ```c++
//     virtual void AddString(const string & value)
//     {
//         m_Value += NStr::StringToUInt8(value);
//     }
// ```
// The number of sequences (`NSEQ`) and the total length (`LENGTH`) of a
// database are such walks: a node that has the key contributes its value and
// is not descended into; any other node contributes the walks of its sub-nodes
// and its own volumes. A value equal to the sum the walk gives below its node
// therefore changes nothing, and LOSAT reads the key only in that case.
/// The number of sequences and the total length of `node`, as the walks of
/// `NSEQ` and `LENGTH` (quoted above) compute them over the volumes below it
/// when no node has the keys. `volume_sums` gives both for a volume path.
/// Every alias node below `node` that has `NSEQ` or `LENGTH` is checked on the
/// way: the value must be plain decimal digits (1 to 19) and equal the sum of
/// that node. NCBI's parse of other values is not reproduced.
pub(super) fn check_sums(
    node: &AliasNode,
    volume_sums: &dyn Fn(&str) -> Option<(u64, u64)>,
) -> Result<(u64, u64), SeqDbError> {
    let mut nseq: u64 = 0;
    let mut length: u64 = 0;
    let overflow = || SeqDbError::Unsupported {
        path: node
            .this_name
            .as_deref()
            .map_or_else(Default::default, |name| Path::new(name).to_path_buf()),
        message: "the number of sequences or the total length of the databases of an alias \
                  file overflows"
            .to_string(),
    };
    for subnode in &node.subnodes {
        let (sub_nseq, sub_length) = check_sums(subnode, volume_sums)?;
        nseq = nseq.checked_add(sub_nseq).ok_or_else(overflow)?;
        length = length.checked_add(sub_length).ok_or_else(overflow)?;
    }
    for volume in &node.vol_names {
        if let Some((vol_nseq, vol_length)) = volume_sums(volume) {
            nseq = nseq.checked_add(vol_nseq).ok_or_else(overflow)?;
            length = length.checked_add(vol_length).ok_or_else(overflow)?;
        }
    }
    let path = node
        .this_name
        .as_deref()
        .map_or_else(Default::default, |name| Path::new(name).to_path_buf());
    for (key, value, sum, what) in [
        ("NSEQ", &node.nseq, nseq, "the number of sequences"),
        ("LENGTH", &node.length, length, "the total length"),
    ] {
        let Some(value) = value else { continue };
        let digits = (1..=19).contains(&value.len()) && value.iter().all(u8::is_ascii_digit);
        let parsed = if digits {
            std::str::from_utf8(value)
                .ok()
                .and_then(|text| text.parse::<u64>().ok())
        } else {
            None
        };
        let shown = String::from_utf8_lossy(value);
        match parsed {
            Some(parsed) if parsed == sum => {}
            Some(_) => {
                return Err(SeqDbError::Unsupported {
                    path,
                    message: format!(
                        "alias file key {key} ({shown}) differs from {what} of its databases \
                         ({sum}); LOSAT reads NSEQ and LENGTH only when they equal those sums"
                    ),
                })
            }
            None => {
                return Err(SeqDbError::Unsupported {
                    path,
                    message: format!(
                        "alias file key {key} ({shown}) is not a plain decimal number of 1 to \
                         19 digits; LOSAT reads NSEQ and LENGTH only when they equal {what} of \
                         their databases ({sum})"
                    ),
                })
            }
        }
    }
    Ok((nseq, length))
}

fn join_delim(a: &mut Vec<u8>, b: &[u8]) {
    if b.is_empty() {
        return;
    }
    if !a.is_empty() {
        a.extend_from_slice(b"; ");
    }
    a.extend_from_slice(b);
}

#[cfg(test)]
pub(super) fn test_read_values(bytes: &[u8]) -> Vec<(Vec<u8>, Vec<u8>)> {
    read_values(bytes)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:929-958
// ```c++
// void CSeqDBAliasNode::FindVolumePaths(vector<string> & vols, vector<string> * alias, bool recursive) const
// {
//
//     set<string> volset;
//     set<string> aliset;
//
//     if (recursive) {
//         x_FindVolumePaths(volset, aliset);
//     } else {
//         // No alias file list is populated in the non-recursive mode
//         ITERATE(vector<CSeqDB_BasePath>, path, m_VolNames) {
//             volset.insert(path->GetBasePathS());
//         }
//         ITERATE(TSubNodeList, iter, m_SubNodes) {
//             ITERATE(vector<CSeqDB_BasePath>, path, (*iter)->m_VolNames) {
//                 volset.insert(path->GetBasePathS());
//             }
//             ITERATE(TSubNodeList, sub, (*iter)->m_SubNodes) {
//                 volset.insert(((*sub)->m_ThisName).GetPathS());
//             }
//         }
//     }
//
//     vols.clear();
//     ITERATE(set<string>, iter, volset) {
//         vols.push_back(*iter);
//     }
//
//     // Sort to insure deterministic order.
//     sort(vols.begin(), vols.end(), SeqDB_CompareVolume);
// ```
// blastdbcmd prints this list under "Volumes:" (`FindVolumePaths(volumes,
// false)`, blastdbcmd.cpp:568): the volumes of the first two levels of the
// alias tree and the alias files of the third. The tests compare it with what
// NCBI printed.
/// What `blastdbcmd -info` lists under "Volumes:" for the database `dbname`.
#[cfg(test)]
pub(super) fn test_info_volumes(
    dbname: &str,
    seqtype: SeqType,
    search: &SearchPath,
) -> Result<Vec<String>, SeqDbError> {
    let tree = resolve_in(dbname, seqtype, search)?;
    let mut vols: Vec<String> = tree.root.vol_names.clone();
    for node in &tree.root.subnodes {
        vols.extend(node.vol_names.iter().cloned());
        vols.extend(node.subnodes.iter().filter_map(|sub| sub.this_name.clone()));
    }
    vols.sort_by(|a, b| {
        path_base_name(a)
            .cmp(path_base_name(b))
            .then_with(|| a.cmp(b))
    });
    vols.dedup();
    Ok(vols)
}

#[cfg(test)]
pub(super) fn test_split_quoted(dbname: &[u8]) -> Vec<&[u8]> {
    split_quoted(dbname)
}
