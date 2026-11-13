//! Alias files (`.pal` / `.nal`) and the list of volumes behind a database name.
//!
//! LOSAT reads the keys a multi-volume database of makeblastdb uses, `TITLE`
//! and `DBLIST`, the counts NCBI's reader takes from its walks of the alias
//! tree (`NSEQ`, `LENGTH`, `MIN_SEQ_LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`:
//! `walk_totals`) and the OID ranges `FIRST_OID` / `LAST_OID` (`compute_masks`,
//! `range_mask`). The keys that select sequences by a list, an OID mask or a
//! membership bit, and `MASKLIST`, are explicit errors (not supported yet);
//! unknown keys, `VDBLIST` and `MAX_SEQ_LENGTH` among them, are ignored, as in
//! NCBI.
//!
//! Two ways lead from a name to its volumes. `resolve_in` is NCBI's: the
//! name of a search program's `-db` is looked for in the directories of a
//! [`SearchPath`] (the current directory, `BLASTDB`, the registry value
//! `[BLAST] BLASTDB`), names inside alias files are looked for there too when
//! they are found nowhere else, and volumes are recorded under normalised
//! directory names. `resolve_one_name` takes the path of one database as it
//! is and looks nowhere else (`losat makeblastdb`).

use std::collections::VecDeque;
use std::ffi::OsStr;
use std::path::{Path, MAIN_SEPARATOR};

use super::file::DbFile;
use super::volume::with_extension;
use super::{SeqDbError, SeqType};
use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::{string_to_int8, string_to_uint, string_to_uint8};

/// One node of the alias tree: the database name given by the caller (the
/// top node) or one alias file.
#[derive(Debug, Default)]
pub(super) struct AliasNode {
    /// The key / value lines of the alias file (`m_Values`, `read_values`), a
    /// repeated key once with its last value; empty for the top node, whose
    /// only value is the `DBLIST` of the names it was given.
    pub values: Vec<(Vec<u8>, Vec<u8>)>,
    /// `FIRST_OID` / `LAST_OID` of the alias file as `ComputeMasks` reads them
    /// (`compute_masks`); `None` before it ran and for a node without them.
    pub range: Option<OidRange>,
    /// Alias files named in `DBLIST`, in order.
    pub subnodes: Vec<AliasNode>,
    /// Volumes named in `DBLIST` (paths without extension), in order.
    pub vol_names: Vec<String>,
    /// Path of the alias file of this node (`m_ThisName`); `None` for the
    /// top node, which NCBI names "-".
    pub this_name: Option<String>,
}

impl AliasNode {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1578-1583
    // ```c++
    //     TVarList::const_iterator value =
    //         m_Values.find(walker->GetFileKey());
    //
    //     if (value != m_Values.end()) {
    // ```
    /// The value of `key` in the alias file of this node (`m_Values.find(key)`).
    pub fn value(&self, key: &str) -> Option<&[u8]> {
        self.values
            .iter()
            .find(|(name, _)| name == key.as_bytes())
            .map(|(_, value)| value.as_slice())
    }
}

/// The OID range of one alias file: the values of `FIRST_OID` and `LAST_OID`
/// as `NStr::StringToUInt` converts them (`compute_masks`), before NCBI keeps
/// them in an `int`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct OidRange {
    /// `FIRST_OID`, 1-based; `None` without the key.
    pub first_oid: Option<u32>,
    /// `LAST_OID`, 1-based and inclusive; `None` without the key.
    pub last_oid: Option<u32>,
}

impl OidRange {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1975-1998
    // ```c++
    //             int first_oid = 0;
    //             int last_oid  = INT_MAX;
    //             bool has_range = false;
    //
    //             if (f_oid_iter != m_Values.end()) {
    //                 first_oid = NStr::StringToUInt(f_oid_iter->second);
    //
    //                 // Starts at one, adjust to zero-indexed.
    //                 if (first_oid)
    //                     first_oid--;
    //
    //                 has_range = true;
    //             }
    //
    //             if (l_oid_iter != m_Values.end()) {
    //                 // Zero indexing and post notation adjustments cancel.
    //                 last_oid = NStr::StringToUInt(l_oid_iter->second);
    //                 has_range = true;
    //             }
    //
    //             if (has_range) {
    //                 CRef<TMask> mask(new TMask(first_oid, last_oid));
    //                 m_NodeMasks.push_back(mask);
    //             }
    // ```
    // `NStr::StringToUInt` gives a `Uint4` that NCBI stores in an `int`: GCC wraps a value of
    // 2^31 or more to a negative number (defined since C++20, implementation-defined before).
    // `first_oid--` of `FIRST_OID 2147483648` decrements `INT_MIN`, which is undefined; NCBI
    // BLAST+ 2.17.0 gives `INT_MAX` (probe p4a_r2_ao_fo2147483648 of line RG D1a-2), the value of
    // `wrapping_sub(1)` here, which is also 2^31 - 1 without the wrap.
    /// `(begin, end)` as NCBI's `CSeqDB_AliasMask` holds them (`int`): `begin` the 0-based
    /// first OID (0 without `FIRST_OID`), `end` the `LAST_OID` (`INT_MAX` without it), in the
    /// OIDs of the whole volume set. A value of 2^31 or more is negative, as NCBI reads it.
    pub fn bounds(&self) -> (i32, i32) {
        let begin = match self.first_oid {
            None | Some(0) => 0,
            Some(first) => first.wrapping_sub(1) as i32,
        };
        let end = self.last_oid.map_or(i32::MAX, |last| last as i32);
        (begin, end)
    }
}

/// The resolved database: the alias tree and its volumes in NCBI's order.
#[derive(Debug)]
pub(super) struct AliasTree {
    pub root: AliasNode,
    /// Volume paths without extension, sorted as NCBI opens them.
    pub volumes: Vec<String>,
}

/// Alias keys that LOSAT recognises and does not support: they select
/// sequences by a list, an OID mask or a membership bit (alias_files.txt
/// "Type: filter"), or name masks LOSAT does not read (`MASKLIST`). They are
/// the keys NCBI's reader looks up (seqdbalias.cpp) but those LOSAT reads:
/// `TITLE`, `DBLIST`, the counts of the walks (`walk_totals`) and the OID
/// ranges (`compute_masks`). `MAX_SEQ_LENGTH`, which NCBI never applies
/// (seqdbimpl.cpp:972-976 takes the maximum from the volumes), and `VDBLIST`,
/// which only NCBI's VDB tools read (vdbalias.cpp; `CSeqDBAliasNode` looks it
/// up nowhere), are ignored as unknown keys are.
const UNSUPPORTED_KEYS: [&str; 8] = [
    "GILIST",
    "TILIST",
    "SEQIDLIST",
    "TAXIDLIST",
    "OIDLIST",
    "MEMB_BIT",
    "MASKLIST",
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

    /// The search path that holds the current directory only (`SeqDb::open`: `BLASTDB` and
    /// the registry are not consulted).
    pub(crate) fn current_dir() -> Result<Self, SeqDbError> {
        Self::new(&std::env::current_dir().unwrap_or_default(), None, None)
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
// `keep_quote` is false for a `DBLIST` and for the names the reader opens
// (seqdbalias.cpp:85-88), true for the names of the reports (`report_names`).
// With `keep_quote`, the closing quote stays in the name and NCBI's `i++` also
// skips the character after it (`"a"b` gives `"a"` alone).
fn split_quoted(dbname: &[u8], keep_quote: bool) -> Vec<&[u8]> {
    let mut dbs = Vec::new();
    let mut quoted = false;
    let mut begin = 0;
    let mut i = 0;
    while i < dbname.len() {
        let ch = dbname[i];
        if quoted {
            if ch == b'"' {
                if begin < i {
                    if keep_quote {
                        i += 1;
                    }
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
            begin = if keep_quote { i } else { i + 1 };
            quoted = true;
        }
        i += 1;
    }
    if begin < dbname.len() {
        dbs.push(&dbname[begin..]);
    }
    dbs
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:479-500
// ```c++
//     vector<CTempString> dbs;
//     SeqDB_SplitQuoted(blastdb_names, dbs, true);
//     retval.reserve(dbs.size());
//
//     ITERATE(vector<CTempString>, i, dbs) {
//         CAlignFormatUtil::SDbInfo info;
//         info.is_protein = is_protein;
//         bool success = false;
// 	// Unsafe OK as kDbName only used in this loop.
//         const string kDbName = NStr::TruncateSpaces_Unsafe(*i);
//         if (kDbName.empty())
//             continue;
//
//         success = s_FillDbInfoLocally(kDbName, info, dbfilt_algorithm);
// ```
// `TruncateSpaces` drops `isspace` bytes at both ends (`db_argument_name`).
/// The names of the `-db` value `dbname` as NCBI's reports take them, one per database of
/// the list, each of which they open alone (`CAlignFormatUtil::GetBlastDbInfo`): split at
/// spaces with the double quotes kept, without white space at the ends, empty names left
/// out.
pub(crate) fn report_names(dbname: &str) -> Vec<String> {
    split_quoted(dbname.as_bytes(), true)
        .into_iter()
        .map(|name| super::db_argument_name(&String::from_utf8_lossy(name)).to_string())
        .filter(|name| !name.is_empty())
        .collect()
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
        let mut dblist: &[u8] = &[];
        for (key, value) in &values {
            match key.as_slice() {
                b"DBLIST" => dblist = value,
                other => {
                    if let Some(key) = UNSUPPORTED_KEYS.iter().find(|k| k.as_bytes() == other) {
                        return Err(unsupported(format!(
                            "alias file key {key} is not supported; LOSAT reads TITLE, DBLIST, NSEQ, \
                             LENGTH, STATS_NSEQ, STATS_TOTLEN, MIN_SEQ_LENGTH, FIRST_OID and LAST_OID"
                        )));
                    }
                }
            }
        }
        let mut entries = Vec::new();
        for entry in split_quoted(dblist, false) {
            let entry = std::str::from_utf8(entry)
                .map_err(|_| unsupported("a DBLIST name is not valid UTF-8".to_string()))?;
            entries.push(fix_delimiters(entry));
        }

        self.recurse.push(path.to_string_lossy().into_owned());
        let mut node = AliasNode {
            this_name: Some(path.to_string_lossy().into_owned()),
            values,
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
// A volume reached from several nodes (an alias file and a name of `-db`, or
// two alias files) is one entry of the set, so it is opened once and its OIDs
// are numbered once; NCBI's walks still count it once per node that names it
// (`walk_totals`), and its OIDs are those of every path to it (`range_mask`).
// A node lists a volume once however often its `DBLIST` names it
// (`push_unique`, seqdbalias.cpp:789-798, 864-873, 905-914). Volumes of the
// same base name in different directories (two databases of one name, `-db
// "d1/g d2/g"`) are opened in the order above, each with the LMDB file of its
// own directory (`lmdb_entries`), as NCBI does.
/// The volumes of the alias tree `root` in NCBI's order.
fn volumes_in_order(root: AliasNode) -> AliasTree {
    let mut volumes = Vec::new();
    collect_volumes(&root, &mut volumes);
    volumes.sort_by(|a, b| {
        path_base_name(a)
            .cmp(path_base_name(b))
            .then_with(|| a.cmp(b))
    });
    volumes.dedup();
    AliasTree { root, volumes }
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
//         if (m_DBList[i].GetBasePathS().size() < common) {
//             common = m_DBList[i].GetBasePathS().size();
//         }
//     }
//
//     if (common) {
//         --common;
//     }
//
//     // Reduce common length to largest universal prefix.
//     const string & first = m_DBList[0].GetBasePathS();
//
//     for(i = 1; common && (i < m_DBList.size()); i++) {
//         // Reduce common prefix length until match is found.
//         while(memcmp(first.c_str(),
//                      m_DBList[i].GetBasePathS().c_str(),
//                      common)) {
//             --common;
//         }
//     }
//
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
// Several names are the databases of one search, the entries of the top
// alias node "-" as the names of a `DBLIST` are; a volume that two names give
// is listed once by the top node (`push_unique`), one that a name and an
// alias file give is opened once and counted by both (`volumes_in_order`). A
// name carries no extension: `db.pin` is looked for as `db.pin.pal` and
// `db.pin.pin`. A name with a path separator is joined to every
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
    let names = split_quoted(dbname.as_bytes(), false);
    // `x_ExpandAliases` on an empty list (quoted at `Resolver::expand`).
    if names.is_empty() {
        return Err(SeqDbError::NotFound(
            "No database names were passed to CSeqDB::CSeqDB().".to_string(),
        ));
    }
    let mut resolver = Resolver {
        seqtype,
        recurse: Vec::new(),
        local_dir: Some(search.cwd.clone()),
        search: Some(search),
    };

    // Every name is resolved before any is expanded; the first one not found
    // is the error.
    let mut resolved_names = Vec::with_capacity(names.len());
    let mut skip_local = Vec::with_capacity(names.len());
    for name in names {
        let name = fix_delimiters(&String::from_utf8_lossy(name));
        skip_local.push(name.contains(MAIN_SEPARATOR));

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
        resolved_names.push(format!(
            "{}{MAIN_SEPARATOR}{}",
            normalize_path(dir_name(&resolved_path))?,
            remove_extn(base_name(&resolved_path))
        ));
    }

    // The common prefix of the resolved names, cut back to a whole directory
    // (a separator included), is the directory of the top node.
    let first = resolved_names[0].as_bytes();
    let mut common = resolved_names.iter().map(String::len).min().unwrap_or(0);
    common = common.saturating_sub(1);
    for other in &resolved_names[1..] {
        while common > 0 && first[..common] != other.as_bytes()[..common] {
            common -= 1;
        }
    }
    while common > 0 && char::from(first[common - 1]) != MAIN_SEPARATOR {
        common -= 1;
    }
    let db_dir = if common > 1 {
        &resolved_names[0][..common]
    } else {
        "."
    };

    let mut root = AliasNode::default();
    // The top node is named "-" (seqdbalias.cpp:106, 127). One `expand` per
    // name: each has its own `m_SkipLocal`, and `expand` keeps no state from
    // one entry to the next but the node's volume list.
    for (name, skip_local) in resolved_names.iter().zip(skip_local) {
        let entry = if common > 1 {
            &name[common..]
        } else {
            name.as_str()
        };
        resolver.expand(&mut root, db_dir, "-", &[entry.to_string()], skip_local)?;
    }
    Ok(volumes_in_order(root))
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
// a list in the syntax above (`losat makeblastdb`, whose output name is one name
// by makeblastdb's own check): the name is not split, so a space or a double
// quote in it is a character of the path. NCBI has no such entry point: its
// reader splits every name it is given.
/// Resolves the one database `db_base` (a path without extension, taken as
/// it is): an alias file when `db_base.{p|n}al` exists, else the single
/// volume `db_base.{p|n}in`. `local_copy` says where the volumes named in an
/// alias file are looked for; `search` is the search path that a database which is not there
/// names in NCBI's message (the files are looked for at `db_base` only).
pub(super) fn resolve_one_name(
    db_base: &str,
    seqtype: SeqType,
    local_copy: LocalCopy<'_>,
    search: &SearchPath,
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
            "No alias or index file found for {} database [{db_base}] in search path [{}]. \
             Please verify the spelling of the BLAST database and its molecule type.",
            seqtype.name(),
            search.as_str()
        )));
    }
    let mut root = AliasNode::default();
    // The top node is named "-" (seqdbalias.cpp:106, 127).
    resolver.expand(&mut root, "", "-", std::slice::from_ref(&db_base), true)?;
    Ok(volumes_in_order(root))
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
pub(super) fn collect_alias_files(node: &AliasNode, out: &mut std::collections::BTreeSet<String>) {
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
    let tree = resolve_one_name(db_base, seqtype, LocalCopy::No, &SearchPath::current_dir()?)?;
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
    if let Some(value) = node.value("TITLE") {
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
// Every volume a node names is in the volume set (`volumes_in_order` takes them all), so
// `GetVol` finds it; a volume that several nodes name is accumulated once per node.
/// What NCBI's `WalkNodes` hands its walker: the value of a node that has the walker's key
/// (`AddString`; the walk does not go below that node) or a volume of a node without it
/// (`Accumulate`).
enum Reached<'a> {
    Value(&'a [u8]),
    Volume(&'a str),
}

/// NCBI's `WalkNodes` with a walker whose key is `key`: `visit` takes what the walk reaches in
/// NCBI's order. The first exception of `visit` ends the walk, as the exception of `AddString`
/// leaves NCBI's.
fn walk_nodes<'a>(
    node: &'a AliasNode,
    key: &str,
    visit: &mut dyn FnMut(Reached<'a>) -> Result<(), NcbiException>,
) -> Result<(), NcbiException> {
    if let Some(value) = node.value(key) {
        return visit(Reached::Value(value));
    }
    for subnode in &node.subnodes {
        walk_nodes(subnode, key, visit)?;
    }
    for volume in &node.vol_names {
        visit(Reached::Volume(volume))?;
    }
    Ok(())
}

/// The keys of NCBI's count walks (`walk_totals`). A value that does not convert fails the
/// walk of its key, where NCBI's reader runs that walk (`SeqDb::open_tree`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum WalkKey {
    Nseq,
    Length,
    MinSeqLength,
    StatsNseq,
    StatsTotlen,
}

/// What NCBI's walks of the alias tree give (`walk_totals`), each with the exception of the
/// first value of its key that does not convert.
#[derive(Debug)]
pub(super) struct AliasTotals {
    /// `CSeqDBAliasNode::GetNumSeqs`: `NSEQ`, else the volumes' OIDs (`int` of it is
    /// `GetNumSeqs`).
    pub nseq: Result<i64, NcbiException>,
    /// `CSeqDBAliasNode::GetTotalLength`: `LENGTH`, else the volumes' lengths.
    pub length: Result<u64, NcbiException>,
    /// `CSeqDBAliasNode::GetMinLength`: the smallest `MIN_SEQ_LENGTH`, 0 for a branch of
    /// volumes alone.
    pub min_seq_length: Result<i32, NcbiException>,
    /// `CSeqDBAliasNode::GetNumSeqsStats`: `STATS_NSEQ`, 0 for a branch of volumes alone.
    pub stats_nseq: Result<i64, NcbiException>,
    /// `CSeqDBAliasNode::GetTotalLengthStats`: `STATS_TOTLEN`, 0 for a branch of volumes
    /// alone.
    pub stats_totlen: Result<u64, NcbiException>,
}

impl Default for AliasTotals {
    /// The walks of a tree without the keys and without volumes (a placeholder until
    /// `walk_totals` runs).
    fn default() -> Self {
        Self {
            nseq: Ok(0),
            length: Ok(0),
            min_seq_length: Ok(0),
            stats_nseq: Ok(0),
            stats_totlen: Ok(0),
        }
    }
}

impl AliasTotals {
    /// The exception of the walk of `key`, if its walk failed.
    pub fn error(&self, key: WalkKey) -> Option<&NcbiException> {
        match key {
            WalkKey::Nseq => self.nseq.as_ref().err(),
            WalkKey::Length => self.length.as_ref().err(),
            WalkKey::MinSeqLength => self.min_seq_length.as_ref().err(),
            WalkKey::StatsNseq => self.stats_nseq.as_ref().err(),
            WalkKey::StatsTotlen => self.stats_totlen.as_ref().err(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1189-1204
// ```c++
//     virtual void Accumulate(const CSeqDBVol & vol)
//     {
//         m_Value += vol.GetNumOIDs();
//     }
//     ...
//     virtual void AddString(const string & value)
//     {
//         m_Value += NStr::StringToInt8(value);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1257-1260
// ```c++
//     virtual void Accumulate(const CSeqDBVol & vol)
//     {
//         // only alias file data is included.
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1290-1305
// ```c++
//     virtual void Accumulate(const CSeqDBVol & vol)
//     {
//         m_Value += vol.GetVolumeLength();
//     }
//     ...
//     virtual void AddString(const string & value)
//     {
//         m_Value += NStr::StringToUInt8(value);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1107-1147
// ```c++
//     CSeqDB_MinLengthWalker()
//     {
//         m_Value = INT4_MAX;
//     }
//     ...
//     virtual void Accumulate(const CSeqDBVol & vol)
//     {
//         int new_min = vol.GetMinLength();
//
//         if (new_min < m_Value) {
//             m_Value = new_min;
//         }
//     }
//     ...
//     virtual void AddString(const string & value)
//     {
//         int new_min = NStr::StringToUInt(value);
//
//         if (new_min < m_Value) {
//             m_Value = new_min;
//         }
//     }
// ```
// `CSeqDB_NSeqsStatsWalker` and `CSeqDB_TotalLengthStatsWalker` are the `NSEQ` and `LENGTH`
// walkers with the keys `STATS_NSEQ` and `STATS_TOTLEN` and the empty `Accumulate` quoted
// above (seqdbalias.cpp:1241-1262, 1342-1363). A volume's `GetMinLength` is its index file's
// `m_MinLen`, which no index file sets (seqdbfile.cpp:142 starts it at 0; the formats store
// none), so a branch of volumes alone makes the minimum 0. NCBI adds to an `Int8` with
// signed overflow undefined; only the `int` of the sum is ever read (`x_GetNumSeqs`,
// `GetNumSeqsStats`), whose 32 bits a wrapping sum gives as any wider sum would. The
// unsigned sums of the lengths wrap as NCBI's `Uint8` does.
/// NCBI's walks of the counts of the alias tree `root`; `volume` gives the number of OIDs and
/// the total length of a volume path (`GetNumOIDs`, `GetVolumeLength`). The values parse
/// with NCBI's `NStr` functions and flags (`NSEQ`, `STATS_NSEQ`: `StringToInt8`; `LENGTH`,
/// `STATS_TOTLEN`: `StringToUInt8`; `MIN_SEQ_LENGTH`: `StringToUInt`, read as an `int`).
pub(super) fn walk_totals(root: &AliasNode, volume: &dyn Fn(&str) -> (u32, u64)) -> AliasTotals {
    let mut nseq: i64 = 0;
    let nseq = walk_nodes(root, "NSEQ", &mut |reached| {
        let add = match reached {
            Reached::Value(value) => string_to_int8(value, 0)?,
            Reached::Volume(path) => i64::from(volume(path).0),
        };
        nseq = nseq.wrapping_add(add);
        Ok(())
    })
    .map(|()| nseq);
    let mut length: u64 = 0;
    let length = walk_nodes(root, "LENGTH", &mut |reached| {
        let add = match reached {
            Reached::Value(value) => string_to_uint8(value)?,
            Reached::Volume(path) => volume(path).1,
        };
        length = length.wrapping_add(add);
        Ok(())
    })
    .map(|()| length);
    let mut min = i32::MAX;
    let min_seq_length = walk_nodes(root, "MIN_SEQ_LENGTH", &mut |reached| {
        // `int new_min = NStr::StringToUInt(value)`: a value of 2^31 or more is a negative `int`.
        let new_min = match reached {
            Reached::Value(value) => string_to_uint(value)? as i32,
            Reached::Volume(_) => 0,
        };
        min = min.min(new_min);
        Ok(())
    })
    .map(|()| min);
    let mut stats_nseq: i64 = 0;
    let stats_nseq = walk_nodes(root, "STATS_NSEQ", &mut |reached| {
        if let Reached::Value(value) = reached {
            stats_nseq = stats_nseq.wrapping_add(string_to_int8(value, 0)?);
        }
        Ok(())
    })
    .map(|()| stats_nseq);
    let mut stats_totlen: u64 = 0;
    let stats_totlen = walk_nodes(root, "STATS_TOTLEN", &mut |reached| {
        if let Reached::Value(value) = reached {
            stats_totlen = stats_totlen.wrapping_add(string_to_uint8(value)?);
        }
        Ok(())
    })
    .map(|()| stats_totlen);
    AliasTotals {
        nseq,
        length,
        min_seq_length,
        stats_nseq,
        stats_totlen,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1963-1980
// ```c++
//     if (! m_DBList.empty()) {
//         if (oid_iter   != m_Values.end() ||
//             gil_iter   != m_Values.end() ||
//             til_iter   != m_Values.end() ||
//             sil_iter   != m_Values.end() ||
//             f_oid_iter != m_Values.end() ||
//             l_oid_iter != m_Values.end() ||
//             mbit_iter  != m_Values.end() ||
//             taxid_iter != m_Values.end()) {
//
//             has_filters = true;
//
//             int first_oid = 0;
//             int last_oid  = INT_MAX;
//             bool has_range = false;
//
//             if (f_oid_iter != m_Values.end()) {
//                 first_oid = NStr::StringToUInt(f_oid_iter->second);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:2094-2097
// ```c++
//     NON_CONST_ITERATE(TSubNodeList, sn, m_SubNodes) {
//         (**sn).ComputeMasks(has_filters);
//     }
// }
// ```
// The rest of the range is quoted at `OidRange::bounds`. NCBI's constructor runs this first
// (`m_Aliases.HasFilters()`, seqdbimpl.cpp:96), before any count walk. The other filter keys
// are refused when the alias file is read, so a node's filters are its range alone.
/// `ComputeMasks` of the top node: `FIRST_OID`, then `LAST_OID`, of every node from the top
/// down, a node before its sub-nodes (`NStr::StringToUInt`); the first value that does not
/// convert is the error.
pub(super) fn compute_masks(node: &mut AliasNode) -> Result<(), NcbiException> {
    let first_oid = node.value("FIRST_OID").map(string_to_uint).transpose()?;
    let last_oid = node.value("LAST_OID").map(string_to_uint).transpose()?;
    if first_oid.is_some() || last_oid.is_some() {
        node.range = Some(OidRange {
            first_oid,
            last_oid,
        });
    }
    for subnode in &mut node.subnodes {
        compute_masks(subnode)?;
    }
    Ok(())
}

/// Whether some node of the tree has an OID range (`CSeqDB_FilterTree::HasFilter`, the filters
/// being ranges alone).
pub(super) fn has_ranges(node: &AliasNode) -> bool {
    node.range.is_some() || node.subnodes.iter().any(has_ranges)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:91-99
// ```c++
//     for(int i = 0; i < volset.GetNumVols(); i++) {
//         const CSeqDBVolEntry * v1 = volset.GetVolEntry(i);
//
//         CRef<CSeqDB_BitSet> vol_bits =
//             x_ComputeFilters(filters, *v1, gi_list_set, locked, lmdb_set.IsBlastDBVersion5());
//
//         m_AllBits->UnionWith(*vol_bits, true);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:131-177
// ```c++
//     const string & vn = vol.Vol()->GetVolName();
//     CRef<CSeqDB_FilterTree> ft = filters.Specialize(vn);
//
//     int vol_start = vol.OIDStart();
//     int vol_end   = vol.OIDEnd();
//     ...
//     if (vols > 0) {
//         ...
//         volume_map.Reset(new CSeqDB_BitSet(vol_start,
//                                      vol_end,
//                                      CSeqDB_BitSet::eAllSet));
//     } else {
//         ...
//         ITERATE(vector< CRef< CSeqDB_FilterTree > >, sub, ft->GetNodes()) {
//             CRef<CSeqDB_BitSet> sub_bits =
//                 x_ComputeFilters(**sub, vol, gis, locked, isBlastDBv5);
//
//             volume_map->UnionWith(*sub_bits, true);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:193-206
// ```c++
//     CRef<CSeqDB_BitSet> filter(new CSeqDB_BitSet(vol_start,
//                                                  vol_end,
//                                                  CSeqDB_BitSet::eAllSet));
//
//     // First, apply any 'range' filters, because they can be combined
//     // very efficiently.
//
//     typedef CSeqDB_FilterTree::TFilters TFilters;
//
//     ITERATE(TFilters, range, ft->GetFilters()) {
//         const CSeqDB_AliasMask & mask = **range;
//
//         if (mask.GetType() == CSeqDB_AliasMask::eOidRange) {
//             CSeqDB_BitSet r2(mask.GetBegin(),
//                                 mask.GetEnd(),
//                                 CSeqDB_BitSet::eAllSet);
//             filter->IntersectWith(r2, true);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfilter.cpp:38-55
// ```c++
// CRef<CSeqDB_FilterTree> CSeqDB_FilterTree::Specialize(string volname) const
// {
//     CRef<CSeqDB_FilterTree> clone(new CSeqDB_FilterTree);
//     clone->SetName(m_Name);
//     clone->AddFilters(m_Filters);
//
//     ITERATE(vector< CRef<CSeqDB_FilterTree> >, iter, m_SubNodes) {
//         CRef<CSeqDB_FilterTree> sub = (**iter).Specialize(volname);
//
//         if (sub.NotEmpty()) {
//             if (! sub->GetFilters().empty()) {
//                 clone->AddNode(sub);
//             } else {
//                 clone->AddNodes(sub->GetNodes());
//                 clone->AddVolumes(sub->GetVolumes());
//             }
//         }
//     }
// ```
// The tree of one volume keeps the branches that reach it and lifts a node without a filter
// into its parent. A node's OIDs of the volume are its range (the whole volume without one)
// intersected with the volume when the node names the volume, else with the union of its
// sub-nodes' OIDs; a sub-node's OIDs below a node that names the volume are inside the whole
// volume, so ignoring them changes nothing. By distribution the OIDs of a volume are therefore
// the union, over the paths from the top node to each node that names the volume, of the
// intersection of the ranges on the path with the volume, which is what is computed here.
// The ranges count the OIDs of the whole volume set (`vol_start`, `vol_end` are global).
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.hpp:104-114
// ```c++
//     CSeqDB_BitSet(size_t start, size_t end, ESpecialCase sp = eNone)
//         : m_Start  (start),
//           m_End    (end),
//           m_Special(sp)
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:208-221
// ```c++
// void CSeqDB_BitSet::IntersectWith(const CSeqDB_BitSet & other, bool consume)
// {
//     if (m_Special == eAllClear || other.m_Special == eAllSet) {
//         return;
//     }
//     ...
//     if (m_Special == eAllSet && other.m_Special == eAllSet) {
//         m_Start = max(m_Start, other.m_Start);
//         m_End = min(m_End, other.m_End);
//         return;
//     }
// ```
// The `int` bounds of a range become `size_t`: a negative one is sign-extended past every OID
// (line RG D1a-2, survey A M2), so a negative `FIRST_OID` keeps no OID and a negative
// `LAST_OID` is no upper bound, as `bound` reads them.
/// The OIDs that the OID ranges of the tree include, one bit per OID of the volume set (bit
/// `oid % 64` of word `oid / 64`); `volume` gives the global OIDs `[start, end)` of a volume
/// path. `None` when no node has a range: NCBI builds no OID list then, and every OID is
/// included.
pub(super) fn range_mask(
    root: &AliasNode,
    volume: &dyn Fn(&str) -> (u32, u32),
    num_oids: u32,
) -> Option<Vec<u64>> {
    if !has_ranges(root) {
        return None;
    }
    /// NCBI's `size_t` of an `int` bound.
    fn bound(value: i32) -> u64 {
        i64::from(value) as u64
    }
    fn mark(
        node: &AliasNode,
        lo: u64,
        hi: u64,
        volume: &dyn Fn(&str) -> (u32, u32),
        words: &mut [u64],
    ) {
        let (lo, hi) = match node.range {
            Some(range) => {
                let (begin, end) = range.bounds();
                (lo.max(bound(begin)), hi.min(bound(end)))
            }
            None => (lo, hi),
        };
        for path in &node.vol_names {
            let (start, end) = volume(path);
            set_bits(words, lo.max(u64::from(start)), hi.min(u64::from(end)));
        }
        for subnode in &node.subnodes {
            mark(subnode, lo, hi, volume, words);
        }
    }
    let mut words = vec![0u64; (num_oids as usize).div_ceil(64)];
    mark(root, 0, u64::MAX, volume, &mut words);
    Some(words)
}

/// Sets the bits `[from, to)` of `words` (nothing when `to <= from`).
fn set_bits(words: &mut [u64], from: u64, to: u64) {
    if to <= from {
        return;
    }
    let mut at = from;
    while at < to {
        let word = (at / 64) as usize;
        let bit = at % 64;
        let span = (64 - bit).min(to - at);
        let ones = if span == 64 {
            u64::MAX
        } else {
            ((1u64 << span) - 1) << bit
        };
        words[word] |= ones;
        at += span;
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:247-258
// ```c++
//         // Handle the case where FIRST_OID and LAST_OID is set on a top level
//         // alias file and that's the only alias file present
//         if (ft->HasFilter()) {
//             const vector< CRef<CSeqDB_FilterTree> >& nodes = ft->GetNodes();
//             if (nodes.size() == 1) {
//                 const CSeqDB_FilterTree::TFilters& filters = nodes.front()->GetFilters();
//                 if (filters.size() == 1 && filters.front()->GetType() == CSeqDB_AliasMask::eOidRange) {
//                     const CSeqDB_AliasMask& alias_mask = *filters.front();
//                     SetIterationRange(alias_mask.GetBegin(), alias_mask.GetEnd());
//                 }
//             }
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:183-201
// ```c++
// void CSeqDBImpl::SetIterationRange(int oid_begin, int oid_end)
// {
//     ...
//     m_RestrictBegin = (oid_begin < 0) ? 0 : oid_begin;
//     m_RestrictEnd   = (oid_end   < 0) ? 0 : oid_end;
//
//     if ((oid_begin == 0) && (oid_end == 0)) {
//         m_RestrictEnd = m_VolSet.GetNumOIDs();
//     } else {
//         if ((oid_end == 0) || (m_RestrictEnd > m_VolSet.GetNumOIDs())) {
//             m_RestrictEnd = m_VolSet.GetNumOIDs();
//         }
//         if (m_RestrictBegin > m_RestrictEnd) {
//             m_RestrictBegin = m_RestrictEnd;
//         }
//     }
// }
// ```
// The filter tree's top node has one sub-node per alias file of the top node
// (`BuildFilterTree`, seqdbalias.cpp:1880-1894), and a node's only filter is its range (the
// other filter keys are refused). The bounds are NCBI's `int`s (`OidRange::bounds`), negative
// for a value of 2^31 or more: a negative `oid_begin` starts at 0, a negative `oid_end` ends at
// 0, and `LAST_OID 0` (`oid_end` 0) ends at the last OID.
/// The iteration range NCBI sets on its first OID-list build when the top node has exactly one
/// alias file and that file has an OID range: `[begin, end)` as `SetIterationRange` clamps it
/// to the `num_oids` OIDs of the volume set. `None` in every other case.
pub(super) fn top_iteration_range(root: &AliasNode, num_oids: u32) -> Option<(u32, u32)> {
    let [only] = root.subnodes.as_slice() else {
        return None;
    };
    let (oid_begin, oid_end) = only.range?.bounds();
    let num_oids = i64::from(num_oids);
    let mut restrict_begin = i64::from(oid_begin.max(0));
    let mut restrict_end = i64::from(oid_end.max(0));
    if oid_begin == 0 && oid_end == 0 {
        restrict_end = num_oids;
    } else {
        if oid_end == 0 || restrict_end > num_oids {
            restrict_end = num_oids;
        }
        if restrict_begin > restrict_end {
            restrict_begin = restrict_end;
        }
    }
    Some((restrict_begin as u32, restrict_end as u32))
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
pub(super) fn info_volumes(
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
pub(super) fn test_split_quoted(dbname: &[u8], keep_quote: bool) -> Vec<&[u8]> {
    split_quoted(dbname, keep_quote)
}
