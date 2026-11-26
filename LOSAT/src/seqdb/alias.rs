//! Alias files (`.pal` / `.nal`) and the list of volumes behind a database name.
//!
//! LOSAT reads the keys a multi-volume database of makeblastdb uses, `TITLE`
//! and `DBLIST`, the counts NCBI's reader takes from its walks of the alias
//! tree (`NSEQ`, `LENGTH`, `MIN_SEQ_LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN`:
//! `walk_totals`), the OID ranges `FIRST_OID` / `LAST_OID`, the OID lists of
//! `OIDLIST` and the mask type of each, `OID_MASK_TYPE` (`compute_masks`,
//! `alias_oids`), the id lists `GILIST`, `TILIST`, `SEQIDLIST` and `TAXIDLIST`
//! (`compute_masks`), the membership bit `MEMB_BIT` (`compute_masks`,
//! `alias_oids`), and `MASKLIST`, the names of the GI-indexed masks of the top
//! alias file (`gi_mask_list`, read by `super::gi_mask`). Unknown keys,
//! `VDBLIST` and `MAX_SEQ_LENGTH` among them, are ignored, as in NCBI.
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
use std::path::{Path, PathBuf, MAIN_SEPARATOR};

use crate::utils::os_bytes::{os_str_bytes, path_bytes, path_from_bytes, path_of_bytes};

use super::file::DbFile;
use super::id_list::{AliasIdList, IdListKind, PackedStrings, SeqIdSlice};
use super::oid_list::{self, OidListError, OidListFile};
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
    /// The path of the OID list file of the alias file's `OIDLIST` as `ComputeMasks`
    /// resolves it (`compute_masks`); `None` before it ran and for a node without the key.
    pub oid_list: Option<Vec<u8>>,
    /// The `OID_MASK_TYPE` of that OID list (`EOidMaskType` bits, `NStr::StringToUInt` of the
    /// value: `compute_masks`); 0 without the key, before `compute_masks` ran and for a node
    /// without `OIDLIST`.
    pub oid_mask_type: u32,
    /// The id lists of the alias file's `GILIST`, `TILIST`, `SEQIDLIST` and `TAXIDLIST`, in that
    /// order, with their paths as `ComputeMasks` resolves them (`compute_masks`); empty before it
    /// ran.
    pub lists: Vec<AliasList>,
    /// `MEMB_BIT` of the alias file as `ComputeMasks` converts it (`NStr::StringToUInt`,
    /// `compute_masks`); `None` before it ran and for a node without the key. `Some(0)` stays a
    /// filter of the node, as NCBI pushes an alias mask of 0.
    pub memb_bit: Option<u32>,
    /// The directory of the alias file (`m_DBPath`), which a relative `OIDLIST` or id list is
    /// joined to; empty for the top node.
    pub db_path: Vec<u8>,
    /// Alias files named in `DBLIST`, in order.
    pub subnodes: Vec<AliasNode>,
    /// Volumes named in `DBLIST` (paths without extension), in order.
    pub vol_names: Vec<Vec<u8>>,
    /// Path of the alias file of this node (`m_ThisName`); `None` for the
    /// top node, which NCBI names "-".
    pub this_name: Option<Vec<u8>>,
}

impl AliasNode {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1578-1583
    /// The value of `key` in the alias file of this node (`find` of `key` in `m_Values`).
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
    // `NStr::StringToUInt` gives a `Uint4` that NCBI stores in an `int`: GCC wraps a value of 2^31
    // or more to a negative number (defined since C++20, implementation-defined before). For
    // `FIRST_OID 2147483648`, `first_oid` holds `INT_MIN` when NCBI decrements it, which is
    // undefined; NCBI BLAST+ 2.17.0 gives `INT_MAX` (probe p4a_r2_ao_fo2147483648 of line RG
    // D1a-2), the value of `wrapping_sub(1)` here, which is also 2^31 - 1 without the wrap.
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

/// An id list of an alias file (`CSeqDB_AliasMask` of `eGiList`, `eTiList`, `eSiList` or
/// `eTaxIdList`): its kind and the path of its file, the value joined to the alias file's
/// directory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct AliasList {
    pub kind: IdListKind,
    pub path: Vec<u8>,
}

/// The resolved database: the alias tree and its volumes in NCBI's order.
#[derive(Debug)]
pub(super) struct AliasTree {
    pub root: AliasNode,
    /// Volume paths without extension, sorted as NCBI opens them.
    pub volumes: Vec<Vec<u8>>,
    /// The names of the GI-indexed masks that NCBI's reader uses for the database
    /// (`gi_mask_list`), their bytes; `None` where it uses the mask columns of the volumes.
    pub gi_masks: Option<Vec<Vec<u8>>>,
}

/// The alias keys whose value is a path (`compute_masks`).
const PATH_KEYS: [&str; 5] = ["OIDLIST", "GILIST", "TILIST", "SEQIDLIST", "TAXIDLIST"];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:232-262
/// Extensions that `s_SeqDB_DBExists` (`SearchPath::find`) checks for a database name, in NCBI's
/// lookup order: alias file (`.nal`/`.pal`) before index file (`.nin`/`.pin`).
const DB_LOOKUP_EXTNS: [&str; 2] = ["al", "in"];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:48
/// Group alias file name checked before an alias file or volume (`kSeqDBGroupAliasFileName`).
const GROUP_ALIAS_FILE: &[u8] = b"index.alx";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:160-182
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3204-3216
/// `CSeqDBAtlas::DoesFileExist`: `path` is a regular file, or a symbolic link
/// to one.
pub(super) fn exists(path: &Path) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_file())
}

/// The path separator of the host (`CDirEntry::GetPathSeparator`) as a byte: NCBI's paths are
/// the bytes of a `std::string`.
const SEPARATOR: u8 = MAIN_SEPARATOR as u8;

fn dir_name(path: &[u8]) -> &[u8] {
    path.iter()
        .rposition(|&byte| byte == SEPARATOR)
        .map_or(&[], |at| &path[..at])
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:50-59
/// `SeqDB_RemoveDirName` (`CSeqDB_Path::FindFileName`): `path` after its last path separator.
pub(crate) fn base_name(path: &[u8]) -> &[u8] {
    path.iter()
        .rposition(|&byte| byte == SEPARATOR)
        .map_or(path, |at| &path[at + 1..])
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:88-107
// `CSeqDB_Path::FindBaseName` is `SeqDB_RemoveExtn` of `SeqDB_RemoveDirName` of the path
// (seqdbgeneral.hpp:975-979): the file name without a final `.pal`, `.nal`,
// `.pin`, `.nin`, `.pdb` or `.ndb`. `SeqDB_CompareVolume` compares these.
pub(crate) fn path_base_name(path: &[u8]) -> &[u8] {
    remove_extn(base_name(path))
}

/// `SeqDB_RemoveExtn` (cited above): `s` without a final `.pal`, `.nal`,
/// `.pin`, `.nin`, `.pdb` or `.ndb`, when something precedes it.
fn remove_extn(s: &[u8]) -> &[u8] {
    let len = s.len();
    if len > 4
        && s[len - 4] == b'.'
        && (s[len - 3] == b'n' || s[len - 3] == b'p')
        && matches!(&s[len - 2..], b"al" | b"in" | b"db")
    {
        &s[..len - 4]
    } else {
        s
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:284-295
pub(super) fn fix_delimiters(name: &[u8]) -> Vec<u8> {
    name.iter()
        .map(|&byte| {
            if byte == b'/' || byte == b'\\' {
                SEPARATOR
            } else {
                byte
            }
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:140-181
pub(super) fn combine_path(one: &[u8], two: &[u8]) -> Vec<u8> {
    if two.is_empty() {
        return one.to_vec();
    }
    let mut only_two = one.is_empty() || two.first() == Some(&SEPARATOR);
    if SEPARATOR == b'\\'
        && two.len() > 3
        && two[0].is_ascii_alphabetic()
        && two[1] == b':'
        && two[2] == b'\\'
    {
        only_two = true;
    }
    if only_two {
        return two.to_vec();
    }
    let mut out = one.to_vec();
    if out.last() != Some(&SEPARATOR) {
        out.push(SEPARATOR);
    }
    out.extend_from_slice(two);
    out
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:465-472
// `DIR_SEPARATORS` is "/" on Unix and "/\\" on Windows (ncbifile.cpp:110, 114).
fn delete_trailing_path_separator(path: &[u8]) -> &[u8] {
    let separator = |byte: &u8| *byte == b'/' || (cfg!(windows) && *byte == b'\\');
    let end = path
        .iter()
        .rposition(|byte| !separator(byte))
        .map_or(0, |at| at + 1);
    &path[..end]
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-972
// `DIR_SEPARATOR` is '/' on Unix and '\\' on Windows, where `DISK_SEPARATOR`
// is ':' (ncbifile.cpp:102-116). `NStr::Split` without flags keeps empty
// pieces. NCBI applies the function to a whole value such as `BLASTDB`, the
// separators of a list included, so `/a/x:/../b` becomes `/a/b`.
//
// The symbolic links are read as NCBI reads them, one path component at a
// time and only outside Windows; a relative path is read relative to the
// current directory of the process. The warning NCBI posts at the depth
// limit is not printed. The path is the bytes of NCBI's `std::string`, a
// link's target too (`CDirEntry::LookupLink`).
/// `CDirEntry::NormalizePath` of `path` with `eFollowLinks`.
fn normalize_path(path: &[u8]) -> Vec<u8> {
    if path.is_empty() {
        return Vec::new();
    }
    let windows = cfg!(windows);
    let sep = SEPARATOR;
    let disk_separator = b':';
    let has_disk = |piece: &[u8]| piece.contains(&disk_separator);

    let mut head: Vec<Vec<u8>> = Vec::new();
    let mut tail: VecDeque<Vec<u8>> = VecDeque::new();
    let mut link_depth = 0;
    let mut follow_links = true;

    let mut current = if windows && has_disk(path) {
        path.to_vec()
    } else {
        delete_trailing_path_separator(path).to_vec()
    };
    if current.is_empty() {
        return vec![sep];
    }
    if windows {
        for byte in &mut current {
            if *byte == b'/' {
                *byte = b'\\';
            }
        }
    }

    while !current.is_empty() || !tail.is_empty() {
        if !current.is_empty() {
            let mut pretail: VecDeque<Vec<u8>> = current
                .split(|&byte| byte == sep)
                .map(<[u8]>::to_vec)
                .collect();
            current.clear();
            let absolute = pretail
                .front()
                .is_some_and(|front| front.is_empty() || (windows && has_disk(front)));
            if absolute {
                head.clear();
                if windows {
                    const UNC: [&[u8]; 4] = [b"", b"", b"?", b"UNC"];
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
                            head.push(Vec::new());
                            head.push(Vec::new());
                        }
                        1 => head.push(Vec::new()),
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
            if last.as_slice() == b"." {
                if !next.is_empty() {
                    head.pop();
                }
            } else if windows && last.last() == Some(&disk_separator) {
                // Allow almost anything right after a volume specification
            } else if next.as_slice() == b"." || next.is_empty() {
                continue;
            } else if next.as_slice() == b".." {
                if last.is_empty() {
                    continue;
                } else if windows && has_disk(last) {
                    // `pos` receives the result of the comparison, 1, so
                    // erasing `last` from the byte after `pos` keeps the first two bytes.
                    last.truncate(2);
                } else if last.as_slice() != b".." {
                    head.pop();
                    continue;
                }
            }
        }
        if follow_links && !windows {
            let link = if head.is_empty() {
                next.clone()
            } else {
                let mut link = head.join(&sep);
                link.push(sep);
                link.extend_from_slice(&next);
                link
            };
            if let Ok(target) = std::fs::read_link(path_of_bytes(&link)) {
                let target = path_bytes(&target).into_owned();
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

    if head.is_empty() || (head.len() == 2 && head[0].as_slice() == b"." && head[1].is_empty()) {
        return b".".to_vec();
    }
    if head.len() == 1 && head[0].is_empty() {
        return vec![sep];
    }
    if windows
        && has_disk(&head[0])
        && ((head.len() == 2 && head[1].as_slice() == b".")
            || (head.len() == 3 && head[1].as_slice() == b"." && head[2].is_empty()))
    {
        let mut root = head[0].clone();
        root.push(sep);
        return root;
    }
    head.join(&sep)
}

/// The separator of the directories of a search path: `s_GetPathSplitter`
/// (seqdbcommon.cpp:270-281), ":" on Unix and ";" elsewhere.
fn path_splitter() -> u8 {
    if cfg!(windows) {
        b';'
    } else {
        b':'
    }
}

/// The directories where NCBI BLAST+ looks for a database name: the current
/// directory, the environment variable `BLASTDB`, and the registry value
/// `[BLAST] BLASTDB` (`.ncbirc`, the program's `.ini` file, or the variable
/// `NCBI_CONFIG__BLAST__BLASTDB`), in this order.
///
/// A relative directory of `BLASTDB` or of the registry value is used as it
/// is written, so the operating system resolves it against the current
/// directory of the process, as for NCBI. The directories are the bytes of
/// NCBI's `std::string` (`utils::os_bytes`), UTF-8 or not.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SearchPath {
    /// `CDir::GetCwd()`: where a volume named in an alias file is looked for
    /// first ("local copy").
    cwd: Vec<u8>,
    /// The string of `CSeqDBAtlas::GenerateSearchPath`.
    path: Vec<u8>,
}

impl SearchPath {
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
    // `env.Get` gives the empty string for a variable that is not set
    // (ncbienv.cpp:109-130), and `NormalizePath` returns an empty string as
    // it is, so the string is `CWD::` without `BLASTDB` and without a
    // registry entry.
    /// The search path of a process whose current directory is `cwd`, whose
    /// environment variable `BLASTDB` is `blastdb` (`None`: not set) and
    /// whose registry has the entry `[BLAST] BLASTDB` with the value
    /// `registry` (`None`: no such entry;
    /// `blastinput::ncbi_environment::ApplicationSettings::blastdb`).
    pub fn new(cwd: &Path, blastdb: Option<&OsStr>, registry: Option<&[u8]>) -> Self {
        let cwd = path_bytes(cwd).into_owned();
        let blastdb = blastdb.map(os_str_bytes).unwrap_or_default();
        let splitter = path_splitter();
        let mut path = normalize_path(&cwd);
        path.push(splitter);
        path.extend_from_slice(&normalize_path(&blastdb));
        path.push(splitter);
        if let Some(registry) = registry {
            path.extend_from_slice(&normalize_path(registry));
            path.push(splitter);
        }
        Self { cwd, path }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3708-3716
    /// The search path of this process: its current directory (the empty
    /// string when it cannot be read, as `CDir::GetCwd` gives), its
    /// environment variable `BLASTDB`, and `registry` as in [`Self::new`].
    pub fn from_environment(registry: Option<&[u8]>) -> Self {
        let cwd = std::env::current_dir().unwrap_or_default();
        let blastdb = std::env::var_os("BLASTDB");
        Self::new(&cwd, blastdb.as_deref(), registry)
    }

    /// The search path that holds the current directory only (`SeqDb::open`: `BLASTDB` and
    /// the registry are not consulted).
    pub(crate) fn current_dir() -> Self {
        Self::new(&std::env::current_dir().unwrap_or_default(), None, None)
    }

    /// The string NCBI names in its "No alias or index file found" message:
    /// the directories, each followed by the separator.
    pub fn as_bytes(&self) -> &[u8] {
        &self.path
    }

    /// `CDir::GetCwd()` held by this search path.
    pub(super) fn cwd_bytes(&self) -> &[u8] {
        &self.cwd
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:330-374
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:232-262
    // `fSplit_Tokenize` is `fSplit_MergeDelimiters` and `fSplit_Truncate` together
    // (ncbistr.hpp:2500-2510): empty directories are left out. With
    // `seqtype` the name is a database name and `s_SeqDB_DBExists` looks for
    // its alias file, then for its index file, inside one directory before
    // the next directory is tried; without it the name is a file name
    // ("exact"). A file exists when it is a regular file (`exists`).
    /// The first directory of the path that holds `dbname`, joined with it.
    fn find(&self, dbname: &[u8], seqtype: Option<SeqType>) -> Option<Vec<u8>> {
        self.find_recording(dbname, seqtype, None)
    }

    /// [`Self::find`], optionally recording every candidate file path tested (`record`).
    /// When `record` is `None` (the default engine path), nothing is recorded or allocated.
    fn find_recording(
        &self,
        dbname: &[u8],
        seqtype: Option<SeqType>,
        mut record: Option<&mut Vec<Vec<u8>>>,
    ) -> Option<Vec<u8>> {
        let mut record_candidate = |path: Vec<u8>| {
            if let Some(log) = record.as_deref_mut() {
                push_unique(log, path);
            }
        };
        for road in split_search_dirs(&self.path) {
            let attempt = combine_path(&fix_delimiters(road), dbname);
            let base = path_of_bytes(&attempt);
            let found = match seqtype {
                Some(seqtype) => DB_LOOKUP_EXTNS.iter().any(|ext| {
                    let candidate = with_extension(&base, seqtype, ext);
                    record_candidate(path_bytes(&candidate).into_owned());
                    exists(&candidate)
                }),
                None => {
                    record_candidate(attempt.clone());
                    exists(&base)
                }
            };
            if found {
                return Some(attempt);
            }
        }
        None
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:453-462
    // `exact` is true, so `s_SeqDB_TryPaths` takes the file name as it is
    // (the `dbtype` `'-'` is not used) and returns the first directory of the
    // search path in which `CFile::GetLength` finds a regular file of that
    // name; the search path is `CSeqDBAtlas::GenerateSearchPath()` (the
    // string of this `SearchPath`).
    /// `SeqDB_ResolveDbPath`: the first directory of the path that holds a
    /// regular file `filename`, joined with it (`None`: NCBI's empty string).
    pub fn resolve_file(&self, filename: &[u8]) -> Option<Vec<u8>> {
        self.find(filename, None)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:464-470
    /// `SeqDB_ResolveDbPathNoExtension` of `dbname` and the type's letter: the first directory of
    /// the path that holds its alias file or its index file, joined with it (`None`: NCBI's
    /// empty string).
    pub(crate) fn resolve_db(&self, dbname: &[u8], seqtype: SeqType) -> Option<Vec<u8>> {
        self.find(dbname, Some(seqtype))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1762-1807
// `keep_quote` is false for a `DBLIST` and for the names the reader opens
// (seqdbalias.cpp:85-88), true for the names of the reports (`report_names`).
// With `keep_quote`, the closing quote stays in the name and NCBI's step of `i` also
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:476-500
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
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:406-436
// `SEQDB_ISEOL` of `x` is true for a newline or a carriage return (seqdbgeneral.hpp:192). A
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:117-139
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1743-1757
// `m_HasGiMask` starts true (seqdbalias.cpp:106). `MASKLIST` anywhere else (an alias file over
// the database, a database among several names, a deeper alias file) is passed over without a
// word: the volumes' mask columns are then the database's masks.
/// The names of the GI-indexed masks that NCBI's reader uses for the top node `root` of the
/// database given as `names` names (`CSeqDBAliasFile::GetMaskList` when `HasGiMask`): one name,
/// whose alias file (the root's one subnode) has `MASKLIST`, whose value is split as a list of
/// database names is. `None` otherwise. Each name is its bytes (NCBI's `std::string`: the
/// file names of the mask and its line of the names table; item K10a of line RG, 2026-11-24).
fn gi_mask_list(root: &AliasNode, names: usize) -> Option<Vec<Vec<u8>>> {
    let [subnode] = root.subnodes.as_slice() else {
        return None;
    };
    if names != 1 {
        return None;
    }
    let value = subnode.value("MASKLIST")?;
    Some(
        split_quoted(value, false)
            .into_iter()
            .map(<[u8]>::to_vec)
            .collect(),
    )
}

struct Resolver<'a, 'r> {
    seqtype: SeqType,
    /// Alias files being expanded, outermost first.
    recurse: Vec<Vec<u8>>,
    /// The directory where a volume named in an alias file is looked for
    /// before the alias file's own directory: NCBI's "local copy", in the
    /// current directory. `None`: the alias file's directory only.
    local_dir: Option<Vec<u8>>,
    /// NCBI's search path, when the name is resolved as a search program
    /// resolves it (`resolve_in`). `None`: the name is the path of one
    /// database (`resolve_one_name`), and names are looked for where they
    /// point only.
    search: Option<&'a SearchPath>,
    /// Optional recorder of every path checked or opened during resolution
    /// (`database_lookup_files`); `None` on the default engine path.
    record: Option<&'r mut Vec<Vec<u8>>>,
}

/// Where the volumes that an alias file names are looked for.
#[derive(Clone, Copy)]
pub(super) enum LocalCopy<'a> {
    /// In the current directory, then in the alias file's (NCBI's rule).
    CurrentDir,
    /// In this directory, then in the alias file's: the rule above with the
    /// directory given by the caller.
    #[cfg_attr(not(test), allow(dead_code))]
    Dir(&'a [u8]),
    /// In the alias file's directory only.
    No,
}

impl Resolver<'_, '_> {
    fn alias_path(&self, base: &[u8]) -> PathBuf {
        with_extension(&path_of_bytes(base), self.seqtype, "al")
    }

    fn index_path(&self, base: &[u8]) -> PathBuf {
        with_extension(&path_of_bytes(base), self.seqtype, "in")
    }

    #[inline]
    fn record_path(&mut self, path: &Path) {
        if let Some(log) = self.record.as_deref_mut() {
            push_unique(log, path_bytes(path).into_owned());
        }
    }

    #[inline]
    fn exists_recorded(&mut self, path: &Path) -> bool {
        self.record_path(path);
        exists(path)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:851-862
    // `m_ExpandLinks` is true for the alias file of a database that is opened
    // (seqdbalias.hpp:1005-1008, the default of `CSeqDBAliasFile`;
    // seqdbimpl.cpp:56 passes no value). The same lines stand at the three
    // other places where NCBI records a volume (seqdbalias.cpp:777-786,
    // 894-903) or a resolved name (seqdbalias.cpp:325-335).
    /// The name under which NCBI records the volume `base`: its directory
    /// normalised, with symbolic links followed. `resolve_one_name` keeps
    /// the name as it is.
    fn normal_name(&self, base: &[u8]) -> Vec<u8> {
        if self.search.is_none() {
            return base.to_vec();
        }
        let mut name = normalize_path(dir_name(base));
        name.push(SEPARATOR);
        name.extend_from_slice(base_name(base));
        name
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:203-222
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:453-460
    // `kSeqDBGroupAliasFileName` is "index.alx" (seqdbcommon.cpp:48), and
    // `FindBlastDBPath` of a `CSeqDB_Path` looks for that file name in the
    // directories of the search path (seqdbalias.hpp:285-300, "exact").
    // NCBI looks a name up in the group alias file before the alias file or
    // the volume of that name, and goes on when the group file has no entry
    // for it. LOSAT does not read group alias files, so the group file NCBI
    // would read is an explicit error whatever it holds.
    /// Rejects the group alias file that NCBI would consult for the alias
    /// file or volume `base`.
    fn reject_group_alias(&mut self, base: &[u8]) -> Result<(), SeqDbError> {
        let aset_path = combine_path(dir_name(base), GROUP_ALIAS_FILE);
        let found = match self.search {
            Some(search) => search.find_recording(&aset_path, None, self.record.as_deref_mut()),
            None => self
                .exists_recorded(&path_of_bytes(&aset_path))
                .then_some(aset_path),
        };
        match found {
            Some(group) => Err(SeqDbError::Unsupported {
                path: path_from_bytes(group),
                message: "the directory has a group alias file; LOSAT does not read group \
                          alias files"
                    .to_string(),
            }),
            None => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:169-176
    /// Reads the alias file of `base` and expands its `DBLIST`.
    fn alias_node(&mut self, base: &[u8]) -> Result<AliasNode, SeqDbError> {
        let path = self.alias_path(base);
        self.record_path(&path);
        // NCBI maps the alias file (`CSeqDBFileMemMap` made from `m_Atlas` and `path.GetPathS()`,
        // seqdbalias.cpp:662-667); a file that cannot be opened gives the
        // errors of `DbFile::open`.
        let bytes = DbFile::open(&path)?.read_all()?;
        let values = read_values(&bytes);
        let mut dblist: &[u8] = &[];
        let mut oidlist: Option<Vec<u8>> = None;
        let dir = dir_name(base);
        for (key, value) in &values {
            if key == b"DBLIST" {
                dblist = value;
            } else if self.record.is_some()
                && PATH_KEYS.iter().any(|k| k.as_bytes() == key.as_slice())
                && !value.is_empty()
            {
                let target = combine_path(dir, &fix_delimiters(value));
                if key == b"OIDLIST" {
                    oidlist = Some(value.clone());
                }
                if let Some(log) = self.record.as_deref_mut() {
                    push_unique(log, target);
                }
            }
        }
        let entries: Vec<Vec<u8>> = split_quoted(dblist, false)
            .into_iter()
            .map(fix_delimiters)
            .collect();

        let this_name = path_bytes(&path).into_owned();
        self.recurse.push(this_name.clone());
        let mut node = AliasNode {
            this_name: Some(this_name),
            values,
            db_path: dir.to_vec(),
            ..AliasNode::default()
        };
        let result = self.expand(&mut node, dir, base, &entries, false);
        self.recurse.pop();
        if let Some(value) = oidlist {
            let direct = combine_path(dir, &fix_delimiters(&value));
            if !exists(&path_of_bytes(&direct)) {
                if let Some(first) = node.vol_names.first() {
                    let fallback = combine_path(dir_name(first), &fix_delimiters(&value));
                    if let Some(log) = self.record.as_deref_mut() {
                        push_unique(log, fallback);
                    }
                }
            }
        }
        result.map(|()| node)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:740-926
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:877-917
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
        db_path: &[u8],
        this_name: &[u8],
        db_list: &[Vec<u8>],
        skip_local: bool,
    ) -> Result<(), SeqDbError> {
        if db_list.is_empty() {
            return Err(SeqDbError::NotFound(
                [
                    b"No database names were found in alias file [".as_slice(),
                    this_name,
                    b"].",
                ]
                .concat(),
            ));
        }
        for entry in db_list {
            if dir_name(entry).is_empty() && base_name(entry) == base_name(this_name) {
                let vol = self.normal_name(this_name);
                if self.record.is_some() {
                    let idx = self.index_path(&vol);
                    self.record_path(&idx);
                }
                push_unique(&mut node.vol_names, vol);
                continue;
            }

            let base = combine_path(db_path, entry);
            let new_db_path = self.alias_path(&base);
            if self
                .recurse
                .iter()
                .any(|name| **name == *path_bytes(&new_db_path))
            {
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
            if self.exists_recorded(&new_db_path) {
                let subnode = self.alias_node(&base)?;
                node.subnodes.push(subnode);
                continue;
            }

            let mut found = None;
            if !skip_local {
                if let Some(local_dir) = &self.local_dir {
                    let local_base = combine_path(local_dir, base_name(entry));
                    let local_idx = self.index_path(&local_base);
                    if let Some(log) = self.record.as_deref_mut() {
                        push_unique(log, path_bytes(&local_idx).into_owned());
                    }
                    if exists(&local_idx) {
                        found = Some(local_base);
                    }
                }
            }
            if found.is_none() {
                let base_idx = self.index_path(&base);
                if self.exists_recorded(&base_idx) {
                    found = Some(base);
                }
            }
            if let Some(volume) = found {
                push_unique(&mut node.vol_names, self.normal_name(&volume));
                continue;
            }

            if let Some(result) = self.search.and_then(|s| {
                s.find_recording(entry, Some(self.seqtype), self.record.as_deref_mut())
            }) {
                let new_alias = self.alias_path(&result);
                if self.exists_recorded(&new_alias) {
                    if self
                        .recurse
                        .iter()
                        .any(|name| **name == *path_bytes(&new_alias))
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
                } else {
                    let idx = self.index_path(&result);
                    if self.exists_recorded(&idx) {
                        push_unique(&mut node.vol_names, self.normal_name(&result));
                    }
                }
                continue;
            }

            return Err(SeqDbError::NotFound(
                [
                    b"Could not find volume or alias file (".as_slice(),
                    entry,
                    b") referenced in alias file (",
                    this_name,
                    b").",
                ]
                .concat(),
            ));
        }
        Ok(())
    }
}

fn push_unique(names: &mut Vec<Vec<u8>>, name: Vec<u8>) {
    if !names.contains(&name) {
        names.push(name);
    }
}

/// Non-empty directory entries of a `BLASTDB`-formatted path list (`SearchPath::find`,
/// `database_lookup_files`, `taxdb::lookup_files_for_use`).
pub(super) fn split_search_dirs(raw: &[u8]) -> impl Iterator<Item = &[u8]> {
    raw.split(|&byte| byte == path_splitter())
        .filter(|road| !road.is_empty())
}

/// Normalised directory path (`normalize_path`) for relativising paths under `cwd`
/// (`database_lookup_files`, `taxdb::lookup_files_for_use`).
pub(super) fn normalize_dir(path: &[u8]) -> Vec<u8> {
    normalize_path(path)
}

/// Strips directory prefix `dir` (with its path separator) from `path`, returning the non-empty
/// relative tail when `path` lies strictly inside `dir`.
pub(super) fn strip_dir_prefix<'a>(path: &'a [u8], dir: &[u8]) -> Option<&'a [u8]> {
    if dir.is_empty() {
        return None;
    }
    let rest = path.strip_prefix(dir)?;
    if dir.last() == Some(&SEPARATOR) {
        (!rest.is_empty()).then_some(rest)
    } else {
        let rest = rest.strip_prefix(&[SEPARATOR])?;
        (!rest.is_empty()).then_some(rest)
    }
}

// No NCBI counterpart: LOSAT Web reads `-logfile` and `windowmasker -out` against the file names
// that `resolve_in` checks or opens for `dbname`, so a write name equal to a missing
// alias/index/group-alias file is refused instead of diverging from the command line, which
// creates the file and then reads it back (line WQA D19b).
/// File names that [`resolve_in`] checks or opens for `dbname` (`seqtype`) using the engine's own
/// [`SearchPath::from_environment`] and [`resolve_in_recording`] walk (returning recorded paths on
/// both `Ok` and `Err`). Paths inside the working directory are returned relative to it so
/// `run_files` matches them without making missing parent folders present.
pub(crate) fn database_lookup_files(
    dbname: &[u8],
    seqtype: SeqType,
    blastdb: Option<&[u8]>,
) -> Vec<Vec<u8>> {
    let search = SearchPath::from_environment(blastdb);
    let mut raw = Vec::new();
    let _ = resolve_in_recording(dbname, seqtype, &search, Some(&mut raw));
    let cwd = fix_delimiters(&search.cwd);
    let norm_cwd = normalize_path(&cwd);
    let mut out = Vec::with_capacity(raw.len());
    for path in raw {
        let stripped = strip_dir_prefix(&path, &norm_cwd)
            .or_else(|| strip_dir_prefix(&path, &cwd))
            .unwrap_or(path.as_slice());
        push_unique(&mut out, stripped.to_vec());
    }
    out
}

fn collect_volumes(node: &AliasNode, out: &mut Vec<Vec<u8>>) {
    out.extend(node.vol_names.iter().cloned());
    for subnode in &node.subnodes {
        collect_volumes(subnode, out);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:190-197
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:932-958
// The volumes are opened, and ordinal ids are given, in the sorted order of
// the volume names, not in `DBLIST` order.
//
// A volume reached from several nodes (an alias file and a name of `-db`, or
// two alias files) is one entry of the set, so it is opened once and its OIDs
// are numbered once; NCBI's walks still count it once per node that names it
// (`walk_totals`), and its OIDs are those of every path to it (`alias_oids`).
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
    AliasTree {
        root,
        volumes,
        gi_masks: None,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:155-165
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:249-337
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:343-383
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:112-127
// The name is a list (`x_Tokenize`, cited at `resolve_one_name`): spaces
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
    dbname: &[u8],
    seqtype: SeqType,
    search: &SearchPath,
) -> Result<AliasTree, SeqDbError> {
    resolve_in_recording(dbname, seqtype, search, None)
}

/// [`resolve_in`], optionally recording every file path checked or opened (`record`) for
/// [`database_lookup_files`]. When `record` is `None`, behavior and cost are unchanged.
fn resolve_in_recording(
    dbname: &[u8],
    seqtype: SeqType,
    search: &SearchPath,
    record: Option<&mut Vec<Vec<u8>>>,
) -> Result<AliasTree, SeqDbError> {
    if dbname.is_empty() {
        return Err(SeqDbError::NotFound(
            "Database name is required.".to_string().into(),
        ));
    }
    let names = split_quoted(dbname, false);
    // `x_ExpandAliases` on an empty list (cited at `Resolver::expand`).
    if names.is_empty() {
        return Err(SeqDbError::NotFound(
            "No database names were passed to CSeqDB::CSeqDB()."
                .to_string()
                .into(),
        ));
    }
    let mut resolver = Resolver {
        seqtype,
        recurse: Vec::new(),
        local_dir: Some(search.cwd.clone()),
        search: Some(search),
        record,
    };

    // Every name is resolved before any is expanded; the first one not found
    // is the error.
    let mut resolved_names = Vec::with_capacity(names.len());
    let mut skip_local = Vec::with_capacity(names.len());
    for name in names {
        let name = fix_delimiters(name);
        skip_local.push(name.contains(&SEPARATOR));

        // `db_path` is the name with the extension of an alias file, and
        // `db_path.FindBasePath()` takes that extension off again.
        let db_path = resolver.alias_path(&name);
        let db_path = path_bytes(&db_path);
        resolver.reject_group_alias(&db_path)?;
        let Some(resolved) = search.find_recording(
            remove_extn(&db_path),
            Some(seqtype),
            resolver.record.as_deref_mut(),
        ) else {
            return Err(not_found(seqtype, &name, search));
        };
        // `resolved_path` is the resolved name with the extension of an alias
        // file; its directory is normalised and its base name loses the extension.
        let resolved_path = resolver.alias_path(&resolved);
        let resolved_path = path_bytes(&resolved_path);
        let mut resolved_name = normalize_path(dir_name(&resolved_path));
        resolved_name.push(SEPARATOR);
        resolved_name.extend_from_slice(remove_extn(base_name(&resolved_path)));
        resolved_names.push(resolved_name);
    }

    // The common prefix of the resolved names, cut back to a whole directory
    // (a separator included), is the directory of the top node.
    let first = resolved_names[0].as_slice();
    let mut common = resolved_names.iter().map(Vec::len).min().unwrap_or(0);
    common = common.saturating_sub(1);
    for other in &resolved_names[1..] {
        while common > 0 && first[..common] != other[..common] {
            common -= 1;
        }
    }
    while common > 0 && first[common - 1] != SEPARATOR {
        common -= 1;
    }
    let db_dir: &[u8] = if common > 1 { &first[..common] } else { b"." };

    let mut root = AliasNode::default();
    // The top node is named "-" (seqdbalias.cpp:106, 127). One `expand` per
    // name: each has its own `m_SkipLocal`, and `expand` keeps no state from
    // one entry to the next but the node's volume list.
    for (name, skip_local) in resolved_names.iter().zip(skip_local) {
        let entry = if common > 1 { &name[common..] } else { name };
        resolver.expand(&mut root, db_dir, b"-", &[entry.to_vec()], skip_local)?;
    }
    let gi_masks = gi_mask_list(&root, resolved_names.len());
    Ok(AliasTree {
        gi_masks,
        ..volumes_in_order(root)
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:311-320
/// The exception of a database name `name` of the molecule `seqtype` that `search` does not hold:
/// the name and the search path are the bytes of NCBI's strings.
fn not_found(seqtype: SeqType, name: &[u8], search: &SearchPath) -> SeqDbError {
    SeqDbError::NotFound(
        [
            b"No alias or index file found for ".as_slice(),
            seqtype.name().as_bytes(),
            b" database [",
            name,
            b"] in search path [",
            search.as_bytes(),
            b"]. Please verify the spelling of the BLAST database and its molecule type.",
        ]
        .concat(),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:85-97
// This entry point is for a caller that holds the path of one database, not
// a list in the syntax above (`losat makeblastdb`, whose output name is one name
// by makeblastdb's own check): the name is not split, so a space or a double
// quote in it is a character of the path. NCBI has no such entry point: its
// reader splits every name it is given.
/// Resolves the one database `db_base` (a path without extension, taken as
/// it is): an alias file when `db_base.{p|n}al` exists, else the single
/// volume `db_base.{p|n}in`. `local_copy` says where the volumes named in an
/// alias file are looked for; `search` is the search path that a database which is not there
/// names in NCBI's message (the files are looked for at `db_base` only). The callers read no
/// masks: `MASKLIST` is passed over (`AliasTree::gi_masks` is `None`).
pub(super) fn resolve_one_name(
    db_base: &[u8],
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
                .map(|cwd| path_bytes(&cwd).into_owned()),
            LocalCopy::Dir(dir) => Some(dir.to_vec()),
            LocalCopy::No => None,
        },
        search: None,
        record: None,
    };
    let db_base = fix_delimiters(db_base);
    if !exists(&resolver.alias_path(&db_base)) && !exists(&resolver.index_path(&db_base)) {
        return Err(not_found(seqtype, &db_base, search));
    }
    let mut root = AliasNode::default();
    // The top node is named "-" (seqdbalias.cpp:106, 127).
    resolver.expand(&mut root, b"", b"-", std::slice::from_ref(&db_base), true)?;
    Ok(volumes_in_order(root))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:970-984
pub(super) fn collect_alias_files(node: &AliasNode, out: &mut std::collections::BTreeSet<Vec<u8>>) {
    if let Some(name) = &node.this_name {
        out.insert(name.clone());
    }
    for subnode in &node.subnodes {
        collect_alias_files(subnode, out);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1517-1532
// The only caller is the removal of an existing database before a build
// (`DeleteBlastDb`), and for that use LOSAT departs from NCBI on purpose in
// two points:
//
// * `db_base` is one name (`resolve_one_name`). NCBI splits it at spaces and
//   double quotes, fails to resolve the pieces of a path that has one, and
//   so leaves the existing database in place.
// * A volume named in an alias file is looked for in the alias file's own
//   directory only (an absolute name as it is). NCBI looks in the current
//   directory first ("local copy", seqdbalias.cpp:831-840, cited at
//   `Resolver::expand`), so `makeblastdb -out sub/db` over an existing
//   `sub/db.pal` removes the files `./db.00.*` of the current directory when
//   they exist, although they are not the files of `sub/db`. The files to
//   remove are those the alias file sits beside.
/// The files of the one database `db_base` for its removal: its volumes
/// (paths without extension, in NCBI's order) and the paths of the alias
/// files that name them (`CSeqDB::FindVolumePaths` with `dbname`, `seqtype`, `paths` and
/// `alias_paths`, with the two differences above).
pub(crate) fn find_volume_paths(
    db_base: &[u8],
    seqtype: SeqType,
) -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>), SeqDbError> {
    let tree = resolve_one_name(db_base, seqtype, LocalCopy::No, &SearchPath::current_dir())?;
    let mut alias_files = std::collections::BTreeSet::new();
    collect_alias_files(&tree.root, &mut alias_files);
    Ok((tree.volumes, alias_files.into_iter().collect()))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1578-1592
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:480-506
// `CSeqDB_TitleWalker` (seqdbalias.cpp:1011-1026) passes a volume's title and
// an alias file's `TITLE` to `SeqDB_JoinDelim`, which joins `m_Value` and `value` with `"; "`.
/// Appends the title of `node` to `title`: its `TITLE` when it has one, else
/// the titles below it joined with "; ".
pub(super) fn walk_title(
    node: &AliasNode,
    volume_title: &dyn Fn(&[u8]) -> Option<Vec<u8>>,
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
// Every volume a node names is in the volume set (`volumes_in_order` takes them all), so
// `GetVol` finds it; a volume that several nodes name is accumulated once per node.
/// What NCBI's `WalkNodes` hands its walker: the value of a node that has the walker's key
/// (`AddString`; the walk does not go below that node) or a volume of a node without it
/// (`Accumulate`).
enum Reached<'a> {
    Value(&'a [u8]),
    Volume(&'a [u8]),
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
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1257-1260
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1290-1305
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1107-1147
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
pub(super) fn walk_totals(root: &AliasNode, volume: &dyn Fn(&[u8]) -> (u32, u64)) -> AliasTotals {
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
        // `new_min` (an `int`) is set from `NStr::StringToUInt` of `value`: a value of 2^31 or more
        // is a negative `int`.
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1963-2014
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:2094-2097
// The rest of the range is cited at `OidRange::bounds`. NCBI's constructor runs this first
// (`m_Aliases.HasFilters()`, seqdbimpl.cpp:96), before any count walk. The path of `OIDLIST` is
// its value joined to the directory of the alias file (`combine_path`: an absolute value stays as
// it is, an empty one is the directory); when that is not a regular file (`exists`) and the node
// names volumes itself, it is the value joined to the directory of the first of them (probes I5,
// I7, I8 and I1b of line RG D1b: a volume found through `BLASTDB` or in the current directory
// moves it there). The file is read where NCBI builds its OID list (`alias_oids`). `OID_MASK_TYPE`
// is converted only beside `OIDLIST` (probe E_abc_no_oidlist of line RG D1b: `OID_MASK_TYPE abc`
// without it is never read), after its path; its errors are those of `FIRST_OID` (probes E_t_abc,
// E_t_4294967296, E_t_empty; a conversion error of `FIRST_OID` comes first, E_order_first).
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:2016-2090
// After `OIDLIST` come `GILIST`, `TILIST`, `SEQIDLIST`, then `MEMB_BIT` (`NStr::StringToUInt`),
// then `TAXIDLIST`, in that exact order. A list value with a space names several lists, which
// NCBI refuses with a `CSeqDBException` naming the alias file's directory (probes b11 f4 and g2 of
// line RG D1c: before any file is read, also under `NSEQ` and `LENGTH`). The path is the value
// joined to that directory, with no other place tried (probes b10 of D1c).
/// `ComputeMasks` of the top node: `FIRST_OID`, then `LAST_OID` (`NStr::StringToUInt`), then the
/// path of `OIDLIST` and its `OID_MASK_TYPE` (`NStr::StringToUInt`), then `GILIST`, `TILIST` and
/// `SEQIDLIST`, then `MEMB_BIT` (`NStr::StringToUInt`), then `TAXIDLIST`, of every node from the
/// top down, a node before its sub-nodes; the first value that does not convert is the error
/// (`SeqDbError::Ncbi`), as is a value naming several lists (`SeqDbError::NotFound`).
pub(super) fn compute_masks(node: &mut AliasNode) -> Result<(), SeqDbError> {
    let convert = |value: &[u8]| string_to_uint(value).map_err(SeqDbError::Ncbi);
    let first_oid = node.value("FIRST_OID").map(convert).transpose()?;
    let last_oid = node.value("LAST_OID").map(convert).transpose()?;
    if first_oid.is_some() || last_oid.is_some() {
        node.range = Some(OidRange {
            first_oid,
            last_oid,
        });
    }
    if let Some(value) = node.value("OIDLIST") {
        let mut path = combine_path(&node.db_path, value);
        if !exists(&path_of_bytes(&path)) {
            if let Some(first) = node.vol_names.first() {
                path = combine_path(dir_name(first), value);
            }
        }
        node.oid_list = Some(path);
        node.oid_mask_type = node
            .value("OID_MASK_TYPE")
            .map(convert)
            .transpose()?
            .unwrap_or(0);
    }
    node.lists.clear();
    for kind in IdListKind::ALL {
        if kind == IdListKind::TaxId {
            node.memb_bit = node.value("MEMB_BIT").map(convert).transpose()?;
        }
        let Some(value) = node.value(kind.key()) else {
            continue;
        };
        if value.contains(&b' ') {
            return Err(SeqDbError::NotFound(
                [
                    b"Alias file (".as_slice(),
                    &node.db_path,
                    b") has multiple ",
                    kind.multiple_words().as_bytes(),
                    b" (",
                    value,
                    b").",
                ]
                .concat(),
            ));
        }
        let path = combine_path(&node.db_path, value);
        node.lists.push(AliasList { kind, path });
    }
    for subnode in &mut node.subnodes {
        compute_masks(subnode)?;
    }
    Ok(())
}

/// Whether some node of the tree has a filter, an OID range, an OID list, an id list or a
/// membership bit (`CSeqDB_FilterTree::HasFilter`).
pub(super) fn has_filters(node: &AliasNode) -> bool {
    node.range.is_some()
        || node.oid_list.is_some()
        || !node.lists.is_empty()
        || node.memb_bit.is_some()
        || node.subnodes.iter().any(has_filters)
}

/// Whether the alias tree under `node` names an id list that NCBI's OID list translates
/// (`GILIST`, `TILIST`, `SEQIDLIST`, and on a version 5 database `TAXIDLIST`): its translation
/// then depends on the search's list ([`super::filter::DeflineRule`]).
pub(super) fn reads_id_lists(node: &AliasNode, version5: bool) -> bool {
    node.lists
        .iter()
        .any(|list| version5 || !matches!(list.kind, IdListKind::TaxId))
        || node
            .subnodes
            .iter()
            .any(|sub| reads_id_lists(sub, version5))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1510-1556
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1598-1611
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:109-112
// The explorer looks at the keys, not at their values (probe K1: `NSEQ abc` beside `OIDLIST`
// is never converted), and the walk does not go below a node whose `Explore` ends it. A
// `TAXIDLIST` asks for the scan on a version 4 database too, where no pass reads it (probe
// b2_v4_tx_info of line RG D1c).
/// Whether NCBI's constructor scans the totals of the included OIDs for the alias tree under
/// `node` (`NeedTotalsScan`): a node with `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST` or
/// `OIDLIST` is reached without passing a node that has both `NSEQ` and `LENGTH`.
pub(super) fn needs_totals_scan(node: &AliasNode) -> bool {
    if node.value("NSEQ").is_some() && node.value("LENGTH").is_some() {
        return false;
    }
    if PATH_KEYS.iter().any(|key| node.value(key).is_some()) {
        return true;
    }
    node.subnodes.iter().any(needs_totals_scan)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:80-99
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:338-348
// `x_Setup` ORs the bit set of each volume (`volume_set`) into `m_AllBits`, a bit array over
// every OID, all clear. That set becomes the special all-set only by the copy of a volume's
// special all-set that covers every OID: the `-db` names reach one volume (that holds OIDs), and
// at least one path of the alias tree to it crosses no `FIRST_OID`-`LAST_OID` range that cuts
// the volume and no OID list (a union of special sets stays special only when one covers the
// others; two ranges whose union is the volume make a bit array; an OID list is a bit array).
// Its `m_Bits` then takes the empty vector of the copy, and keeps the array of the volume set as
// its capacity (`vector::operator=` reuses the storage).
/// What the filters of an alias tree make of NCBI's OID list (`CSeqDBOIDList::x_Setup`,
/// [`alias_oids`]).
#[derive(Debug)]
pub(super) struct AliasOids {
    /// The OIDs the tree includes, one bit per OID of the volume set (bit `oid % 64` of word
    /// `oid / 64`). `None` when no node has a filter: NCBI builds no OID list then, and every
    /// OID is included.
    pub words: Option<Vec<u64>>,
    /// Whether NCBI's bit set of the volume set (`m_AllBits`) is the special all-set one.
    pub all_bits_special_set: bool,
    /// The `OID_MASK_TYPE` that each volume of `volumes` (in their order) is left with
    /// (`CSeqDBVol::SetOidMaskType`, [`volume_set`]): 0 for a volume that no OID list reaches.
    pub mask_types: Vec<u32>,
    /// The `MEMB_BIT` that each volume of `volumes` (in their order) is left with
    /// (`CSeqDBVol::SetMemBit`, [`volume_set`]): 0 for a volume that no membership bit reaches.
    pub memb_bits: Vec<u32>,
    /// The id lists attached to each volume of `volumes` (in their order) for NCBI's filtered
    /// header read (`CSeqDBVol::AttachVolumeGiList`, [`ListSet`]): none for a volume that no
    /// list reaches.
    pub lists: Vec<Vec<std::sync::Arc<AliasIdList>>>,
}

/// A volume of the alias tree as [`alias_oids`] reads it: its path without extension and its
/// global OIDs `[start, end)`.
#[derive(Debug, Clone, Copy)]
pub(super) struct TreeVolume<'a> {
    pub name: &'a [u8],
    pub start: u32,
    pub end: u32,
}

/// The OIDs that the alias tree `root` includes and the shape of NCBI's bit set of them:
/// `volumes` are the volumes in NCBI's order, of the sequence type `seqtype`, `version5` whether
/// NCBI reads the database as version 5 (`IsBlastDBVersion5`), `user` the search's positive list
/// as Pass 1 (`CSeqDBGiListSet`) holds it ([`UserList`]), `user_rule` the search's list rule as
/// Pass 2 reads its Seq-ids and taxids ([`super::filter::DeflineRule`]), and `lmdb` the database's
/// LMDB entries (opened only when version 5 Pass 2 reads a `SEQIDLIST` or `TAXIDLIST`). The OID
/// list files and the id lists are read once each (`files`, [`ListSet`]), in NCBI's order
/// ([`volume_set`], then Pass 2 [`v5_list_filters`]): the first that fails is the error.
#[allow(clippy::too_many_arguments)]
pub(super) fn alias_oids(
    root: &AliasNode,
    volumes: &[TreeVolume<'_>],
    num_oids: u32,
    seqtype: SeqType,
    version5: bool,
    user: UserList<'_>,
    user_rule: Option<&super::filter::DeflineRule>,
    lmdb: &dyn Fn() -> Result<Vec<super::filter::LmdbEntry>, SeqDbError>,
) -> Result<AliasOids, OidListError> {
    let mut words = has_filters(root).then(|| vec![0u64; (num_oids as usize).div_ceil(64)]);
    let mut all = Shape::new(0, u64::from(num_oids), Special::None);
    let mut files = std::collections::HashMap::new();
    let mut lists = ListSet {
        seqtype,
        version5,
        user,
        loaded: Vec::new(),
    };
    let mut mask_types = vec![0u32; volumes.len()];
    let mut memb_bits = vec![0u32; volumes.len()];
    let mut attached = vec![Vec::new(); volumes.len()];
    let mut v5_root_lists: Vec<&[AliasList]> = vec![&[]; volumes.len()];
    for (vol_idx, (((volume, mask_type), memb_bit), attached)) in volumes
        .iter()
        .zip(&mut mask_types)
        .zip(&mut memb_bits)
        .zip(&mut attached)
        .enumerate()
    {
        let Some(node) = specialize(root, volume.name) else {
            continue;
        };
        if version5 {
            v5_root_lists[vol_idx] = node.lists;
        }
        let frame = Frame {
            start: u64::from(volume.start),
            end: u64::from(volume.end),
        };
        let mut reach = VolumeReach {
            volume,
            files: &mut files,
            lists: &mut lists,
            mask_type,
            memb_bit,
            attached,
        };
        let set = volume_set(&node, frame, &mut reach)?;
        all = match all.union_with(set.shape) {
            Union::Left => all,
            Union::Right => set.shape,
            Union::Array(shape) => shape,
        };
        if let Some(words) = &mut words {
            set.or_into(words, frame);
        }
    }
    let made: Vec<std::sync::Arc<AliasIdList>> = lists
        .loaded
        .into_iter()
        .map(|list| std::sync::Arc::new(list.attached()))
        .collect();
    let mut vol_lists: Vec<Vec<std::sync::Arc<AliasIdList>>> = attached
        .into_iter()
        .map(|mut at: Vec<usize>| {
            at.sort_unstable();
            at.dedup();
            at.into_iter().map(|at| made[at].clone()).collect()
        })
        .collect();
    if version5 {
        if let Some(words) = &mut words {
            if v5_list_filters(
                volumes,
                &v5_root_lists,
                num_oids,
                user_rule,
                lmdb,
                words,
                &mut vol_lists,
            )? {
                all.special = Special::None;
            }
        }
    }
    Ok(AliasOids {
        words,
        all_bits_special_set: all.special == Special::AllSet,
        mask_types,
        memb_bits,
        lists: vol_lists,
    })
}

/// What the walk of one volume ([`volume_set`]) reads and leaves: the OID list files read so far,
/// the id lists, the volume's `OID_MASK_TYPE`, its `MEMB_BIT` and the id lists it attaches
/// (indexes of `ListSet::loaded`).
struct VolumeReach<'a, 'u, 'r> {
    volume: &'r TreeVolume<'a>,
    files: &'r mut std::collections::HashMap<&'a [u8], OidListFile>,
    lists: &'r mut ListSet<'a, 'u>,
    mask_type: &'r mut u32,
    memb_bit: &'r mut u32,
    attached: &'r mut Vec<usize>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:394-466
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2501-2551
// `GetNodeIdList` reads a list's file at its first use, one list per kind and path, and at
// each use translates the list through the volume (`IdsToOids`) and attaches it to the volume.
// A list of ids needs the volume's ISAM of their kind (the volume's string ISAM for Seq-ids; the
// volume opens an ISAM only where it has OIDs): without one NCBI throws `<GI|TI|SI> list
// specified but no ISAM file found for <GI|TI|SI> in <volume>`. Whether the volumes translate
// the list depends on the search's positive list (`UserList`).
/// The id lists of the alias files as NCBI's `CSeqDBGiListSet` holds them for the OID list
/// ([`volume_set`]).
struct ListSet<'a, 'u> {
    seqtype: SeqType,
    version5: bool,
    user: UserList<'u>,
    loaded: Vec<LoadedList<'a>>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:417-423
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:447-461
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:472-526
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:588-592
// With any positive user list, `GetNodeIdList` runs `x_TranslateFromUserList` at the list's first
// use, which starts by putting the list in order (`InsureOrder`, `ListIds::insure_order`): a list
// of Seq-ids is then sorted by its ids as the file gives them, before any simplification. For a
// positive GI list, matching GIs of an alias `GILIST` take their OIDs directly from the user
// list (`x_TranslateGisFromUserList`). The volumes then translate the list (`mixed_ids`) only
// when the user list holds Seq-ids or taxids, or when one holds GIs and the other TIs: a positive
// GI list does not translate an alias `GILIST` or v4 `SEQIDLIST` through the volumes' ISAMs, a
// positive IPG list translates none, and a version 5 database's `-seqidlist` whose every id is a
// gi (`s_ProcessSeqIDsForV5` drops them, seqdbgilistset.cpp:148-180) leaves every list's OIDs
// unset and its ISAMs unread.
/// What NCBI's `CSeqDBGiListSet` holds of the search's positive list (`m_UserList`) when it reads
/// the alias files' id lists ([`alias_oids`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum UserList<'a> {
    /// No positive list: the database opened without a list, or with a negative one.
    None,
    /// A positive list (`-seqidlist`, `-taxids`, `-taxidlist`) of `ids` Seq-ids and taxids, as
    /// NCBI counts them after its processing of the list.
    Positive { ids: usize },
    /// A positive GI list (`-gilist`): its sorted GIs with their resolved OIDs.
    PositiveGis(&'a [super::isam_num::NumOid]),
    /// A positive IPG list (`-ipglist`): holds PIGs alone (`GetNumGis`, `GetNumTis`, `GetNumSis`
    /// and `GetNumTaxIds` are all 0).
    PositivePigs,
}

impl UserList<'_> {
    /// Whether the volumes translate an alias file's list (`m_UserList` empty or `mixed_ids` set).
    fn translates(self, list: &ListIds) -> bool {
        match self {
            UserList::None => true,
            UserList::Positive { ids } => ids > 0,
            UserList::PositiveGis(gis) => !gis.is_empty() && matches!(list, ListIds::Ti(_)),
            UserList::PositivePigs => false,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:474-526
/// `x_TranslateGisFromUserList`: copies resolved OIDs from the sorted positive user GI list
/// `source` into matching unset entries of the sorted alias GI list `target`.
fn translate_gis_from_user_list(
    source: &[super::isam_num::NumOid],
    target: &mut [super::isam_num::NumOid],
) {
    let source_num = source.len();
    let target_num = target.len();
    let (mut source_index, mut target_index) = (0usize, 0usize);
    while source_index < source_num && target_index < target_num {
        let source_gi = source[source_index].key;
        let target_gi = target[target_index].key;
        if source_gi == target_gi {
            if target[target_index].oid.is_none() {
                target[target_index].oid = source[source_index].oid;
            }
            target_index += 1;
            source_index += 1;
        } else if source_gi > target_gi {
            target_index += 1;
            let mut jump = 2;
            let mut test = target_index + jump;
            while test < target_num && target[test].key < source_gi {
                target_index = test;
                jump *= 2;
                test = target_index + jump;
            }
        } else {
            source_index += 1;
            let mut jump = 2;
            let mut test = source_index + jump;
            while test < source_num && source[test].key < target_gi {
                source_index = test;
                jump *= 2;
                test = source_index + jump;
            }
        }
    }
}

/// An id list read from its file, with the translation of its ids so far.
struct LoadedList<'a> {
    kind: IdListKind,
    path: &'a [u8],
    ids: ListIds,
    /// The ids the list counts (`GetNumGis`, `GetNumTis`, `GetNumSis`), which decides whether a
    /// volume is asked for its ISAM.
    count: usize,
}

/// The ids of a [`LoadedList`].
enum ListIds {
    /// The GIs, sorted (NCBI puts them in order at their first translation, `InsureOrder`, and
    /// nothing changes them after it), with the OIDs found so far.
    Gi(Vec<super::isam_num::NumOid>),
    /// The TIs, as NCBI's unsigned `TTi`, sorted (as `FindTi` sorts them); a number past `Int8`
    /// is none.
    Ti(Vec<u64>),
    /// The Seq-ids in NCBI's order, as the file gives them until the first translation
    /// simplifies them (`PreprocessIdsForISAMSiLookup`), with the OIDs found so far, and NCBI's
    /// order flag (`m_CurrentOrder` at `eGi`, `isam::insure_order`).
    Si {
        ids: Vec<super::isam::SiOid>,
        sorted: bool,
    },
}

impl ListIds {
    /// `InsureOrder` with `eGi` of the list (`x_TranslateFromUserList`): the GIs and TIs are kept
    /// sorted already.
    fn insure_order(&mut self) {
        if let ListIds::Si { ids, sorted } = self {
            super::isam::insure_order(ids, sorted);
        }
    }
}

impl LoadedList<'_> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:745-767
    // The filtered header read looks the Seq-ids up in the list as the last translation left it
    // (`FindSi`), whose first lookup puts it in order where its flag is not set yet (no volume
    // translated it): every lookup comes after the last translation.
    /// The list as the volumes it reached hold it for the filtered header read.
    fn attached(self) -> AliasIdList {
        let mut list = AliasIdList::default();
        match self.ids {
            ListIds::Gi(ids) => list.gis = ids.into_iter().map(|id| id.key).collect(),
            ListIds::Ti(tis) => list.tis = tis,
            ListIds::Si {
                mut ids,
                mut sorted,
            } => {
                super::isam::insure_order(&mut ids, &mut sorted);
                list.sis = PackedStrings::from_vecs(ids.into_iter().map(|id| id.si));
            }
        }
        list
    }
}

impl<'a> ListSet<'a, '_> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:86-131
    // The file is mapped as an OID list file is (`not_a_file`), then read as its kind's list
    // (`id_list`; a seqid list as text, `seqidlist::read_memory_si_list`, its ids as the file
    // gives them, not in order). LOSAT refuses a token that NCBI's text reader reports (it writes
    // a warning).
    /// The index in `loaded` of the list of `kind` at `path`, read at its first use.
    fn load(&mut self, kind: IdListKind, path: &'a [u8]) -> Result<usize, OidListError> {
        if let Some(at) = self
            .loaded
            .iter()
            .position(|list| list.kind == kind && list.path == path)
        {
            return Ok(at);
        }
        let at = &*path_of_bytes(path);
        if !exists(at) {
            return Err(OidListError::NotAFile(path.to_vec()));
        }
        let keep = OidListError::reader;
        let file = DbFile::open(at).map_err(keep)?;
        let bytes = if file.is_empty() {
            Vec::new()
        } else {
            file.read_all().map_err(keep)?
        };
        let list_error =
            |text: String| OidListError::reader(SeqDbError::NotFound(text.into_bytes()));
        let (mut ids, count) = match kind {
            IdListKind::Gi => {
                let mut gis: Vec<super::isam_num::NumOid> = super::id_list::read_gi_list(&bytes)
                    .map_err(list_error)?
                    .into_iter()
                    .map(|gi| super::isam_num::NumOid {
                        key: i64::from(gi),
                        oid: None,
                    })
                    .collect();
                gis.sort_by_key(|id| id.key);
                let count = gis.len();
                (ListIds::Gi(gis), count)
            }
            IdListKind::Ti => {
                let read = super::id_list::read_ti_list(&bytes).map_err(list_error)?;
                let count = read.len();
                let mut tis: Vec<u64> = read
                    .into_iter()
                    .filter_map(|ti| match ti {
                        super::id_list::ListNumber::Value(ti) => Some(ti as u64),
                        super::id_list::ListNumber::Past => None,
                    })
                    .collect();
                tis.sort_unstable();
                (ListIds::Ti(tis), count)
            }
            IdListKind::Si => {
                let mut warnings = Vec::new();
                let read = super::seqidlist::read_memory_si_list(&bytes, &mut warnings);
                if let Some(warning) = warnings.first() {
                    return Err(OidListError::reader(SeqDbError::Unsupported {
                        path: at.to_path_buf(),
                        message: format!(
                            "the seqid list of an alias file's SEQIDLIST has a token NCBI reports \
                             (`{}`)",
                            warning.trim()
                        ),
                    }));
                }
                let ids: Vec<super::isam::SiOid> =
                    read.into_iter().map(super::isam::SiOid::new).collect();
                let count = ids.len();
                (ListIds::Si { ids, sorted: false }, count)
            }
            IdListKind::TaxId => unreachable!("no pass of a version 4 database reads TAXIDLIST"),
        };
        if !matches!(self.user, UserList::None) {
            ids.insure_order();
            if let (UserList::PositiveGis(source), ListIds::Gi(target)) = (self.user, &mut ids) {
                translate_gis_from_user_list(source, target);
            }
        }
        self.loaded.push(LoadedList {
            kind,
            path,
            ids,
            count,
        });
        Ok(self.loaded.len() - 1)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:407-460
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1388-1419
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:155
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:128
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:447-451
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:620-630
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:548-562
    // A TI list is translated by the TI ISAM's `x_TranslateGiList` of the ISAM's own `TTi`, a
    // signed `Int8`, while the list's specialisations of its size and keys are of the global
    // `TTi`, an unsigned `Uint8`: the translation takes the generic size, the count of the list's
    // GIs, which a TI list has none of, and returns before it looks at the ISAM. So no TI of any
    // list finds an OID on any TI ISAM (an uninitialised one included: no batch-mode error), and
    // only the ISAM's opening can fail (`NumericIsam::open_for_list`; probes b4b, b4c and tiself
    // of line RG D1c agree). A seqid list is simplified again before each translation
    // (`PreprocessIdsForISAMSiLookup`): LOSAT refuses an id that NCBI's simplification throws on.
    /// The OIDs of `volume` that the list `at` holds (`IdsToOids`, then `x_IdsToBitSet`), as bits
    /// over the volume's OIDs.
    fn volume_bits(
        &mut self,
        at: usize,
        volume: &TreeVolume<'_>,
    ) -> Result<Vec<u64>, OidListError> {
        let keep = OidListError::reader;
        let no_isam = |what: &str| {
            OidListError::reader(SeqDbError::NotFound(
                [
                    format!("{what} list specified but no ISAM file found for {what} in ")
                        .as_bytes(),
                    volume.name,
                ]
                .concat(),
            ))
        };
        let num_oids = volume.end - volume.start;
        let frame = Frame {
            start: u64::from(volume.start),
            end: u64::from(volume.end),
        };
        let list = &mut self.loaded[at];
        if list.count == 0 {
            return Ok(frame.clear_words());
        }
        let bits_of = |oids: &mut dyn Iterator<Item = Option<u32>>| {
            let mut bits = frame.clear_words();
            for oid in oids.flatten() {
                if oid >= volume.start && oid < volume.end {
                    let at = oid - volume.start;
                    bits[(at / 64) as usize] |= 1u64 << (at % 64);
                }
            }
            bits
        };
        if !self.user.translates(&list.ids) {
            return Ok(match &list.ids {
                ListIds::Gi(gis) => bits_of(&mut gis.iter().map(|gi| gi.oid)),
                ListIds::Ti(_) | ListIds::Si { .. } => frame.clear_words(),
            });
        }
        let base = &*path_of_bytes(volume.name);
        match &mut list.ids {
            ListIds::Gi(_) | ListIds::Ti(_) => {
                let (kind, what) = match list.ids {
                    ListIds::Gi(_) => (super::isam_num::ListIsam::Gi, "GI"),
                    _ => (super::isam_num::ListIsam::Ti, "TI"),
                };
                let isam =
                    super::isam_num::NumericIsam::open_for_list(base, self.seqtype, kind, num_oids)
                        .map_err(keep)?
                        .ok_or_else(|| no_isam(what))?;
                match (&mut list.ids, isam) {
                    (ListIds::Gi(gis), super::isam_num::Opened::Ready(isam)) => {
                        isam.ids_to_oids(volume.start, gis).map_err(keep)?;
                        Ok(bits_of(&mut gis.iter().map(|gi| gi.oid)))
                    }
                    (ListIds::Gi(_), super::isam_num::Opened::Uninitialized) => {
                        Err(keep(super::isam::batch_mode_error()))
                    }
                    _ => Ok(frame.clear_words()),
                }
            }
            ListIds::Si { ids, sorted } => {
                for id in ids.iter_mut() {
                    id.si = super::isam::simplify_accession_for_isam(&id.si).map_err(|error| {
                        OidListError::reader(SeqDbError::NcbiThrows {
                            path: path_from_bytes(list.path),
                            message: format!(
                                "the id {} of the seqid list of an alias file's SEQIDLIST throws \
                                 in NCBI's SeqDB_SimplifyAccession: {}",
                                String::from_utf8_lossy(&id.si),
                                error.to_string().trim_end()
                            ),
                        })
                    })?;
                }
                let protein = matches!(self.seqtype, SeqType::Protein);
                let isam = if num_oids == 0 {
                    None
                } else {
                    super::isam::StringIsam::open_existing(base, protein).map_err(keep)?
                };
                let Some(isam) = isam else {
                    return Err(no_isam("SI"));
                };
                isam.ids_to_oids(volume.start, ids, sorted).map_err(keep)?;
                Ok(bits_of(&mut ids.iter().map(|id| id.oid)))
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:130-177
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:191-207
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:219-234
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:290-291
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:464-505
// For each volume in order, the tree specialised to it (`specialize`): a node that names the
// volume starts from all of its OIDs and does not look at its sub-nodes, any other node from the
// union of its sub-nodes; the node's filter, every OID of the volume ANDed with its range and
// then with its OID list, is then ANDed in. So the OID lists are read in that order, a
// sub-node's before its parent's, and only on the branches that reach the volume without passing
// a node that names it: the first file that fails is NCBI's error (`oid_list::OidListError`).
// An OID list is read relative to each volume (bit `i` is the volume's OID `i`, probe F1 of line
// RG D1b), a bit array over the volume's OIDs (`x_GetOidMask`); the ranges count the OIDs of the
// whole volume set. Each OID list sets the volume's `OID_MASK_TYPE` to its node's
// (`SetOidMaskType`, which overwrites): the last one read wins, a parent's over its sub-nodes', a
// later sub-node's over an earlier one's, and a node that names the volume reads none below it
// (probes J2, J4a, J4b, J3 of line RG D1b). The `int` bounds of a range become `size_t`: a
// negative one is sign-extended past every OID (line RG D1a-2, survey A M2), so a negative
// `FIRST_OID` keeps no OID and a negative `LAST_OID` is no upper bound, as `bound` reads them.
// After the OID list come the node's id lists, `GILIST`, `TILIST` and `SEQIDLIST` in that order
// (seqdboidlist.cpp:211-269): each a bit array of the OIDs its ids are found in (`ListSet`), the
// list attached to the volume. On a version 5 database, `SEQIDLIST` and `TAXIDLIST` are read in
// Pass 2 ([`v5_list_filters`]); on a version 4 database, NCBI reads `TAXIDLIST` nowhere (probe
// b2_v4_tx_info of line RG D1c).
/// NCBI's bit set of the OIDs of one volume (`frame`) that the specialised tree `node` keeps
/// (`CSeqDBOIDList::x_ComputeFilters`); `reach` holds the files read so far and takes the
/// volume's `OID_MASK_TYPE`, set by each OID list read, and the id lists the volume attaches.
fn volume_set<'a>(
    node: &VolumeNode<'a>,
    frame: Frame,
    reach: &mut VolumeReach<'a, '_, '_>,
) -> Result<OidSet, OidListError> {
    let map = if node.volumes > 0 {
        OidSet::special(frame.start, frame.end, Special::AllSet)
    } else {
        let mut map = OidSet::special(frame.start, frame.end, Special::AllClear);
        for sub in &node.nodes {
            let sub = volume_set(sub, frame, reach)?;
            map = map.union_with(sub, frame);
        }
        map
    };
    let mut filter = OidSet::special(frame.start, frame.end, Special::AllSet);
    if let Some(range) = node.range {
        let (first, last) = range.bounds();
        let range = OidSet::special(bound(first), bound(last), Special::AllSet);
        filter = filter.intersect_with(range, frame);
    }
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:208-215
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:879-885
    // `SetMemBit` runs in the first filter loop with `eOidRange`, before any OID list or id list
    // is read: a non-zero bit already stored on the volume conflicts with any different bit
    // (including 0 from a parent), while 0 allows any later bit to replace it.
    if let Some(bit) = node.memb_bit {
        if *reach.memb_bit != 0 && bit != *reach.memb_bit {
            return Err(OidListError::reader(SeqDbError::NotFound(
                b"MEMB_BIT error: conflicting bit found.".to_vec(),
            )));
        }
        *reach.memb_bit = bit;
    }
    if let Some(path) = node.oid_list {
        let file = match reach.files.entry(path) {
            std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::hash_map::Entry::Vacant(entry) => entry.insert(oid_list::read(path)?),
        };
        let bits = file.volume_words((frame.end - frame.start) as u32)?;
        filter = filter.intersect_with(OidSet::array(frame, bits), frame);
        *reach.mask_type = node.oid_mask_type;
    }
    for list in node.lists {
        let read = match list.kind {
            IdListKind::Gi | IdListKind::Ti => true,
            IdListKind::Si => !reach.lists.version5,
            IdListKind::TaxId => false,
        };
        if !read {
            continue;
        }
        let at = reach.lists.load(list.kind, &list.path)?;
        let bits = reach.lists.volume_bits(at, reach.volume)?;
        reach.attached.push(at);
        filter = filter.intersect_with(OidSet::array(frame, bits), frame);
    }
    if map.clears_outside(&filter) {
        return Err(OidListError::RangeOutsideBits(path_from_bytes(
            node.name.unwrap_or(b"-"),
        )));
    }
    Ok(map.intersect_with(filter, frame))
}

/// The global OIDs `[start, end)` of the volume whose bit sets are computed ([`volume_set`]).
#[derive(Debug, Clone, Copy)]
struct Frame {
    start: u64,
    end: u64,
}

impl Frame {
    /// An empty array of bits over the volume's OIDs (bit `oid - start`).
    fn clear_words(self) -> Vec<u64> {
        vec![0u64; ((self.end - self.start) as usize).div_ceil(64)]
    }

    /// The bits over the volume's OIDs of the global OIDs `[from, to)`.
    fn range_words(self, from: u64, to: u64) -> Vec<u64> {
        let mut words = self.clear_words();
        let from = from.max(self.start) - self.start;
        let to = to.min(self.end).saturating_sub(self.start);
        set_bits(&mut words, from, to);
        words
    }
}

/// ORs the volume-relative bits `bits` into the global bits `words` from OID `start`.
fn or_at(words: &mut [u64], bits: &[u64], start: u32) {
    let shift = start % 64;
    let first = (start / 64) as usize;
    for (at, &word) in bits.iter().enumerate() {
        if word == 0 {
            continue;
        }
        words[first + at] |= word << shift;
        if shift > 0 && first + at + 1 < words.len() {
            words[first + at + 1] |= word >> (64 - shift);
        }
    }
}

/// NCBI's `size_t` of an `int` bound of an OID range (the `CSeqDB_BitSet` constructor, given
/// `mask.GetBegin()`).
fn bound(value: i32) -> u64 {
    i64::from(value) as u64
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:746
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:376-408
/// Clears the bits `[from, to)` of `words` (nothing when `to <= from`), matching
/// `CSeqDB_BitSet::AssignBitRange` with `false`.
fn clear_bits(words: &mut [u64], from: u64, to: u64) {
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
        words[word] &= !ones;
        at += span;
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:544-563
/// `s_AddFilterFile`: records that `vol` uses the filter file `name`, preserving first-seen file
/// order.
fn add_filter_file<'a>(
    name: &'a [u8],
    vol: &'a [u8],
    fnames: &mut Vec<&'a [u8]>,
    fnames_vols: &mut Vec<Vec<&'a [u8]>>,
) {
    if let Some(pos) = fnames.iter().position(|&existing| existing == name) {
        fnames_vols[pos].push(vol);
    } else {
        fnames.push(name);
        fnames_vols.push(vec![vol]);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:520-542
/// `s_GetFilteredOidRange`: attaches `list` to every volume named in `vols`, and returns the
/// `[start, end)` OID ranges of all other volumes (`excluded_vols`).
fn attach_and_excluded(
    volumes: &[TreeVolume<'_>],
    vols: &[&[u8]],
    list: std::sync::Arc<AliasIdList>,
    vol_lists: &mut [Vec<std::sync::Arc<AliasIdList>>],
) -> Vec<(u32, u32)> {
    let mut excluded = Vec::new();
    for (idx, volume) in volumes.iter().enumerate() {
        if vols.contains(&volume.name) {
            vol_lists[idx].push(std::sync::Arc::clone(&list));
        } else {
            excluded.push((volume.start, volume.end));
        }
    }
    excluded
}

const ALIAS_SEQIDLIST_KIND: &str = "an alias file's SEQIDLIST";
const USER_SEQIDLIST_KIND: &str = "the search's -seqidlist";
const USER_NEG_SEQIDLIST_KIND: &str = "the search's -negative_seqidlist";

fn optional_path(list_path: Option<&[u8]>) -> PathBuf {
    list_path.map(path_from_bytes).unwrap_or_default()
}

fn display_seqid_token(bytes: &[u8]) -> String {
    if bytes.iter().all(|&b| (0x20..0x7f).contains(&b)) {
        String::from_utf8_lossy(bytes).into_owned()
    } else {
        bytes.escape_ascii().to_string()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:565-576
/// Parses one Seq-id string as `s_CompareSeqId` parses it (`fParse_AnyRaw` and
/// `fParse_ValidLocal`), returning NCBI's `CSeqIdException` on parse failure and refusing tokens
/// that post a Seq-id warning.
fn parse_compare_seq_id(
    id: &[u8],
    list_path: Option<&[u8]>,
    list_kind: &str,
) -> Result<crate::seqid::ParsedSeqId, OidListError> {
    let flags = crate::seqid::accguide::parse_flags::ANY_RAW
        | crate::seqid::accguide::parse_flags::VALID_LOCAL;
    let mut warnings = Vec::new();
    let parsed = crate::seqid::parse::seq_id_from_text(id, flags, &mut |w| warnings.push(w))
        .map_err(|error| OidListError::reader(SeqDbError::Ncbi(error.exception())))?;
    if let Some(warning) = warnings.first() {
        return Err(OidListError::reader(SeqDbError::Unsupported {
            path: optional_path(list_path),
            message: format!(
                "comparing seqids of {list_kind} with a search list posts a Seq-id warning (`{}`)",
                String::from_utf8_lossy(&warning.message)
            ),
        }));
    }
    Ok(parsed)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:565-576
/// `s_CompareSeqId`: compares two Seq-id strings using `CSeq_id` with flags `fParse_AnyRaw` and
/// `fParse_ValidLocal`; identical strings or matching `CSeq_id`s compare equivalent (`false`),
/// otherwise bytewise lexicographical order.
fn compare_seq_id(
    id1: &[u8],
    path1: Option<&[u8]>,
    kind1: &str,
    id2: &[u8],
    path2: Option<&[u8]>,
    kind2: &str,
) -> Result<bool, OidListError> {
    if id1 == id2 {
        return Ok(false);
    }
    let seq1 = parse_compare_seq_id(id1, path1, kind1)?;
    let seq2 = parse_compare_seq_id(id2, path2, kind2)?;
    if crate::seqid::defline_set::compare(&seq1, &seq2) == crate::seqid::defline_set::Sic::Yes {
        return Ok(false);
    }
    Ok(id1 < id2)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:565-576,592,611
/// Scans all unique entries of a bytewise-sorted list `items` (with at least two distinct entries)
/// as `s_CompareSeqId` parses them, collecting compact `(match_key, index)` pairs while checking
/// in order: (1) any Seq-id warning (refused), (2) two or more distinct exception texts (refused),
/// (3) one common exception text (raised as NCBI's `CSeqIdException`).
fn scan_compare_seqids<S: SeqIdSlice + ?Sized>(
    items: &S,
    list_path: Option<&[u8]>,
    list_kind: &str,
) -> Result<Vec<(u64, u32)>, OidListError> {
    let flags = crate::seqid::accguide::parse_flags::ANY_RAW
        | crate::seqid::accguide::parse_flags::VALID_LOCAL;
    let mut first_warning: Option<Vec<u8>> = None;
    let mut first_exc: Option<(usize, crate::ncbi_exception::NcbiException, Vec<u8>)> = None;
    let mut second_exc_idx: Option<usize> = None;
    let mut key_entries: Vec<(u64, u32)> = Vec::new();
    let mut tmp_keys: Vec<u64> = Vec::new();
    for idx in 0..items.len() {
        let bytes = items.get(idx);
        if idx > 0 && bytes == items.get(idx - 1) {
            continue;
        }
        let mut item_warning: Option<Vec<u8>> = None;
        let parsed = crate::seqid::parse::seq_id_from_text(bytes, flags, &mut |w| {
            if item_warning.is_none() {
                item_warning = Some(w.message);
            }
        });
        if first_warning.is_none() {
            if let Some(warn_msg) = item_warning {
                first_warning = Some(warn_msg);
                key_entries = Vec::new();
            }
        }
        match parsed {
            Ok(seq) if first_warning.is_none() && first_exc.is_none() => {
                tmp_keys.clear();
                crate::seqid::defline_set::match_keys(&seq, &mut tmp_keys);
                let idx_u32 = u32::try_from(idx).map_err(|_| {
                    OidListError::reader(SeqDbError::Unsupported {
                        path: optional_path(list_path),
                        message: format!("{list_kind} has more than u32::MAX entries"),
                    })
                })?;
                for &k in &tmp_keys {
                    key_entries.push((k, idx_u32));
                }
            }
            Ok(_) => {}
            Err(error) => {
                let exc = error.exception();
                let c_what = crate::cli::c_string(&exc.ncbi_what()).to_vec();
                match &first_exc {
                    None => {
                        key_entries = Vec::new();
                        first_exc = Some((idx, exc, c_what));
                    }
                    Some((_, _, prev_what)) => {
                        if second_exc_idx.is_none() && *prev_what != c_what {
                            second_exc_idx = Some(idx);
                        }
                    }
                }
            }
        }
    }
    if let Some(warn_msg) = first_warning {
        return Err(OidListError::reader(SeqDbError::Unsupported {
            path: optional_path(list_path),
            message: format!(
                "comparing seqids of {list_kind} with a search list posts a Seq-id warning (`{}`)",
                String::from_utf8_lossy(&warn_msg)
            ),
        }));
    }
    if let (Some((first_idx, _, _)), Some(second_idx)) = (&first_exc, second_exc_idx) {
        return Err(OidListError::reader(SeqDbError::Unsupported {
            path: optional_path(list_path),
            message: format!(
                "{list_kind} holds Seq-ids (`{}` and `{}`) that fail NCBI's `s_CompareSeqId` with \
                 different exceptions; combining a version 5 alias SEQIDLIST with a search seqid \
                 list then depends on the C++ library's sort order",
                display_seqid_token(items.get(*first_idx)),
                display_seqid_token(items.get(second_idx)),
            ),
        }));
    }
    if let Some((_, exc, _)) = first_exc {
        return Err(OidListError::reader(SeqDbError::Ncbi(exc)));
    }
    Ok(key_entries)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:565-576,592,611
/// Rejects a bytewise-sorted list `items` when it contains warnings, multiple distinct
/// `s_CompareSeqId` exceptions, or two differently spelled Seq-ids that `s_CompareSeqId` treats
/// as equivalent (`CSeq_id::Match`). Only key-colliding runs are re-parsed, keeping peak memory
/// linear in the number of `(u64, u32)` keys rather than holding all parsed Seq-ids at once.
fn validate_compare_seqids<S: SeqIdSlice + ?Sized>(
    items: &S,
    list_path: Option<&[u8]>,
    list_kind: &str,
) -> Result<(), OidListError> {
    if items.len() <= 1 || items.get(0) == items.get(items.len() - 1) {
        return Ok(());
    }
    let mut key_entries = scan_compare_seqids(items, list_path, list_kind)?;
    key_entries.sort_unstable();
    let mut run_start = 0;
    while run_start < key_entries.len() {
        let key = key_entries[run_start].0;
        let mut run_end = run_start + 1;
        while run_end < key_entries.len() && key_entries[run_end].0 == key {
            run_end += 1;
        }
        if run_end - run_start >= 2 {
            let mut parsed_run: Vec<(usize, crate::seqid::ParsedSeqId)> =
                Vec::with_capacity(run_end - run_start);
            for &(_, idx_u32) in &key_entries[run_start..run_end] {
                let idx = idx_u32 as usize;
                if parsed_run
                    .last()
                    .is_some_and(|&(prev_idx, _)| prev_idx == idx)
                {
                    continue;
                }
                let seq = parse_compare_seq_id(items.get(idx), list_path, list_kind)?;
                for &(prev_idx, ref prev_seq) in &parsed_run {
                    if crate::seqid::defline_set::compare(prev_seq, &seq)
                        == crate::seqid::defline_set::Sic::Yes
                    {
                        return Err(OidListError::reader(SeqDbError::Unsupported {
                            path: optional_path(list_path),
                            message: format!(
                                "{list_kind} holds two differently spelled Seq-ids (`{}` and \
                                 `{}`) that NCBI's `s_CompareSeqId` treats as equivalent; \
                                 combining a version 5 alias SEQIDLIST with a search seqid list \
                                 then depends on the C++ library's sort order",
                                display_seqid_token(items.get(prev_idx)),
                                display_seqid_token(items.get(idx)),
                            ),
                        }));
                    }
                }
                parsed_run.push((idx, seq));
            }
        }
        run_start = run_end;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:612-622
/// Sorted set intersection of `lhs` and `rhs` under fallible directional comparators `less_lr`
/// and `less_rl` (keeps matching elements from `lhs`, matching `std::set_intersection` in
/// `s_ProcessSeqIdFilters`).
fn set_intersection_by<'a, L, R, F, G>(
    lhs: &'a L,
    rhs: &R,
    less_lr: &mut F,
    less_rl: &mut G,
) -> Result<Vec<&'a [u8]>, OidListError>
where
    L: SeqIdSlice + ?Sized,
    R: SeqIdSlice + ?Sized,
    F: FnMut(&[u8], &[u8]) -> Result<bool, OidListError>,
    G: FnMut(&[u8], &[u8]) -> Result<bool, OidListError>,
{
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < lhs.len() && j < rhs.len() {
        let a = lhs.get(i);
        let b = rhs.get(j);
        if less_lr(a, b)? {
            i += 1;
        } else if less_rl(b, a)? {
            j += 1;
        } else {
            out.push(a);
            i += 1;
            j += 1;
        }
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:623-633
/// Sorted set difference of `lhs` and `rhs` under fallible directional comparators `less_lr` and
/// `less_rl` (matching `std::set_difference` in `s_ProcessSeqIdFilters`).
fn set_difference_by<'a, L, R, F, G>(
    lhs: &'a L,
    rhs: &R,
    less_lr: &mut F,
    less_rl: &mut G,
) -> Result<Vec<&'a [u8]>, OidListError>
where
    L: SeqIdSlice + ?Sized,
    R: SeqIdSlice + ?Sized,
    F: FnMut(&[u8], &[u8]) -> Result<bool, OidListError>,
    G: FnMut(&[u8], &[u8]) -> Result<bool, OidListError>,
{
    set_difference_by_get(lhs.len(), |idx| lhs.get(idx), rhs, less_lr, less_rl)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:623-633
fn set_difference_by_get<'a, R, LGet, F, G>(
    lhs_len: usize,
    lhs_get: LGet,
    rhs: &R,
    less_lr: &mut F,
    less_rl: &mut G,
) -> Result<Vec<&'a [u8]>, OidListError>
where
    R: SeqIdSlice + ?Sized,
    LGet: Fn(usize) -> &'a [u8],
    F: FnMut(&[u8], &[u8]) -> Result<bool, OidListError>,
    G: FnMut(&[u8], &[u8]) -> Result<bool, OidListError>,
{
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < lhs_len {
        if j >= rhs.len() {
            out.extend((i..lhs_len).map(&lhs_get));
            break;
        }
        let a = lhs_get(i);
        let b = rhs.get(j);
        if less_lr(a, b)? {
            out.push(a);
            i += 1;
        } else if less_rl(b, a)? {
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:636,707
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:288-333
/// Returns a reference to cached LMDB entries, opening them lazily on first Pass 2 lookup
/// (`AccessionsToOids` or `TaxIdsToOids`).
fn ensure_lmdb_entries<'e>(
    slot: &'e mut Option<Vec<super::filter::LmdbEntry>>,
    lmdb: &dyn Fn() -> Result<Vec<super::filter::LmdbEntry>, SeqDbError>,
) -> Result<&'e [super::filter::LmdbEntry], OidListError> {
    if slot.is_none() {
        *slot = Some(lmdb().map_err(OidListError::reader)?);
    }
    Ok(slot.as_deref().expect("lmdb entries loaded"))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:610-635
/// Validates a sorted alias `SEQIDLIST` against `s_CompareSeqId` and applies the search's
/// `-seqidlist` intersection and/or `-negative_seqidlist` difference, borrowing slices from
/// `sorted_accs`.
fn filter_alias_seqids<'a>(
    sorted_accs: &'a PackedStrings,
    fname: &[u8],
    user_accs: &[Vec<u8>],
    neg_user_accs: &[Vec<u8>],
) -> Result<Vec<&'a [u8]>, OidListError> {
    validate_compare_seqids(sorted_accs, Some(fname), ALIAS_SEQIDLIST_KIND)?;
    let mut matched: Option<Vec<&'a [u8]>> = None;
    if !user_accs.is_empty() {
        let common = set_intersection_by(
            sorted_accs,
            user_accs,
            &mut |a, b| {
                compare_seq_id(
                    a,
                    Some(fname),
                    ALIAS_SEQIDLIST_KIND,
                    b,
                    None,
                    USER_SEQIDLIST_KIND,
                )
            },
            &mut |b, a| {
                compare_seq_id(
                    b,
                    None,
                    USER_SEQIDLIST_KIND,
                    a,
                    Some(fname),
                    ALIAS_SEQIDLIST_KIND,
                )
            },
        )?;
        if common.is_empty() {
            return Ok(common);
        }
        matched = Some(common);
    }
    if !neg_user_accs.is_empty() {
        let mut less_lr = |a: &[u8], b: &[u8]| {
            compare_seq_id(
                a,
                Some(fname),
                ALIAS_SEQIDLIST_KIND,
                b,
                None,
                USER_NEG_SEQIDLIST_KIND,
            )
        };
        let mut less_rl = |b: &[u8], a: &[u8]| {
            compare_seq_id(
                b,
                None,
                USER_NEG_SEQIDLIST_KIND,
                a,
                Some(fname),
                ALIAS_SEQIDLIST_KIND,
            )
        };
        let diff = match &matched {
            Some(prev) => set_difference_by_get(
                prev.len(),
                |idx| prev[idx],
                neg_user_accs,
                &mut less_lr,
                &mut less_rl,
            )?,
            None => set_difference_by(sorted_accs, neg_user_accs, &mut less_lr, &mut less_rl)?,
        };
        return Ok(diff);
    }
    Ok(matched.unwrap_or_default())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:578-654
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1476-1491
/// `s_ProcessSeqIdFilters`: reads each version 5 alias `SEQIDLIST` file, attaches it to its
/// volumes before any intersection/difference with the search's seqid list, resolves remaining
/// accessions through LMDB `AccessionsToOids`, and sets matching OIDs of included volumes in
/// `f_bits`.
#[allow(clippy::too_many_arguments)]
fn process_seqid_filters(
    fnames: &[&[u8]],
    fnames_vols: &[Vec<&[u8]>],
    volumes: &[TreeVolume<'_>],
    num_oids: u32,
    user: Option<&super::filter::DeflineRule>,
    lmdb: &dyn Fn() -> Result<Vec<super::filter::LmdbEntry>, SeqDbError>,
    entries_slot: &mut Option<Vec<super::filter::LmdbEntry>>,
    f_bits: &mut [u64],
    vol_lists: &mut [Vec<std::sync::Arc<AliasIdList>>],
) -> Result<(), OidListError> {
    if fnames.is_empty() {
        return Ok(());
    }
    let mut user_accs = match user {
        Some(super::filter::DeflineRule::Seqids(ids)) if !ids.is_empty() => ids.clone(),
        Some(
            super::filter::DeflineRule::Seqids(_)
            | super::filter::DeflineRule::NegativeSeqids(_)
            | super::filter::DeflineRule::Taxids(_)
            | super::filter::DeflineRule::NegativeTaxids(_)
            | super::filter::DeflineRule::Gis(_)
            | super::filter::DeflineRule::NegativeGis(_)
            | super::filter::DeflineRule::Pigs
            | super::filter::DeflineRule::NegativePigs,
        )
        | None => Vec::new(),
    };
    if !user_accs.is_empty() {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:592
        // Sort bytewise first so validation runs order-independently without extra hash sets;
        // once `validate_compare_seqids` passes, `s_CompareSeqId` reduces to bytewise `<`.
        user_accs.sort_unstable();
        validate_compare_seqids(user_accs.as_slice(), None, USER_SEQIDLIST_KIND)?;
    }
    let neg_user_accs: &[Vec<u8>] = match user {
        Some(super::filter::DeflineRule::NegativeSeqids(ids)) if !ids.is_empty() => ids.as_slice(),
        Some(
            super::filter::DeflineRule::NegativeSeqids(_)
            | super::filter::DeflineRule::Seqids(_)
            | super::filter::DeflineRule::Taxids(_)
            | super::filter::DeflineRule::NegativeTaxids(_)
            | super::filter::DeflineRule::Gis(_)
            | super::filter::DeflineRule::NegativeGis(_)
            | super::filter::DeflineRule::Pigs
            | super::filter::DeflineRule::NegativePigs,
        )
        | None => &[],
    };
    for (fname, vols) in fnames.iter().zip(fnames_vols) {
        let file = {
            let bytes = super::seqidlist::read_memory_file(fname).map_err(OidListError::reader)?;
            let mut warnings = Vec::new();
            let file = super::seqidlist::read_si_list_packed(&bytes, &mut warnings)
                .map_err(|text| OidListError::reader(SeqDbError::NotFound(text.into_bytes())))?;
            if let Some(warning) = warnings.first() {
                return Err(OidListError::reader(SeqDbError::Unsupported {
                    path: path_from_bytes(*fname),
                    message: format!(
                        "the seqid list of an alias file's SEQIDLIST has a token NCBI reports \
                         (`{}`)",
                        warning.trim()
                    ),
                }));
            }
            file
        };
        let is_v4 = file.info.is_v4;
        let mut ids = file.ids;
        // Binary `.bsl` keeps `in_order = true` (stored order); text lists sort on first `FindSi`.
        if is_v4 {
            ids.sort_unstable();
        }
        let has_user_filter = !user_accs.is_empty() || !neg_user_accs.is_empty();
        // When a binary `.bsl` is intersected/differenced with a search list, `attached` keeps
        // stored file order while the sorted view clones only the span index (sharing the bytes).
        let mut sorted_binary = if !is_v4 && has_user_filter && !ids.is_empty() {
            let mut sorted = ids.clone();
            sorted.sort_unstable();
            Some(sorted)
        } else {
            None
        };
        let attached = std::sync::Arc::new(AliasIdList {
            sis: ids,
            ..Default::default()
        });
        let excluded =
            attach_and_excluded(volumes, vols, std::sync::Arc::clone(&attached), vol_lists);
        if attached.sis.is_empty() {
            continue;
        }
        if has_user_filter {
            let sorted_ref = sorted_binary.as_ref().unwrap_or(&attached.sis);
            let matched = filter_alias_seqids(sorted_ref, fname, &user_accs, neg_user_accs)?;
            if matched.is_empty() {
                continue;
            }
            let entries = ensure_lmdb_entries(entries_slot, lmdb)?;
            super::filter::entries_accessions_to_bits(
                entries,
                matched.as_slice(),
                num_oids,
                &excluded,
                f_bits,
            )
            .map_err(OidListError::reader)?;
            drop(matched);
            sorted_binary.take();
        } else {
            let entries = ensure_lmdb_entries(entries_slot, lmdb)?;
            super::filter::entries_accessions_to_bits(
                entries,
                &attached.sis,
                num_oids,
                &excluded,
                f_bits,
            )
            .map_err(OidListError::reader)?;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:656-716
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1080-1121
/// `s_ProcessTaxIdFilters`: reads each version 5 alias `TAXIDLIST` file, attaches it to its
/// volumes before any intersection/difference with the search's taxid list, resolves remaining
/// taxids through LMDB `TaxIdsToOids` (failing with `Taxonomy ID(s) not found...` when no OID
/// matches), and sets matching OIDs of included volumes in `f_bits`.
#[allow(clippy::too_many_arguments)]
fn process_taxid_filters(
    fnames: &[&[u8]],
    fnames_vols: &[Vec<&[u8]>],
    volumes: &[TreeVolume<'_>],
    num_oids: u32,
    user: Option<&super::filter::DeflineRule>,
    lmdb: &dyn Fn() -> Result<Vec<super::filter::LmdbEntry>, SeqDbError>,
    entries_slot: &mut Option<Vec<super::filter::LmdbEntry>>,
    f_bits: &mut [u64],
    vol_lists: &mut [Vec<std::sync::Arc<AliasIdList>>],
) -> Result<(), OidListError> {
    if fnames.is_empty() {
        return Ok(());
    }
    let user_taxids = match user {
        Some(super::filter::DeflineRule::Taxids(found)) if !found.is_empty() => Some(found),
        Some(
            super::filter::DeflineRule::Taxids(_)
            | super::filter::DeflineRule::NegativeTaxids(_)
            | super::filter::DeflineRule::Seqids(_)
            | super::filter::DeflineRule::NegativeSeqids(_)
            | super::filter::DeflineRule::Gis(_)
            | super::filter::DeflineRule::NegativeGis(_)
            | super::filter::DeflineRule::Pigs
            | super::filter::DeflineRule::NegativePigs,
        )
        | None => None,
    };
    let neg_user_taxids = match user {
        Some(super::filter::DeflineRule::NegativeTaxids(found)) if !found.is_empty() => Some(found),
        Some(
            super::filter::DeflineRule::NegativeTaxids(_)
            | super::filter::DeflineRule::Taxids(_)
            | super::filter::DeflineRule::Seqids(_)
            | super::filter::DeflineRule::NegativeSeqids(_)
            | super::filter::DeflineRule::Gis(_)
            | super::filter::DeflineRule::NegativeGis(_)
            | super::filter::DeflineRule::Pigs
            | super::filter::DeflineRule::NegativePigs,
        )
        | None => None,
    };
    for (fname, vols) in fnames.iter().zip(fnames_vols) {
        let bytes = super::seqidlist::read_memory_file(fname).map_err(OidListError::reader)?;
        let raw = super::id_list::read_tax_id_list(&bytes)
            .map_err(|text| OidListError::reader(SeqDbError::NotFound(text.into_bytes())))?;
        let raw_non_empty = !raw.is_empty();
        let mut has_past = false;
        let mut taxids = std::collections::BTreeSet::new();
        for number in raw {
            match number {
                super::id_list::ListNumber::Value(value) => {
                    taxids.insert(value as i32);
                }
                super::id_list::ListNumber::Past => has_past = true,
            }
        }
        let attached = std::sync::Arc::new(AliasIdList {
            taxids: taxids.clone(),
            ..Default::default()
        });
        let excluded = attach_and_excluded(volumes, vols, attached, vol_lists);
        if !raw_non_empty {
            continue;
        }
        if let Some(user_set) = user_taxids {
            taxids.retain(|taxid| user_set.contains(taxid));
            has_past = false;
            if taxids.is_empty() {
                continue;
            }
        }
        if let Some(neg_set) = neg_user_taxids {
            taxids.retain(|taxid| !neg_set.contains(taxid));
            if taxids.is_empty() && !has_past {
                continue;
            }
        }
        let entries = ensure_lmdb_entries(entries_slot, lmdb)?;
        let (oids, _) = super::filter::entries_tax_oids(entries, false, &taxids)
            .map_err(OidListError::reader)?;
        if oids.is_empty() {
            return Err(OidListError::reader(SeqDbError::NotFound(
                super::filter::tax_ids_not_found(entries),
            )));
        }
        for oid in oids {
            if oid >= num_oids || excluded.iter().any(|&(s, e)| oid >= s && oid < e) {
                continue;
            }
            f_bits[(oid / 64) as usize] |= 1u64 << (oid % 64);
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:101-107
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:718-761
/// Version 5 Pass 2 (`CSeqDBOIDList::x_ComputeFilters` overload for `CSeqDBLMDBSet`): scans only
/// the root filters (`v5_root_lists`) of each volume's specialized tree for `SEQIDLIST` and
/// `TAXIDLIST`, clears `f_bits` for those volumes, unions the OIDs resolved by
/// `process_seqid_filters` and `process_taxid_filters`, and intersects `words` with `f_bits`.
/// Returns `true` when at least one v5 `SEQIDLIST` or `TAXIDLIST` was applied (making NCBI's
/// `m_AllBits` a bit array, not `Special::AllSet`).
fn v5_list_filters<'a>(
    volumes: &[TreeVolume<'a>],
    v5_root_lists: &[&'a [AliasList]],
    num_oids: u32,
    user: Option<&super::filter::DeflineRule>,
    lmdb: &dyn Fn() -> Result<Vec<super::filter::LmdbEntry>, SeqDbError>,
    words: &mut [u64],
    vol_lists: &mut [Vec<std::sync::Arc<AliasIdList>>],
) -> Result<bool, OidListError> {
    let mut seqid_fnames: Vec<&'a [u8]> = Vec::new();
    let mut seqid_vols: Vec<Vec<&'a [u8]>> = Vec::new();
    let mut taxid_fnames: Vec<&'a [u8]> = Vec::new();
    let mut taxid_vols: Vec<Vec<&'a [u8]>> = Vec::new();
    let mut cleared: Vec<(u64, u64)> = Vec::new();
    for (volume, root_lists) in volumes.iter().zip(v5_root_lists) {
        let mut vol_cleared = false;
        for list in *root_lists {
            match list.kind {
                IdListKind::Si => {
                    add_filter_file(&list.path, volume.name, &mut seqid_fnames, &mut seqid_vols);
                    vol_cleared = true;
                }
                IdListKind::TaxId => {
                    add_filter_file(&list.path, volume.name, &mut taxid_fnames, &mut taxid_vols);
                    vol_cleared = true;
                }
                IdListKind::Gi | IdListKind::Ti => {}
            }
        }
        if vol_cleared {
            cleared.push((u64::from(volume.start), u64::from(volume.end)));
        }
    }
    if seqid_fnames.is_empty() && taxid_fnames.is_empty() {
        return Ok(false);
    }
    let mut f_bits = vec![0u64; (num_oids as usize).div_ceil(64)];
    set_bits(&mut f_bits, 0, u64::from(num_oids));
    for (start, end) in cleared {
        clear_bits(&mut f_bits, start, end);
    }
    let mut entries_slot = None;
    process_seqid_filters(
        &seqid_fnames,
        &seqid_vols,
        volumes,
        num_oids,
        user,
        lmdb,
        &mut entries_slot,
        &mut f_bits,
        vol_lists,
    )?;
    process_taxid_filters(
        &taxid_fnames,
        &taxid_vols,
        volumes,
        num_oids,
        user,
        lmdb,
        &mut entries_slot,
        &mut f_bits,
        vol_lists,
    )?;
    for (word, kept) in words.iter_mut().zip(&f_bits) {
        *word &= kept;
    }
    Ok(true)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.hpp:104-114
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.hpp:248-251
/// The shape of an NCBI `CSeqDB_BitSet`, which its operations decide by: the OIDs
/// `[start, end)` and whether it is a special set (all of them set, `eAllSet`, or none,
/// `eAllClear`) or a bit array (`eNone`) of one bit per OID, rounded up to whole bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Shape {
    start: u64,
    end: u64,
    special: Special,
}

/// `CSeqDB_BitSet::ESpecialCase`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Special {
    /// `eNone`: a bit array.
    None,
    /// `eAllSet`.
    AllSet,
    /// `eAllClear`.
    AllClear,
}

/// The result of `UnionWith` ([`Shape::union_with`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Union {
    /// The set itself.
    Left,
    /// A copy of the other set.
    Right,
    /// A bit array of this shape (`x_Normalize` widens the range to both).
    Array(Shape),
}

impl Shape {
    fn new(start: u64, end: u64, special: Special) -> Self {
        Self {
            start,
            end,
            special,
        }
    }

    /// The bytes of the bit array (`x_Alloc`).
    fn bytes(self) -> u64 {
        (self.end - self.start).div_ceil(8)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:142-189
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:309-334
    /// `UnionWith`: which set the union is, by the shapes.
    fn union_with(self, other: Self) -> Union {
        if other.special == Special::AllClear {
            return Union::Left;
        }
        if self.special == Special::AllClear {
            return Union::Right;
        }
        if other.start >= self.start && other.end <= self.end && self.special == Special::AllSet {
            return Union::Left;
        }
        if other.start <= self.start && other.end >= self.end && other.special == Special::AllSet {
            return Union::Right;
        }
        Union::Array(Self::new(
            self.start.min(other.start),
            self.end.max(other.end),
            Special::None,
        ))
    }
}

/// An NCBI `CSeqDB_BitSet` of OIDs of one volume ([`volume_set`]): its shape and, for a bit
/// array, its OIDs.
#[derive(Debug, Clone)]
struct OidSet {
    shape: Shape,
    /// The OIDs of a bit array (`Special::None`), bit `oid - frame.start` of the volume's
    /// [`Frame`]; empty for a special set, whose OIDs its shape gives.
    bits: Vec<u64>,
}

impl OidSet {
    fn special(start: u64, end: u64, special: Special) -> Self {
        Self {
            shape: Shape::new(start, end, special),
            bits: Vec::new(),
        }
    }

    /// A bit array over the OIDs of the volume (`x_GetOidMask`).
    fn array(frame: Frame, bits: Vec<u64>) -> Self {
        Self {
            shape: Shape::new(frame.start, frame.end, Special::None),
            bits,
        }
    }

    /// The OIDs of the set, as bits over the volume's OIDs.
    fn into_words(self, frame: Frame) -> Vec<u64> {
        match self.shape.special {
            Special::None => self.bits,
            Special::AllClear => frame.clear_words(),
            Special::AllSet => frame.range_words(self.shape.start, self.shape.end),
        }
    }

    /// ORs the OIDs of the set into the global bits `words`.
    fn or_into(self, words: &mut [u64], frame: Frame) {
        match self.shape.special {
            Special::AllClear => {}
            Special::AllSet => set_bits(
                words,
                self.shape.start.max(frame.start),
                self.shape.end.min(frame.end),
            ),
            Special::None => or_at(words, &self.bits, frame.start as u32),
        }
    }

    /// `UnionWith`: the union of the OIDs, in the shape [`Shape::union_with`] gives.
    fn union_with(self, other: Self, frame: Frame) -> Self {
        match self.shape.union_with(other.shape) {
            Union::Left => self,
            Union::Right => other,
            Union::Array(shape) => {
                let mut bits = self.into_words(frame);
                for (word, other) in bits.iter_mut().zip(other.into_words(frame)) {
                    *word |= other;
                }
                Self { shape, bits }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:222-243
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:376-405
    // Of a bit array and a special all-set, NCBI clears the array's bits from its first OID up
    // to the range's start and from the range's end up to its last OID, one bit or one byte at a
    // time at the OID less the array's first. A range that ends before the array's first OID
    // makes the first index wrap below 0: NCBI writes far outside the array and crashes (probes
    // D1a_P, D1a_R of line RG D1b-1, SIGSEGV); one that starts past the array's bytes writes past
    // them, over NCBI's heap: NCBI aborts or crashes (probes RB9_T200 to RB9_T19990 of review B of
    // D1b-1), or, where the cleared bytes stay in the unused end of the array's allocation, keeps
    // no OID (RB_S, RB9_T21), as where it reads past an OID list file's page
    // (`OidListError::PastPage`). Inside the array's last byte the clearing is harmless (probe
    // D1a_Q: no OID).
    /// Whether NCBI's `IntersectWith` of this set with `other` clears bits outside this bit
    /// array (`OidListError::RangeOutsideBits`).
    fn clears_outside(&self, other: &Self) -> bool {
        let (a, b) = (self.shape, other.shape);
        if a.special != Special::None || b.special != Special::AllSet {
            return false;
        }
        b.end < a.start || (b.start > a.end && (b.start - 1 - a.start) / 8 >= a.bytes())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:191-294
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:90-140
    // Of a special all-set and a bit array, the bit array keeps its shape and loses its bits
    // outside the all-set's range (where NCBI's clearing stays inside the array:
    // `clears_outside`). Two bit arrays with the same first OID and the same bytes
    // are ANDed byte by byte. Any other two bit arrays take NCBI's last loop, which walks the set
    // bits of this set and passes its index by reference to the other set's search for its next
    // set bit: where the other set has one, the index jumps there and no bit is cleared, so the
    // only bits cleared are those past the other set's last set bit (probes B1 of review A of line
    // RG D1b-1: OIDs 2-4 and 8-10 under an OID list of OID 2 keep OIDs 1 and 2, under one of OID 8
    // the OIDs 1-3, 7, 8). The result keeps this set's shape. That loop is modeled for an `other`
    // that starts at or before this set (the filter over the volume's OIDs of `volume_set`): of an
    // `other` without set bits that starts later, NCBI keeps this set's bits before it.
    /// `IntersectWith`: this set ANDed with `other` as NCBI computes it.
    fn intersect_with(self, other: Self, frame: Frame) -> Self {
        let (a, b) = (self.shape, other.shape);
        if a.special == Special::AllClear {
            return self;
        }
        if b.special == Special::AllClear {
            return other;
        }
        match (a.special, b.special) {
            (Special::AllSet, Special::AllSet) => {
                let start = a.start.max(b.start);
                let end = a.end.min(b.end);
                if start >= end {
                    Self::special(a.start, a.end, Special::AllClear)
                } else {
                    Self::special(start, end, Special::AllSet)
                }
            }
            (Special::AllSet, _) => other.and_words(frame.range_words(a.start, a.end)),
            (_, Special::AllSet) => self.and_words(frame.range_words(b.start, b.end)),
            _ if a.start == b.start && a.bytes() == b.bytes() => self.and_words(other.bits),
            _ => {
                debug_assert!(b.start <= a.start, "the other set starts after this one");
                let last = other
                    .bits
                    .iter()
                    .rposition(|&word| word != 0)
                    .map(|at| at as u64 * 64 + 63 - u64::from(other.bits[at].leading_zeros()));
                let kept = match last {
                    Some(last) => frame.range_words(frame.start, frame.start + last + 1),
                    None => frame.clear_words(),
                };
                self.and_words(kept)
            }
        }
    }

    /// The bit array with only its OIDs that `kept` (bits over the volume's OIDs) holds.
    fn and_words(mut self, kept: Vec<u64>) -> Self {
        for (word, kept) in self.bits.iter_mut().zip(kept) {
            *word &= kept;
        }
        self
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfilter.cpp:57-78
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfilter.cpp:38-55
// The filters of a node are its OID range, its OID list, its id lists and its membership bit.
/// A node of NCBI's filter tree specialised to one volume (`CSeqDB_FilterTree::Specialize`):
/// its OID range, OID list, id lists and membership bit, its sub-nodes, and how many times it
/// names the volume.
struct VolumeNode<'a> {
    /// The alias file of the node (`m_ThisName`), `None` for the top node.
    name: Option<&'a [u8]>,
    range: Option<OidRange>,
    oid_list: Option<&'a [u8]>,
    /// The `OID_MASK_TYPE` of `oid_list`.
    oid_mask_type: u32,
    lists: &'a [AliasList],
    memb_bit: Option<u32>,
    nodes: Vec<VolumeNode<'a>>,
    volumes: usize,
}

impl VolumeNode<'_> {
    /// Whether the node has a filter of its own (`GetFilters()` not empty).
    fn has_filters(&self) -> bool {
        self.range.is_some()
            || self.oid_list.is_some()
            || !self.lists.is_empty()
            || self.memb_bit.is_some()
    }
}

/// `Specialize` for `volume`, of the tree under `node`; `None` when no branch reaches the volume.
fn specialize<'a>(node: &'a AliasNode, volume: &[u8]) -> Option<VolumeNode<'a>> {
    let mut clone = VolumeNode {
        name: node.this_name.as_deref(),
        range: node.range,
        oid_list: node.oid_list.as_deref(),
        oid_mask_type: node.oid_mask_type,
        lists: &node.lists,
        memb_bit: node.memb_bit,
        nodes: Vec::new(),
        volumes: 0,
    };
    for sub in node
        .subnodes
        .iter()
        .filter_map(|sub| specialize(sub, volume))
    {
        if sub.has_filters() {
            clone.nodes.push(sub);
        } else {
            clone.nodes.extend(sub.nodes);
            clone.volumes += sub.volumes;
        }
    }
    clone.volumes += node.vol_names.iter().filter(|name| *name == volume).count();
    while !clone.has_filters() && clone.volumes == 0 && clone.nodes.len() == 1 {
        clone = clone.nodes.pop().expect("one sub-node");
    }
    (!clone.nodes.is_empty() || clone.volumes > 0).then_some(clone)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:247-258
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:183-201
// The filter tree's top node has one sub-node per alias file of the top node
// (`BuildFilterTree`, seqdbalias.cpp:1880-1894); a node with an OID list or an id list besides
// its range has two filters or more. The bounds are NCBI's `int`s
// (`OidRange::bounds`), negative
// for a value of 2^31 or more: a negative `oid_begin` starts at 0, a negative `oid_end` ends at
// 0, and `LAST_OID 0` (`oid_end` 0) ends at the last OID.
/// The iteration range NCBI sets on its first OID-list build when the top node has exactly one
/// alias file and that file's only filter is an OID range: `[begin, end)` as `SetIterationRange` clamps it
/// to the `num_oids` OIDs of the volume set. `None` in every other case.
pub(super) fn top_iteration_range(root: &AliasNode, num_oids: u32) -> Option<(u32, u32)> {
    let [only] = root.subnodes.as_slice() else {
        return None;
    };
    if only.oid_list.is_some() || !only.lists.is_empty() || only.memb_bit.is_some() {
        return None;
    }
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
// blastdbcmd prints this list under "Volumes:" (`FindVolumePaths` with `volumes` and
// false, blastdbcmd.cpp:568): the volumes of the first two levels of the
// alias tree and the alias files of the third. The tests compare it with what
// NCBI printed.
/// What `blastdbcmd -info` lists under "Volumes:" for the database `dbname`.
pub(super) fn info_volumes(
    dbname: &[u8],
    seqtype: SeqType,
    search: &SearchPath,
) -> Result<Vec<Vec<u8>>, SeqDbError> {
    let tree = resolve_in(dbname, seqtype, search)?;
    let mut vols: Vec<Vec<u8>> = tree.root.vol_names.clone();
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

#[cfg(test)]
#[path = "alias_tests.rs"]
mod tests;
