//! `losat blastdbcmd -list DIR`: the databases of a directory (`FindBlastDBs`) written with
//! `-list_outfmt` (`CBlastDbFormatter`).

use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

use super::args::Args;
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

fn bytes_of(path: &Path) -> Vec<u8> {
    path.as_os_str().as_encoded_bytes().to_vec()
}

fn path_of(bytes: &[u8]) -> std::path::PathBuf {
    std::path::PathBuf::from(String::from_utf8_lossy(bytes).into_owned())
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:3019-3066
// ```c++
// template<class TFindFunc>
// void FindFilesInDir(const CDir&            dir,
//                     const vector<string>&  masks,
//                     const vector<string>&  masks_subdir,
//                     TFindFunc&             find_func,
//                     TFindFiles             flags = fFF_Default)
// {
//     TFindFiles find_type = flags & fFF_All;
//     ...
//     unique_ptr<CDir::TEntries>
//         contents(dir.GetEntriesPtr(kEmptyStr, CDir::fIgnoreRecursive | CDir::fIgnorePath));
//     if (contents.get() == NULL) {
//         // error
//         return;
//     }
//     ...
//     if ( dir.GetPath().length() ) {
//         path = CDirEntry::AddTrailingPathSeparator(dir.GetPath());
//     }
//     ITERATE(CDir::TEntries, it, *contents) {
//         CDirEntry& dir_entry = **it;
//         string name = dir_entry.GetPath();
//         dir_entry.Reset(CDirEntry::MakePath(path, name));
//
//         TFindFiles entry_type = fFF_Dir | fFF_File; // unknown
//         if ( CDirEntry::MatchesMask(name, masks, use_case) ) {
//             if ( find_type != (fFF_Dir | fFF_File) ) {
//                 // need to check actual entry type
//                 entry_type = dir_entry.IsDir()? fFF_Dir: fFF_File;
//             }
//             if ( (entry_type & find_type) != 0 ) {
//                 // entry type matches
//                 find_func(dir_entry);
//             }
//         }
//         if ( (flags & fFF_Recursive) &&
//              (entry_type & fFF_Dir) /*possible dir*/ &&
//              CDirEntry::MatchesMask(name, masks_subdir, use_case) &&
//              (entry_type == fFF_Dir || dir_entry.IsDir()) /*real dir*/ ) {
//             CDir nested_dir(dir_entry.GetPath());
//             FindFilesInDir(nested_dir, masks, masks_subdir, find_func, flags);
//         }
//     }
//     return;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1350-1364
// ```c++
// class CBlastDbFinder {
// public:
//     void operator() (CDirEntry& de) {
//         const string& extn = de.GetPath().substr(de.GetPath().length() - 3, 1);
//         SSeqDBInitInfo value;
//         // rm extension
//         value.m_BlastDbName = de.GetPath().substr(0, de.GetPath().length() - 4);
//         CNcbiOstrstream oss;
//         // Needed for escaping spaces
//         oss << "\"" << value.m_BlastDbName << "\"";
//         value.m_BlastDbName = CNcbiOstrstreamToString(oss);
//         value.m_MoleculeType =
//             (extn == "n" ? CSeqDB::eNucleotide : CSeqDB::eProtein);
//         m_DBs.push_back(value);
//     }
// ```
// The masks are `*.pin`, `*.pal`, `*.nin`, `*.nal` (with case); a file is not opened, so an
// empty or stray `x.pin` is a database, and so is a link that leads nowhere (`IsDir` is
// false for it). The sub-directories (`*`) are entered with `-recursive`, links to
// directories too. The order of the directory does not matter: the list is sorted.
/// `FindFilesInDir(CDir(dir), masks, {"*"}, finder, fFF_File | recursive)`: the databases of
/// the files of `dir` whose names end with one of `extensions` (`.pin` ...).
fn find_files(dir: &[u8], extensions: &[&str], recursive: bool, found: &mut Vec<FoundDb>) {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3879
    // ```c++
    //     string base_path =  AddTrailingPathSeparator(GetPath().empty() ? DIR_CURRENT : GetPath());
    // ```
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
// ```c++
//     string GetFileName(size_t idx) {
//         SSeqDBInitInfo& info = m_DBs[idx];
//         string retval = NStr::Replace(info.m_BlastDbName, "\"", kEmptyStr);
//         if (info.m_MoleculeType == CSeqDB::eNucleotide) {
//             string alias = retval + ".nal", index = retval + ".nin";
//             retval = (CFile(alias).Exists() ? alias : index);
//         } else {
//             string alias = retval + ".pal", index = retval + ".pin";
//             retval = (CFile(alias).Exists() ? alias : index);
//         }
//         return retval;
//     }
// ...
//     bool operator() (const SSeqDBInitInfo& value) const {
//         return (NStr::Find(value.m_BlastDbName, m_Path) != NPOS);
//     }
// ...
// static void s_RemoveAliasComponents(CBlastDbFinder& finder)
// {
//     set<string> dbs2remove;
//     for (size_t i = 0; i < finder.m_DBs.size(); i++) {
//         string path = finder.GetFileName(i);
//         if (path[path.size()-1] != 'l') { // not an alias file
//             continue;
//         }
//         CNcbiIfstream in(path.c_str());
//         if (!in) {
//             continue;
//         }
//         string line;
//         while (getline(in, line)) {
//             if (NStr::StartsWith(line, "DBLIST")) {
//                 vector<string> tokens;
//                 NStr::Split(line, " ", tokens, NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
//                 for (size_t j = 1; j < tokens.size(); j++) {
//                     dbs2remove.insert(tokens[j]);
//                 }
//             }
//         }
//     }
//
//     ITERATE(set<string>, i, dbs2remove) {
//         finder.m_DBs.erase(remove_if(finder.m_DBs.begin(), finder.m_DBs.end(),
//                                      PathFinder(*i)),
//                            finder.m_DBs.end());
//     }
// }
// ```
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1429-1462
// ```c++
// vector<SSeqDBInitInfo>
// FindBlastDBs(const string& path, const string& dbtype, bool recurse,
//              bool include_alias_files /* = false */,
//              bool remove_redundant_dbs /* = false */)
// {
//     // 1. Find every database volume (but not alias files etc).
//     vector<string> fmasks, dmasks;
//
//     // If the type is 'guess' we do both types of databases.
//
//     if (dbtype != "nucl") {
//         fmasks.push_back("*.pin");
//         if (include_alias_files) {
//             fmasks.push_back("*.pal");
//         }
//     }
//     if (dbtype != "prot") {
//         fmasks.push_back("*.nin");
//         if (include_alias_files) {
//             fmasks.push_back("*.nal");
//         }
//     }
//     dmasks.push_back("*");
//
//     EFindFiles flags = (EFindFiles)
//         (fFF_File | (recurse ? fFF_Recursive : 0));
//
//     CBlastDbFinder dbfinder;
//     FindFilesInDir(CDir(path), fmasks, dmasks, dbfinder, flags);
//     if (remove_redundant_dbs) {
//         s_RemoveAliasComponents(dbfinder);
//     }
//     sort(dbfinder.m_DBs.begin(), dbfinder.m_DBs.end());
//     return dbfinder.m_DBs;
// }
// ```
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdb.hpp:1550-1561
// ```c++
//     inline bool operator<(const SSeqDBInitInfo& rhs) const {
//         if (m_BlastDbName < rhs.m_BlastDbName) {
//             return true;
//         } else if (m_BlastDbName > rhs.m_BlastDbName) {
//             return false;
//         } else {
//             return ((int)m_MoleculeType < (int)rhs.m_MoleculeType);
//         }
//     }
// ```
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
// ```c++
// CBlastDbFormatter::CBlastDbFormatter(const string& fmt_spec)
//     : m_FmtSpec(fmt_spec)
// {
//     // Record where the offsets where the replacements must occur
//     for (SIZE_TYPE i = 0; i < m_FmtSpec.size(); i++) {
//         if (m_FmtSpec[i] == '%' && m_FmtSpec[i+1] == '%') {
//             // remove the escape character for '%'
//             m_FmtSpec.erase(i++, 1);
//             continue;
//         }
//
//         if (m_FmtSpec[i] == '%') {
//             m_ReplOffsets.push_back(i);
//             m_ReplacementTypes.push_back(m_FmtSpec[i+1]);
//         }
//     }
//     // Handle %d defline in ASN.1 text format, can only be by itself
//
//     if (m_ReplOffsets.empty() ||
//         m_ReplacementTypes.size() != m_ReplOffsets.size()) {
//         NCBI_THROW(CInvalidDataException, eInvalidInput,
//                    "Invalid format specification");
//     }
// }
// ```
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
    // ```c++
    // string
    // CBlastDbFormatter::x_Replacer(const vector<string>& data2write) const
    // {
    //     ...
    //     SIZE_TYPE fmt_idx = 0;
    //     for (SIZE_TYPE i = 0, kSize = m_ReplOffsets.size(); i < kSize; i++) {
    //         retval.append(&m_FmtSpec[fmt_idx], &m_FmtSpec[m_ReplOffsets[i]]);
    //         retval.append(data2write[i]);
    //         fmt_idx = m_ReplOffsets[i] + 2;
    //     }
    //     if (fmt_idx <= m_FmtSpec.size()) {
    //         retval.append(&m_FmtSpec[fmt_idx], &m_FmtSpec[m_FmtSpec.size()]);
    //     }
    //
    //     return retval;
    // }
    // ```
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
    fn from(error: SeqDbError) -> Self {
        match error.ncbi_message() {
            Some(text) => ListStop::Ncbi(format!("BLAST Database error: {text}").into_bytes(), 2),
            None => ListStop::Refused(refusal(&error.to_string())),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2620-2654
// ```c++
// void SeqDB_GetFileExtensions(bool db_is_protein, vector<string>& extn, EBlastDbVersion dbver)
// {
//     ...
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
//     if (dbver == eBDB_Version5) {
//     	vector<string> lmdbs;
//     	SeqDB_GetLMDBFileExtensions(db_is_protein, lmdbs);
//     	extn.insert(extn.end(), lmdbs.begin(), lmdbs.end());
//     }
//     // Contain masking information
//     extn.push_back(kExtnMol + "aa");   // ISAM mask index file
//     ...
//     extn.push_back(kExtnMol + "xm");   // OID mask file to exclude models
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2657-2666
// ```c++
// void SeqDB_GetLMDBFileExtensions(bool db_is_protein, vector<string>& extn)
// {
//
// 	static const char * ext[]={"db", "os", "ot", "tf", "to", "db-lock", "tf-lock", NULL};
// ```
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
// ```c++
// Int8 CSeqDB::GetDiskUsage() const
// {
//     vector<string> paths;
//     FindVolumePaths(paths);
//     ...
//     ITERATE(vector<string>, path, paths) {
//         ITERATE(vector<string>, ext, extn) {
//             CFile file(*path + "." + *ext);
//             if (file.Exists()) {
//                 Int8 length = file.GetLength();
//                 if (length != -1) {
//                     retval += length;
//                     ...
//                     blastdb_dirname = file.GetDir();
//     ...
//     // For multi-volume databases, take into account files that apply to the
//     // entire BLASTDB
//     if (paths.size() > 1) {
//         _ASSERT( !blastdb_dirname.empty() );
//         auto dbname = GetDBNameList();
//         vector<string> dblist;
//         NStr::Split(dbname, " ", dblist, NStr::fSplit_Tokenize);
//         if (dblist.size() > 1) {
//             CNcbiOstrstream oss;
//             oss << "Cannot compute disk usage for multiple BLASTDBs (i.e.: '"
//                 << dbname << "') at once. Please try again using one BLASTDB "
//                 << "at a time.";
//             NCBI_THROW(CSeqDBException, eArgErr, CNcbiOstrstreamToString(oss));
//         }
//
//         for (const auto& ext: extn) {
//             CFile file(CDirEntry::MakePath(blastdb_dirname, dbname, ext));
//             if (file.Exists()) {
//     ...
// ```
// The name list of a listed database is its name in double quotes, so the files of the
// whole database (`.nal`, `.ndb`, ...) are looked for as `<dir>/"<name>".<ext>` and not
// found: `%U` of a database of several volumes is the sum of its volumes' files (probe L02,
// J161 of line BDC). The volumes are those of the whole alias tree, each once.
/// `GetDiskUsage` of the listed database `db` (of version 5 when `version5`).
fn disk_usage(db: &FoundDb, version5: bool, path: &SearchPath) -> Result<u64, ListStop> {
    let name = String::from_utf8_lossy(&db.quoted).into_owned();
    let (volumes, _) = SeqDb::volume_and_alias_paths(&name, db.seqtype, path)?;
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_formatter.cpp:68-160
// ```c++
// class CBlastDbMetadata {
//     ...
//     string GetFileName() const {
//         return NStr::Replace(m_DbInitInfo.m_BlastDbName, "\"", kEmptyStr);
//     }
//     ...
//     string GetNumberOfSequences() {
//         x_InitBlastDb();
//         // FIXME: should this use CSeqDB::GetTotals?
//         return NStr::IntToString(m_BlastDb->GetNumSeqs());
//     }
//     ...
//     void x_InitBlastDb() {
//         if (m_BlastDb.Empty()) {
//             m_BlastDb = m_DbInitInfo.InitSeqDb();
//         }
//     ...
// string
// CBlastDbFormatter::Write(const SSeqDBInitInfo& db_init_info)
// {
//     CBlastDbMetadata dbmeta(db_init_info);
//     vector<string> data2write;
//     data2write.reserve(m_ReplacementTypes.size());
//     ITERATE(vector<char>, fmt, m_ReplacementTypes) {
//         switch (*fmt) {
//         case 'f':   // file name
//             data2write.push_back(dbmeta.GetFileName());
//             break;
//         ...
//         default:
//             CNcbiOstrstream os;
//             os << "Unrecognized format specification: '%" << *fmt << "'";
//             NCBI_THROW(CInvalidDataException, eInvalidInput,
//                        CNcbiOstrstreamToString(os));
//         }
//     }
//     return x_Replacer(data2write);
// }
// ```
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
        let name = String::from_utf8_lossy(&db.quoted).into_owned();
        *opened = Some(SeqDb::open_in_with_masks(&name, db.seqtype, path)?);
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
                .num_seqs()
                .to_string()
                .into_bytes(),
            b'l' => open(&mut opened, db, path)?
                .total_length()
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
// ```c++
//         } else if (args["list"]) {
//             const string& blastdb_dir = args["list"].AsString();
//             const bool recurse = args["recursive"];
//             const bool remove_redundant_dbs = args["remove_redundant_dbs"];
//             const string dbtype = args[kArgDbType]
//                 ? args[kArgDbType].AsString()
//                 : "guess";
//             const string& kOutFmt = args["list_outfmt"].AsString();
//             const vector<SSeqDBInitInfo> dbs =
//                 FindBlastDBs(blastdb_dir, dbtype, recurse, true,
//                              remove_redundant_dbs);
//             CBlastDbFormatter blastdb_fmt(kOutFmt);
//             ITERATE(vector<SSeqDBInitInfo>, db, dbs) {
//                 out << blastdb_fmt.Write(*db) << NcbiEndl;
//             }
//             return status;
//         }
// ```
// The format is read after the walk, so a format without `%` is an error also for a
// directory without databases (probe J173). NCBI writes each line as it goes; an error of
// NCBI's comes after the lines before it (`CInvalidDataException`: `BLAST options error:`,
// exit status 1; `CSeqDBException`: `BLAST Database error:`, exit status 2), while LOSAT's
// refusal of a database it does not read (an empty `x.pin` on which NCBI crashes: J168)
// writes no line.
/// `-list`: the databases of the directory, one line each, to `out`; the exit status.
pub(super) fn list(args: &Args, path: &SearchPath, diag: &mut Diag, out: &mut dyn Write) -> i32 {
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
        // ```c++
        //     catch (const CInvalidDataException& e) {                                \
        //         LOG_POST(Error << "BLAST options error: " << e.GetMsg());           \
        //         exit_code = BLAST_INPUT_ERROR;                                      \
        //     }                                                                       \
        // ```
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
                diag.post_bytes(&message);
                return status;
            }
            Err(ListStop::Refused(reason)) => return diag.refuse(&reason),
        }
    }
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does this.
    let _ = out.write_all(&text);
    let _ = out.flush();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn found(name: &str, seqtype: SeqType) -> FoundDb {
        FoundDb {
            quoted: format!("\"{name}\"").into_bytes(),
            seqtype,
        }
    }

    // Probes L01-L05, J160-J178 of line BDC.
    #[test]
    fn the_format_is_read_as_ncbi_s_loop_reads_it() {
        let format = ListFormat::parse(b"%f %p").unwrap();
        assert_eq!(format.offsets, [0, 3]);
        assert_eq!(format.types, b"fp");
        let format = ListFormat::parse(b"%%f %f %%").unwrap();
        assert_eq!(format.spec, b"%f %f %");
        assert_eq!(format.offsets, [3]);
        assert_eq!(format.replace(&[b"db/x".to_vec()]), b"%f db/x %");
        let format = ListFormat::parse(b"%f %").unwrap();
        assert_eq!(format.types, b"f\0");
        assert_eq!(format.replace(&[b"a".to_vec(), b"b".to_vec()]), b"a b");
        // After `%%` the loop steps over the byte that follows the second `%`.
        assert!(ListFormat::parse(b"%%%f").is_err());
        assert!(ListFormat::parse(b"abc").is_err());
        assert!(ListFormat::parse(b"").is_err());
        assert!(ListFormat::parse(b"%%").is_err());
    }

    #[test]
    fn the_list_sorts_on_the_quoted_names() {
        let mut dbs = vec![
            found("d/a!", SeqType::Nucleotide),
            found("d/zz", SeqType::Protein),
            found("d/x.00", SeqType::Protein),
            found("d/a b", SeqType::Protein),
            found("d/x", SeqType::Protein),
            found("d/zz", SeqType::Nucleotide),
            found("d/.h", SeqType::Protein),
        ];
        let rank = |seqtype: SeqType| matches!(seqtype, SeqType::Nucleotide) as u8;
        dbs.sort_by(|a, b| {
            a.quoted
                .cmp(&b.quoted)
                .then(rank(a.seqtype).cmp(&rank(b.seqtype)))
        });
        let names: Vec<String> = dbs
            .iter()
            .map(|db| String::from_utf8(db.file_name()).unwrap())
            .collect();
        assert_eq!(
            names,
            ["d/.h", "d/a b", "d/a!", "d/x", "d/x.00", "d/zz", "d/zz"]
        );
        assert_eq!(dbs[5].seqtype, SeqType::Protein);
    }

    #[test]
    fn a_dblist_name_removes_every_name_that_holds_it() {
        assert!(contains(b"\"db/p5\"", b"\"db/p5\""));
        assert!(!contains(b"\"db/p5n\"", b"\"db/p5\""));
        assert!(!contains(b"\"db/j/lst/p5\"", b"\"p5\""));
        assert!(contains(b"\"db/n55\"", b"n5"));
        assert!(contains(b"\"db/s18mv.00\"", b"s18mv.00"));
    }

    #[test]
    fn the_extensions_follow_the_version() {
        assert_eq!(file_extensions(false).len(), 20);
        let v5 = file_extensions(true);
        assert_eq!(v5.len(), 25);
        assert!(!v5.contains(&"si"));
        assert!(v5.contains(&"tf-lock"));
    }
}
