//! Alias files (`.pal` / `.nal`) and the list of volumes behind a database name.
//!
//! LOSAT reads the two keys a multi-volume database of makeblastdb uses,
//! `TITLE` and `DBLIST`. The keys that filter sequences or replace the
//! statistics are explicit errors (not supported yet); unknown keys are
//! ignored, as in NCBI.

use std::path::{Path, MAIN_SEPARATOR};

use super::volume::with_extension;
use super::{SeqDbError, SeqType};

/// One node of the alias tree: the database name given by the caller (the
/// top node) or one alias file.
#[derive(Debug, Default)]
pub(super) struct AliasNode {
    /// `TITLE` of the alias file, when it has that key.
    pub title: Option<Vec<u8>>,
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
/// statistics ("meta-data"), or name data LOSAT does not read. They are the
/// keys NCBI's reader looks up (seqdbalias.cpp), except `MAX_SEQ_LENGTH`,
/// which it never applies (seqdbimpl.cpp:972-976 takes the maximum from the
/// volumes) and which is ignored here too, plus `VDBLIST` of
/// blastdb_aliastool.
const UNSUPPORTED_KEYS: [&str; 16] = [
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
    "NSEQ",
    "LENGTH",
    "STATS_NSEQ",
    "STATS_TOTLEN",
    "MIN_SEQ_LENGTH",
    "OID_MASK_TYPE",
];

fn exists(path: &Path) -> bool {
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
    let name = base_name(path);
    let bytes = name.as_bytes();
    let len = bytes.len();
    if len > 4
        && bytes[len - 4] == b'.'
        && (bytes[len - 3] == b'n' || bytes[len - 3] == b'p')
        && matches!(&bytes[len - 2..], b"al" | b"in" | b"db")
    {
        &name[..len - 4]
    } else {
        name
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

struct Resolver {
    seqtype: SeqType,
    /// Alias files being expanded, outermost first.
    recurse: Vec<String>,
    /// The directory where a volume named in an alias file is looked for
    /// before the alias file's own directory: NCBI's "local copy", in the
    /// current directory. `None`: the alias file's directory only.
    local_dir: Option<String>,
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

impl Resolver {
    fn alias_path(&self, base: &str) -> std::path::PathBuf {
        with_extension(Path::new(base), self.seqtype, "al")
    }

    fn index_path(&self, base: &str) -> std::path::PathBuf {
        with_extension(Path::new(base), self.seqtype, "in")
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
        let bytes = std::fs::read(&path).map_err(|source| SeqDbError::Io {
            path: path.clone(),
            source,
        })?;
        let values = read_values(&bytes);
        let mut title = None;
        let mut dblist: &[u8] = &[];
        for (key, value) in &values {
            match key.as_slice() {
                b"TITLE" => title = Some(value.clone()),
                b"DBLIST" => dblist = value,
                other => {
                    if let Some(key) = UNSUPPORTED_KEYS.iter().find(|k| k.as_bytes() == other) {
                        return Err(unsupported(format!(
                            "alias file key {key} is not supported; LOSAT reads TITLE and DBLIST"
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
    // Not ported: NCBI's last step, which looks the name up again in the
    // BLASTDB search path (seqdbalias.cpp:881-917); a name found only there
    // gets the "Could not find" error. NCBI also follows symbolic links in a
    // volume's directory name (`CDirEntry::NormalizePath`) before it compares
    // volume names; LOSAT keeps the names as written and `resolve` rejects the
    // databases where that could change the result.
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
                push_unique(&mut node.vol_names, this_name.to_string());
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
            let group = combine_path(dir_name(&base), "index.alx");
            if exists(Path::new(&group)) {
                return Err(SeqDbError::Unsupported {
                    path: Path::new(&group).to_path_buf(),
                    message: "the directory has a group alias file; LOSAT does not read group \
                              alias files"
                        .to_string(),
                });
            }
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
            match found {
                Some(volume) => push_unique(&mut node.vol_names, volume),
                None => {
                    return Err(SeqDbError::NotFound(format!(
                        "Could not find volume or alias file ({entry}) referenced in alias file \
                         ({this_name})."
                    )));
                }
            }
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
// LOSAT looks for the database at the given path only; NCBI's search path
// (current directory, `BLASTDB`, `.ncbirc`; seqdbcommon.cpp:378-391) is not
// consulted here. Two cases are rejected because LOSAT does not reproduce
// what NCBI derives from them: a volume reached more than once (NCBI's alias
// walk counts it once per mention, seqdbalias.cpp:1588-1592), and two volumes
// of the same base name in different directories (their order would depend
// on the normalised directory names).
/// Resolves the database `db_base` (a path without extension): an alias file
/// when `db_base.{p|n}al` exists, else the single volume `db_base.{p|n}in`.
pub(super) fn resolve(db_base: &str, seqtype: SeqType) -> Result<AliasTree, SeqDbError> {
    // NCBI reads the name as a list (seqdbalias.cpp:114-116, `x_Tokenize`):
    // spaces separate databases and double quotes keep a name with spaces
    // together. LOSAT searches one database.
    let names = split_quoted(db_base.as_bytes());
    let db_base = match names.as_slice() {
        [] => {
            return Err(SeqDbError::NotFound(
                "No database names were passed to CSeqDB::CSeqDB().".to_string(),
            ));
        }
        [name] => String::from_utf8_lossy(name).into_owned(),
        _ => {
            return Err(SeqDbError::Unsupported {
                path: Path::new(db_base).to_path_buf(),
                message: "the name lists several databases; LOSAT searches one database \
                          (a name with spaces goes between double quotes)"
                    .to_string(),
            });
        }
    };
    resolve_one_name(&db_base, seqtype, LocalCopy::CurrentDir)
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

#[cfg(test)]
pub(super) fn test_split_quoted(dbname: &[u8]) -> Vec<&[u8]> {
    split_quoted(dbname)
}
