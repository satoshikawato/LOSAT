//! NCBI's walk of a directory for the databases in it (`FindBlastDBs`, `CBlastDbFinder`): the
//! names that end in `.pin`, `.pal`, `.nin` or `.nal`, in double quotes, sorted. The one owner of
//! the walk: `losat blastdbcmd -list` and `losat blastdbcheck -dir` read it.

use crate::utils::os_bytes::path_from_bytes;

use std::collections::BTreeSet;

use super::SeqType;

/// The separator NCBI puts after a directory (`CDirEntry::GetPathSeparator`).
const SEPARATOR: u8 = if cfg!(windows) { b'\\' } else { b'/' };

/// A database the walk found (`SSeqDBInitInfo`): its name in double quotes, as NCBI keeps
/// it, and its molecule type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct FoundDb {
    pub(crate) quoted: Vec<u8>,
    pub(crate) seqtype: SeqType,
}

impl FoundDb {
    /// `NStr::Replace` of the double quote by `kEmptyStr` in `m_BlastDbName`: the name without the
    /// double quotes.
    pub(crate) fn file_name(&self) -> Vec<u8> {
        self.quoted.iter().copied().filter(|&b| b != b'"').collect()
    }
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:3019-3076
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1350-1364
// The masks are `*.pin`, `*.pal`, `*.nin`, `*.nal` (with case); a file is not opened, so an
// empty or stray `x.pin` is a database, and so is a link that leads nowhere (`IsDir` is
// false for it). The sub-directories (`*`) are entered with `-recursive`, links to
// directories too. The order of the directory does not matter: the list is sorted.
/// `FindFilesInDir` of `dir` (a `CDir`) with `masks`, the sub-directory mask `"*"`, `finder` and
/// the flags `fFF_File` and `recursive`: the databases of the files of `dir` whose names end with
/// one of `extensions` (`.pin` ...).
fn find_files(dir: &[u8], extensions: &[&str], recursive: bool, found: &mut Vec<FoundDb>) {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3879
    // An empty directory is the current one, its names without a prefix (probes K36, K37).
    let read = if dir.is_empty() {
        std::fs::read_dir(".")
    } else {
        std::fs::read_dir(path_from_bytes(dir))
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
        let is_dir =
            || std::fs::metadata(path_from_bytes(full.as_slice())).is_ok_and(|meta| meta.is_dir());
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
        if !std::fs::metadata(path_from_bytes(path.as_slice())).is_ok_and(|meta| meta.is_file()) {
            path = name;
            path.extend_from_slice(index.as_bytes());
        }
        if path.last() != Some(&b'l') {
            continue;
        }
        let Ok(text) = std::fs::read(path_from_bytes(path.as_slice())) else {
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
/// `FindBlastDBs` with `dir`, `dbtype`, `recursive`, true and `remove_redundant`: blastdbcmd's
/// `-list` and blastdbcheck's `-dir` (without `remove_redundant`).
pub(crate) fn find_blast_dbs(
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

#[cfg(test)]
#[path = "find_tests.rs"]
mod tests;
