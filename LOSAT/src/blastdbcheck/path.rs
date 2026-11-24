//! The paths of blastdbcheck's tests as NCBI's `CDirEntry`, `CFile` and `CDir` take them apart and
//! join them (`ALL_SEPARATORS`, `GetDir`, `GetBase`, `GetExt`, `MakePath`).

use crate::utils::os_bytes::path_from_bytes;

/// The separator NCBI puts after a directory (`CDirEntry::GetPathSeparator`).
const SEPARATOR: u8 = if cfg!(windows) { b'\\' } else { b'/' };

/// Whether `byte` ends the directory part of a path (`ALL_SEPARATORS`: `/` on Unix; `:`, `/` and
/// `\` on Windows).
fn is_separator(byte: u8) -> bool {
    byte == b'/' || (cfg!(windows) && (byte == b'\\' || byte == b':'))
}

/// Whether `byte` is one of `DIR_SEPARATORS` (`/` on Unix; `/` and `\` on Windows).
fn is_dir_separator(byte: u8) -> bool {
    byte == b'/' || (cfg!(windows) && byte == b'\\')
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-377
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:413-462
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:475-483
/// The parts of a `CFile` of a path (`CDirEntry::Reset` drops the trailing separators of all but
/// a lone separator): its directory with the separator after it (`GetDir`: `./` where the path has
/// none), and the name after it cut at its last `.` (`GetBase`, `GetExt`: the extension keeps its
/// `.`; a name without one has an empty extension).
pub(super) struct FileParts<'a> {
    pub(super) path: &'a [u8],
    pub(super) dir: Vec<u8>,
    pub(super) base: &'a [u8],
    pub(super) ext: &'a [u8],
}

impl<'a> FileParts<'a> {
    pub(super) fn of(path: &'a [u8]) -> Self {
        let path = reset(path);
        let name_at = path
            .iter()
            .rposition(|&byte| is_separator(byte))
            .map(|at| at + 1);
        let dir = match name_at {
            Some(at) => path[..at].to_vec(),
            None if path.is_empty() => Vec::new(),
            None => vec![b'.', SEPARATOR],
        };
        let name = &path[name_at.unwrap_or(0)..];
        let (base, ext) = match name.iter().rposition(|&byte| byte == b'.') {
            Some(dot) => name.split_at(dot),
            None => (name, &name[name.len()..]),
        };
        FileParts {
            path,
            dir,
            base,
            ext,
        }
    }

    /// `CFile::Exists` (a regular file, links followed) and `GetLength` (its size; -1 where it
    /// cannot be read): the size of a regular file, `None` where there is none.
    pub(super) fn file_length(&self) -> Option<u64> {
        std::fs::metadata(path_from_bytes(self.path))
            .ok()
            .filter(|metadata| metadata.is_file())
            .map(|metadata| metadata.len())
    }
}

/// `CDirEntry::MakePath` of a directory and a name: the name after the directory and a separator
/// (one is added where the directory does not end with one); the name alone after an empty
/// directory. An absolute name is joined as any other.
pub(super) fn make_path(dir: &[u8], name: &[u8]) -> Vec<u8> {
    let mut path = dir.to_vec();
    if path.last().is_some_and(|&byte| !is_separator(byte)) {
        path.push(SEPARATOR);
    }
    path.extend_from_slice(name);
    path
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:475-483
/// `CDirEntry::Reset` of `path`: without its trailing `DIR_SEPARATORS`, but a lone separator kept
/// (empty when the path has nothing else).
fn reset(path: &[u8]) -> &[u8] {
    if path.len() == 1 && is_separator(path[0]) {
        return path;
    }
    match path.iter().rposition(|&byte| !is_dir_separator(byte)) {
        Some(at) => &path[..=at],
        None => &[],
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3958-3961
/// The directory `dir` as `CDir::GetEntries` lists it and names its entries: the `CDir`'s path
/// (`reset`; `.` where that is empty) with a separator after it; an entry's path is this and its
/// name.
pub(super) fn entries_base(dir: &[u8]) -> Vec<u8> {
    let path = reset(dir);
    let mut base = if path.is_empty() {
        b".".to_vec()
    } else {
        path.to_vec()
    };
    if !base.last().is_some_and(|&byte| is_separator(byte)) {
        base.push(SEPARATOR);
    }
    base
}

#[cfg(test)]
#[path = "path_tests.rs"]
mod tests;
