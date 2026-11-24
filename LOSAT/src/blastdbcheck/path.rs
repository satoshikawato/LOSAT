//! The paths of blastdbcheck's tests as NCBI's `CDirEntry`, `CFile` and `CDir` take them apart and
//! join them, by the crate's one reading of paths (`PathRules::HOST`).

use crate::blastinput::ncbi_environment::PathRules;
use crate::utils::os_bytes::path_from_bytes;

const RULES: PathRules = PathRules::HOST;

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:298-313
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-377
/// The parts of a `CFile` of a path (`CDirEntry::Reset` drops the trailing separators): its
/// directory with the separator after it (`GetDir`: `./` where the path has none), and the base
/// and extension of its name (`GetBase`, `GetExt`).
pub(super) struct FileParts<'a> {
    pub(super) path: &'a [u8],
    pub(super) dir: Vec<u8>,
    pub(super) base: &'a [u8],
    pub(super) ext: &'a [u8],
}

impl<'a> FileParts<'a> {
    pub(super) fn of(path: &'a [u8]) -> Self {
        let path = RULES.reset(path);
        let (base, ext) = RULES.base_and_ext(path);
        FileParts {
            path,
            dir: RULES.dir(path),
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

/// `CDirEntry::MakePath` of a directory and a name (`PathRules::make_path`).
pub(super) fn make_path(dir: &[u8], name: &[u8]) -> Vec<u8> {
    RULES.make_path(dir, name)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3958-3961
/// The directory `dir` as `CDir::GetEntries` lists it and names its entries: the `CDir`'s path
/// (`reset`; `.` where that is empty) with its trailing separator; an entry's path is this and
/// its name.
pub(super) fn entries_base(dir: &[u8]) -> Vec<u8> {
    let path = RULES.reset(dir);
    RULES.add_trailing_separator(if path.is_empty() { b"." } else { path })
}

#[cfg(test)]
#[path = "path_tests.rs"]
mod tests;
