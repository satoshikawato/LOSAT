//! NCBI's `CSeqDBAtlas` as far as the multithreaded scan reads it: the largest file the atlas has
//! measured so far (`m_MaxFileSize`), which sets the budget of the chunks that
//! `CSeqDBImpl::x_FillSeqBuffer` cuts ([`super::mt_scan::mt_chunk_budget`]).
//!
//! The atlas measures a file the first time `GetFileSizeL` asks for it (a cache keeps the answer)
//! and keeps the largest length it measured. NCBI's atlas is one object for the process while a
//! database is open, so the files of every database the run opens count (the repeat database of
//! `-filtering_db` too), and the largest only grows. [`Atlas`] is the one owner of that number in
//! LOSAT: the readers note each file at the point where NCBI's reader sizes it, and the scan asks
//! for the largest when it cuts its chunks.
//!
//! Where NCBI sizes files (`GetFileSizeL`, directly or through `DoesFileExist` and
//! `CSeqDBRawFile::Open`) in a search of a nucleotide database (probes of line RG D2 tail R2, an
//! `LD_PRELOAD` trace of the stats that `GetFileSizeL` makes in NCBI BLAST+ 2.17.0):
//! - the resolution of the name: every alias file and index file it checks, and the index file of
//!   each volume it opens (`SeqDb::open_search_in`, `SeqDb::open_tree`);
//! - the build of the OID list: an `OIDLIST` file, an alias file's `GILIST`, `TILIST` or (version
//!   4) `SEQIDLIST` file and the ISAM files that translate it in each volume (`alias::alias_oids`);
//!   a version 5 `SEQIDLIST` or `TAXIDLIST` is read without the atlas;
//! - the mask columns' files of every volume, with `-db_soft_mask` or `-db_hard_mask`
//!   (`SeqDb::mask_algorithms`);
//! - the set-up of each query batch's search, where the statistics' total length is not positive:
//!   the sequence file of OID 0's volume (`SeqDb::mt_check_past_range_end`);
//! - the scan: the sequence file of each volume it reads, at its first read (the same);
//! - after a query batch's search: the header file of the volume of each subject whose
//!   Seq-align the batch makes, and the ISAM and header files that the data loader's lookups of
//!   the Seq-align ids read (`SeqDb::note_header_reads`, `SeqDb::note_lookup_reads`);
//! - the repeat search of `-filtering_db`, before each batch's search: the files of the repeat
//!   database, measured by its own reader (`Atlas::absorb`).
//!
//! The largest file sets a budget that matters only where the scan's iteration range ends inside
//! a volume after the first (`SeqDb::mt_check_past_range_end`), which an alias file's
//! `FIRST_OID` and `LAST_OID` give. A search's own list (`-gilist`, `-seqidlist`, `-taxids` and
//! their negative forms) sizes the ISAM files that translate it, and makes the traceback's fetches
//! read headers, but NCBI then sets the iteration range back to every OID (`SeqDb::restrict`), so
//! those files are not noted.
//!
//! The lengths are read when the largest is asked for: a file noted twice, or under two names,
//! changes nothing, and a missing file (or a directory, whose length NCBI takes as -1) measures
//! nothing.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:160-195
/// The files NCBI's atlas has measured: the one owner of the largest file's length
/// (`m_MaxFileSize`).
#[derive(Debug, Default)]
pub(crate) struct Atlas {
    state: Mutex<State>,
}

#[derive(Debug, Default)]
struct State {
    /// Every file noted (NCBI's cache of the lengths, by name).
    noted: HashSet<PathBuf>,
    /// The files noted since the lengths were last read.
    pending: Vec<PathBuf>,
    /// The largest length of the files read so far.
    largest: u64,
}

impl Atlas {
    /// `GetFileSizeL` of `path`: the atlas measures the file.
    pub(crate) fn size(&self, path: impl Into<PathBuf>) {
        let path = path.into();
        let mut state = self.state.lock().expect("the atlas's files");
        if !state.noted.contains(&path) {
            state.noted.insert(path.clone());
            state.pending.push(path);
        }
    }

    /// The files of `other` (another database's reader in the same run): NCBI's atlas is the
    /// process's.
    pub(crate) fn absorb(&self, other: &Atlas) {
        let largest = other.largest();
        let noted: Vec<PathBuf> = {
            let state = other.state.lock().expect("the atlas's files");
            state.noted.iter().cloned().collect()
        };
        for path in noted {
            self.size(path);
        }
        let mut state = self.state.lock().expect("the atlas's files");
        state.largest = state.largest.max(largest);
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:496-501,597
    /// `m_MaxFileSize`: the largest length of the files measured so far (0 where none exists).
    pub(crate) fn largest(&self) -> u64 {
        let mut state = self.state.lock().expect("the atlas's files");
        let pending = std::mem::take(&mut state.pending);
        for path in pending {
            if let Ok(meta) = std::fs::metadata(&path) {
                if meta.is_file() {
                    state.largest = state.largest.max(meta.len());
                }
            }
        }
        state.largest
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:58-121
    /// `CSeqDBIsam::x_InitSearch` of the ISAM of `index` and `data` (a string ISAM where
    /// `string`): it measures the index file, and the data file where the index holds the 40
    /// bytes of its header, NCBI's version and the ISAM's type, and is not a memory-only index
    /// (page size 1).
    pub(crate) fn isam(&self, index: &Path, data: &Path, string: bool) {
        self.size(index);
        if isam_measures_data(index, string) {
            self.size(data);
        }
    }
}

/// Whether `x_InitSearch` of the index file `index` goes on to measure the data file
/// ([`Atlas::isam`]).
fn isam_measures_data(index: &Path, string: bool) -> bool {
    // `ISAM_VERSION`; the string type (`eString`), the numeric types (`eNumeric`,
    // `eNumericLongId`); `MEMORY_ONLY_PAGE_SIZE`.
    const VERSION: i32 = 1;
    const STRING: i32 = 2;
    const NUMERIC: i32 = 0;
    const NUMERIC_LONG_ID: i32 = 5;
    const MEMORY_ONLY_PAGE_SIZE: i32 = 1;
    let Ok(bytes) = std::fs::read(index) else {
        return false;
    };
    if bytes.len() < 40 {
        return false;
    }
    let field = |i: usize| {
        i32::from_be_bytes([
            bytes[4 * i],
            bytes[4 * i + 1],
            bytes[4 * i + 2],
            bytes[4 * i + 3],
        ])
    };
    let type_fits = if string {
        field(1) == STRING
    } else {
        matches!(field(1), NUMERIC | NUMERIC_LONG_ID)
    };
    field(0) == VERSION && type_fits && field(5) != MEMORY_ONLY_PAGE_SIZE
}

#[cfg(test)]
mod tests {
    use super::Atlas;

    #[test]
    fn the_largest_is_the_longest_file_noted_and_missing_files_count_nothing() {
        let dir = crate::utils::test_dir::TempDir::unique("atlas");
        let a = dir.join("a");
        let b = dir.join("b");
        std::fs::write(&a, vec![0u8; 120]).unwrap();
        std::fs::write(&b, vec![0u8; 70]).unwrap();
        let atlas = Atlas::default();
        assert_eq!(atlas.largest(), 0);
        atlas.size(dir.join("missing"));
        atlas.size(dir.to_path_buf());
        assert_eq!(atlas.largest(), 0);
        atlas.size(&b);
        assert_eq!(atlas.largest(), 70);
        atlas.size(&a);
        atlas.size(&b);
        assert_eq!(atlas.largest(), 120);
        let other = Atlas::default();
        let c = dir.join("c");
        std::fs::write(&c, vec![0u8; 200]).unwrap();
        other.size(&c);
        atlas.absorb(&other);
        assert_eq!(atlas.largest(), 200);
    }

    #[test]
    fn an_isam_measures_its_data_file_only_after_a_valid_header() {
        let dir = crate::utils::test_dir::TempDir::unique("atlas");
        let index = dir.join("v.nni");
        let data = dir.join("v.nnd");
        std::fs::write(&data, vec![0u8; 500]).unwrap();
        let header = |version: i32, kind: i32, page: i32| {
            let mut bytes = Vec::new();
            for word in [version, kind, 500, 1, 1, page, 0, 0, 0, 0] {
                bytes.extend_from_slice(&word.to_be_bytes());
            }
            bytes
        };
        for (bytes, string, measured) in [
            (header(1, 0, 256), false, 500),
            (header(1, 5, 256), false, 500),
            (header(1, 2, 64), true, 500),
            (header(1, 2, 64), false, 40),
            (header(2, 0, 256), false, 40),
            (header(1, 0, 1), false, 40),
            (header(1, 0, 256)[..36].to_vec(), false, 36),
        ] {
            std::fs::write(&index, &bytes).unwrap();
            let atlas = Atlas::default();
            atlas.isam(&index, &data, string);
            assert_eq!(atlas.largest(), measured, "{bytes:?} {string}");
        }
    }
}
