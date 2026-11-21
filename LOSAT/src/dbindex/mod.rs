//! NCBI's database index for megablast (`algo/blast/dbindex`, `CDbIndex`): the index volumes
//! (`.NN.idx`) and super header (`.shd`) that `losat makembindex` writes (`factory`, NCBI's
//! `dbindex_factory.cpp`) from the sequences of a BLAST database or a FASTA file (`stream`,
//! NCBI's `sequence_istream_bdb.cpp` and `sequence_istream_fasta.cpp`).
//!
//! The files are NCBI's: little-endian 32-bit words (`WriteWord` writes the host's bytes;
//! LOSAT's targets are little-endian), byte-identical to those of NCBI BLAST+ 2.17.0's
//! makembindex.

pub(crate) mod factory;
pub(crate) mod format;
pub(crate) mod search;
pub(crate) mod stream;
pub(crate) mod volume;

use crate::core::subject_split::{DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN};

/// The name in LOSAT's refusals.
pub(crate) const PROGRAM: &str = "MAKEMBINDEX";

/// LOSAT's refusal of `what` in makembindex.
pub(crate) fn unsupported(what: &str) -> String {
    crate::seqmasks_io::run::unsupported(what, PROGRAM)
}

/// `CDbIndex::TWord`: the unit of the index.
pub(crate) type Word = u32;

/// `CDbIndex::TSeqNum` (`CSequenceIStream::TStreamPos`): the number of a sequence of the input.
pub(crate) type SeqNum = u32;

/// `CDbIndex::CR`: the bases in one byte of the packed sequences.
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:441
pub(crate) const CR: u64 = 4;

/// `CDbIndex::STRIDE`: the default stride.
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:447
pub(crate) const STRIDE: u64 = 5;

/// `CDbIndex::VERSION`: the first byte of a legacy index volume (the new format's is one more).
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:463
pub(crate) const VERSION: u8 = 5;

/// `WIDTH_32`, `OFFSET_COMBINED` and `UNCOMPRESSED` of the legacy header.
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:46-54
pub(crate) const WIDTH_32: u64 = 0;
pub(crate) const OFFSET_COMBINED: Word = 1;
pub(crate) const UNCOMPRESSED: Word = 0;

/// `CDbIndex::SOptions`: the parameters of index creation (`report_level` is read by nothing).
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:467-481
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Options {
    /// `idmap`: write `<volume>.map`, the first word of each sequence's title.
    pub(crate) idmap: bool,
    /// `legacy`: the legacy header (version 5) in place of the new one (version 6).
    pub(crate) legacy: bool,
    pub(crate) stride: u64,
    pub(crate) ws_hint: u64,
    /// `hkey_width`: the length of the indexed words (the Nmers).
    pub(crate) hkey_width: u64,
    pub(crate) chunk_size: u64,
    pub(crate) chunk_overlap: u64,
    /// `max_index_size`: the volume size in megabytes.
    pub(crate) max_index_size: u64,
    /// `stat_file_name`; `None` for the empty name.
    pub(crate) stat_file_name: Option<String>,
}

impl Options {
    /// `CDbIndex::DefaultSOptions()`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:383-398
    pub(crate) fn default_options() -> Self {
        Options {
            idmap: false,
            legacy: true,
            stride: STRIDE,
            ws_hint: 28,
            hkey_width: 12,
            // BLAST's chunks of a long sequence (`MAX_DBSEQ_LEN`, `DBSEQ_CHUNK_OVERLAP`).
            chunk_size: MAX_DBSEQ_LEN as u64,
            chunk_overlap: DBSEQ_CHUNK_OVERLAP as u64,
            max_index_size: 1536,
            stat_file_name: None,
        }
    }
}

/// `GetCodeBits(stride)`: the bits of `stride`, at least 1.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:438-443
pub(crate) fn code_bits(stride: u64) -> u64 {
    let mut stride = stride;
    let mut result = 0;
    loop {
        stride >>= 1;
        result += 1;
        if stride == 0 {
            return result;
        }
    }
}

/// `GetMinOffset(stride)`: 4 to the power `code_bits(stride)`, which NCBI computes by shifting
/// an `int`; defined for a stride below 32768 (`code_bits` 15 or less), which `makembindex`
/// checks first (`MAX_STRIDE`).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:446-447
pub(crate) fn min_offset(stride: u64) -> u64 {
    1u64 << (2 * code_bits(stride))
}

/// The strides below this one are those of a defined `GetMinOffset`: for larger ones NCBI
/// shifts an `int` by 32 bits or more, which C++ leaves undefined.
pub(crate) const MAX_STRIDE: u64 = 1 << 15;

/// The index volume `<dir>/<name>.00.idx` of the FASTA text `fasta`, made with `options` by
/// LOSAT's writer (the tests of `volume` and `search`).
#[cfg(test)]
pub(crate) fn test_volume(name: &str, fasta: &str, options: &Options) -> std::path::PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-dbindex-{name}-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let base = dir.join(name).to_string_lossy().into_owned();
    let mut input = stream::FastaStream::new(
        Box::new(std::io::Cursor::new(fasta.as_bytes().to_vec())),
        false,
    );
    let path = format::volume_name(&base, 0);
    let mut stop = SeqNum::MAX;
    let made = factory::make_index(
        &mut input,
        &path,
        0,
        &mut stop,
        options,
        &mut |_| {},
        &mut |_| {},
    );
    assert!(made.is_ok(), "the test index was not made");
    std::path::PathBuf::from(path)
}

/// Letters `A`, `C`, `G`, `T` from a fixed generator (the tests' sequences).
#[cfg(test)]
pub(crate) fn test_letters(len: usize, seed: u64) -> String {
    let mut state = seed | 1;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state >> 20) as usize % 4] as char
        })
        .collect()
}
