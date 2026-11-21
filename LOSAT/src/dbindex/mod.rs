//! NCBI's database index for megablast (`algo/blast/dbindex`, `CDbIndex`): the index volumes
//! (`.NN.idx`) and super header (`.shd`) that `losat makembindex` writes (`factory`, NCBI's
//! `dbindex_factory.cpp`) from the sequences of a BLAST database or a FASTA file (`stream`,
//! NCBI's `sequence_istream_bdb.cpp` and `sequence_istream_fasta.cpp`).
//!
//! The files are NCBI's: little-endian 32-bit words (`WriteWord` writes the host's bytes;
//! LOSAT's targets are little-endian), byte-identical to those of NCBI BLAST+ 2.17.0's
//! makembindex.

pub(crate) mod factory;
pub(crate) mod stream;

use std::io::Write;

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

/// `MAX_DBSEQ_LEN` and `DBSEQ_CHUNK_OVERLAP` of BLAST's core: the chunks of a long sequence.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:192
pub(crate) const MAX_DBSEQ_LEN: u64 = 5_000_000;
pub(crate) const DBSEQ_CHUNK_OVERLAP: u64 = 100;

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
            chunk_size: MAX_DBSEQ_LEN,
            chunk_overlap: DBSEQ_CHUNK_OVERLAP,
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

/// `GetMinOffset(stride)`: `1 << (2 * GetCodeBits(stride))`, a shift of NCBI's `int`; defined
/// for a stride below 32768 (`code_bits` 15 or less), which `makembindex` checks first
/// (`MAX_STRIDE`).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:446-447
pub(crate) fn min_offset(stride: u64) -> u64 {
    1u64 << (2 * code_bits(stride))
}

/// The strides below this one are those of a defined `GetMinOffset` (NCBI shifts its `int` 1
/// by 32 or more bits for larger ones).
pub(crate) const MAX_STRIDE: u64 = 1 << 15;

/// `CIndexSuperHeader<INDEX_FORMAT_VERSION_1>::Save(fname)`: the endianness (0, little), the
/// format version (1), the number of sequences and the number of volumes, each a 32-bit word.
/// `std::ofstream` reports no failure here (`CheckStream` tests `bad()` only, and a stream that
/// did not open is only `fail()`), so none is returned.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:196-203
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:241-252
pub(crate) fn save_super_header(fname: &str, num_seq: u32, num_vol: u32) {
    let Ok(file) = std::fs::File::create(fname) else {
        return;
    };
    let mut out = std::io::BufWriter::new(file);
    for word in [0u32, 1, num_seq, num_vol] {
        let _ = out.write_all(&word.to_le_bytes());
    }
    let _ = out.flush();
}

/// `<base>.<NN>.idx`: the name of index volume `volume` (two digits at least).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:105-112
pub(crate) fn volume_name(base: &str, volume: u32) -> String {
    format!("{base}.{volume:02}.idx")
}
