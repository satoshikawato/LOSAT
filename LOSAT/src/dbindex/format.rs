//! The layout of NCBI's MegaBLAST index files, the one owner for the writer (`losat
//! makembindex`: [`super::factory`], `crate::makembindex::work`) and the reader (BLASTN's index
//! loader, `algorithm::blastn::blast_engine::mb_index`): the names of the index volumes, the
//! super header of a database volume's index (`<volume>.shd`) and the header at the start of
//! each index volume (`<name>.NN.idx`).
//!
//! The files are NCBI's: little-endian 32-bit words (NCBI writes and reads the host's bytes;
//! LOSAT's targets are little-endian).

use std::io::Write;
use std::path::Path;

use super::{SeqNum, Word, OFFSET_COMBINED, UNCOMPRESSED, VERSION, WIDTH_32};

/// `<base>.<NN>.idx`: the name of index volume `volume` (two digits at least).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:105-112
pub(crate) fn volume_name(base: &[u8], volume: u32) -> Vec<u8> {
    [base, format!(".{volume:02}.idx").as_bytes()].concat()
}

/// The `.shd` super header of the index of one database volume (`CIndexSuperHeader<1>`): the
/// volume's number of sequences and the number of its index volumes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SuperHeader {
    pub(crate) num_seq: u32,
    pub(crate) num_vol: u32,
}

/// The size of a version 1 super header (`EXPECTED_SIZE`): four words.
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:197-199
const SUPER_HEADER_SIZE: u64 = 16;

impl SuperHeader {
    /// `Save` of the `CIndexSuperHeader` of index format version 1 (`INDEX_FORMAT_VERSION_1`) to
    /// `fname`: the endianness (0, little), the format version (1), the number of sequences and the
    /// number of volumes, each a 32-bit word. `std::ofstream` reports no failure here
    /// (`CheckStream` tests `bad()` only, and a stream that did not open is only `fail()`), so none
    /// is returned.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:196-203
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:241-252
    pub(crate) fn save(&self, fname: &[u8]) {
        let Ok(file) = std::fs::File::create(crate::utils::os_bytes::path_of_bytes(fname)) else {
            return;
        };
        let mut out = std::io::BufWriter::new(file);
        for word in [0u32, 1, self.num_seq, self.num_vol] {
            let _ = out.write_all(&word.to_le_bytes());
        }
        let _ = out.flush();
    }

    /// `GetIndexSuperHeader` of `fname`: the super header in the file `fname`, or `None` where NCBI
    /// throws (its loader only posts that at the level `Info`, which is not shown, and takes the
    /// database volume as one without an index). NCBI reads a super header only from a readable
    /// regular file (`CFile::Exists`, `CheckAccess`) whose first word is 0 (little-endian, the
    /// host's), whose second is the version 1, and whose size is 16 bytes; a shorter file ends
    /// at one of these checks whatever the words that its reads leave unset.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:128-181
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:208-228
    pub(crate) fn load(fname: &Path) -> Option<SuperHeader> {
        let meta = std::fs::metadata(fname).ok()?;
        if !meta.is_file() || meta.len() != SUPER_HEADER_SIZE {
            return None;
        }
        let bytes = std::fs::read(fname).ok()?;
        let words = words(bytes.get(..16)?);
        (words[0] == 0 && words[1] == 1).then_some(SuperHeader {
            num_seq: words[2],
            num_vol: words[3],
        })
    }
}

/// The little-endian words of `bytes` (a multiple of four bytes).
fn words(bytes: &[u8]) -> Vec<Word> {
    bytes
        .chunks_exact(4)
        .map(|word| Word::from_le_bytes([word[0], word[1], word[2], word[3]]))
        .collect()
}

/// The size of an index volume's header (`HEADER_SIZE`): the version and the width (8 bytes
/// each), then seven words.
// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex_sp.hpp:343
pub(crate) const HEADER_SIZE: usize = 16 + 7 * 4;

/// The header of an index volume: the version byte (5 legacy, 6 new), the Nmer length, the
/// two words after it (the stride and the word size hint of the new format; `OFFSET_COMBINED`
/// and `UNCOMPRESSED` in a legacy one), and the first and last sequences and chunks of the
/// volume.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Header {
    pub(crate) version: u8,
    pub(crate) hkey_width: Word,
    pub(crate) stride_word: Word,
    pub(crate) ws_hint_word: Word,
    pub(crate) start: SeqNum,
    pub(crate) start_chunk: SeqNum,
    pub(crate) stop: SeqNum,
    pub(crate) stop_chunk: SeqNum,
}

impl Header {
    /// The header that `CDbIndex_Factory::SaveHeader` writes for an index of `hkey_width`-mers,
    /// legacy (`OFFSET_COMBINED` and `UNCOMPRESSED`) or not (`stride` and `ws_hint`), of the
    /// sequences and chunks `[start, start_chunk, stop, stop_chunk]`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1604-1634
    pub(crate) fn of(
        legacy: bool,
        hkey_width: Word,
        stride: Word,
        ws_hint: Word,
        [start, start_chunk, stop, stop_chunk]: [SeqNum; 4],
    ) -> Header {
        let (version, stride_word, ws_hint_word) = if legacy {
            (VERSION, OFFSET_COMBINED, UNCOMPRESSED)
        } else {
            (VERSION + 1, stride, ws_hint)
        };
        Header {
            version,
            hkey_width,
            stride_word,
            ws_hint_word,
            start,
            start_chunk,
            stop,
            stop_chunk,
        }
    }

    /// `SaveHeader`'s bytes: the version and seven zero bytes, the width (`WIDTH_32`, 8 bytes),
    /// then the seven words.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1604-1634
    pub(crate) fn write(&self, os: &mut impl Write) -> std::io::Result<()> {
        os.write_all(&[self.version, 0, 0, 0, 0, 0, 0, 0])?;
        os.write_all(&WIDTH_32.to_le_bytes())?;
        for word in [
            self.hkey_width,
            self.stride_word,
            self.ws_hint_word,
            self.start,
            self.start_chunk,
            self.stop,
            self.stop_chunk,
        ] {
            os.write_all(&word.to_le_bytes())?;
        }
        Ok(())
    }

    /// `ReadIndexHeader<false>` of the mapped bytes `map` (the version byte is not read):
    /// the words after the version and the width, or `None` where `CheckIndexEndianness`
    /// throws (an Nmer length over 15). A mapping shorter than the header reads zeros past the
    /// end of the file (the rest of the mapped page), which `map` gets here.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:252-277
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:303-349
    pub(crate) fn read(map: &[u8]) -> Option<Header> {
        let mut bytes = [0u8; HEADER_SIZE];
        let length = map.len().min(HEADER_SIZE);
        bytes[..length].copy_from_slice(&map[..length]);
        let words = words(&bytes[16..]);
        const MAX_KEY_SIZE: Word = 15;
        if words[0] > MAX_KEY_SIZE {
            return None;
        }
        Some(Header {
            version: bytes[0],
            hkey_width: words[0],
            stride_word: words[1],
            ws_hint_word: words[2],
            start: words[3],
            start_chunk: words[4],
            stop: words[5],
            stop_chunk: words[6],
        })
    }

    /// `s_GetMinimumSupportedWordSizeByIndex`: `hkey_width + stride - 1` with the header read
    /// as a new-format one whatever its version, so that a legacy volume gives its Nmer length
    /// (its `OFFSET_COMBINED`, 1, stands where the stride would). NCBI sums `unsigned long`s and
    /// keeps the `int`: the low 32 bits.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:909-916
    pub(crate) fn min_word_size(&self) -> i32 {
        (u64::from(self.hkey_width) + u64::from(self.stride_word)).wrapping_sub(1) as i32
    }
}

/// What `GetIdxVolNumOIDs` of `fname` gives for an index volume.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum VolumeOids {
    /// The number of sequences of the volume, `stop - start` of its header (32-bit, wrapping).
    Count(u32),
    /// The stream ended before the header's words: NCBI's `CheckStream` throws with this place
    /// (`at start oid`, `at end oid`).
    ReadFailure(&'static str),
    /// A read that NCBI does not check leaves the count undefined: a volume that does not open,
    /// or that ends inside the header's stop word.
    Undefined,
}

/// `GetIdxVolNumOIDs` of `fname`: the number of sequences of the index volume `fname` from its
/// header, read from a stream that `CheckStream` checks for its end after the first seven words
/// and after the start word and the next.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:232-245
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:87-93
pub(crate) fn volume_oids(fname: &Path) -> VolumeOids {
    let Ok(mut file) = std::fs::File::open(fname) else {
        return VolumeOids::Undefined;
    };
    let mut bytes = Vec::with_capacity(40);
    if std::io::Read::read_to_end(&mut std::io::Read::take(&mut file, 40), &mut bytes).is_err() {
        return VolumeOids::Undefined;
    }
    match bytes.len() {
        0..=27 => VolumeOids::ReadFailure("at start oid"),
        28..=35 => VolumeOids::ReadFailure("at end oid"),
        36..=39 => VolumeOids::Undefined,
        _ => {
            let words = words(&bytes);
            VolumeOids::Count(words[9].wrapping_sub(words[7]))
        }
    }
}

#[cfg(test)]
#[path = "format_tests.rs"]
mod tests;
