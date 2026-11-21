//! One MegaBLAST index volume (`<name>.NN.idx`) as BLASTN's indexed search reads it
//! (`CDbIndex::Load`, `CDbIndex_Impl`): the file mapped read-only, its header, and the places of
//! its sections, which [`super::search`] reads: the Nmer hash table and the offset lists
//! (`COffsetData`), then the subject map (`CSubjectMap`: the sequence lengths, the local ids, the
//! first chunk of each sequence, the start of each chunk in the packed sequences) and the packed
//! sequences themselves (two bits a base, four bases a byte from the high bits down).
//!
//! The layout is the writer's ([`super::factory`], [`super::format`]). NCBI reads the sections
//! by their sizes without checking them against the file; LOSAT checks that every section lies
//! in the file and that the subject map's words point inside their sections, and refuses a
//! volume that does not (NCBI's search of it reads past its sections).

use std::path::Path;

use super::format::{Header, HEADER_SIZE};
use super::{min_offset, Word, MAX_STRIDE, STRIDE, VERSION};
use crate::utils::mapped_file::MappedFile;

/// The word size hint of a legacy index (`CDbIndex::WS_HINT`, used in place of the header's).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:352-367
const LEGACY_WS_HINT: u64 = 28;

/// The page size of the mapping of NCBI's reference platform (Linux x86-64).
const PAGE_SIZE: usize = 4096;

/// A section: its first word (or byte, for the packed sequences) and its length.
#[derive(Clone, Copy, Debug)]
struct Section {
    at: usize,
    len: usize,
}

/// A mapped index volume and the places of its sections.
pub(crate) struct IndexVolume {
    map: MappedFile,
    /// `hkey_width`: the Nmer length.
    pub(crate) hkey_width: u64,
    /// The stride and the word size hint (5 and 28 for a legacy volume, the header's for the
    /// other format).
    pub(crate) stride: u64,
    pub(crate) ws_hint: u64,
    /// The first and one past the last sequence of the volume (`StartSeq`, `StopSeq`).
    pub(crate) start: Word,
    pub(crate) stop: Word,
    /// `offset_bits_` of the subject map.
    pub(crate) offset_bits: u32,
    hash: Section,
    data: Section,
    lid_map: Section,
    subjects: Section,
    chunks: Section,
    seq_store: Section,
}

/// The words of a section of `map` (checked when the volume was opened).
fn words(map: &MappedFile, section: Section) -> &[Word] {
    map.u32s(section.at * 4, section.len)
        .expect("the section was checked when the volume was opened")
}

impl IndexVolume {
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:413-433
    // NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex_sp.hpp:467-501
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:449-565
    /// `CDbIndex::Load(fname)` of the index volume `path`: its version (5 legacy, 6 the other
    /// format), its header, then the sections in the writer's order. `Err` says why LOSAT does
    /// not read it (a file that does not open, a version NCBI rejects, or a section that does
    /// not lie in the file).
    pub(crate) fn open(path: &Path) -> Result<IndexVolume, String> {
        let map = MappedFile::open(path).map_err(|error| error.to_string())?;
        let bytes = map.bytes();
        if bytes.len() < HEADER_SIZE {
            return Err("shorter than its header".to_string());
        }
        let legacy = match bytes[0] {
            VERSION => true,
            version if version == VERSION + 1 => false,
            _ => return Err("wrong index version".to_string()),
        };
        let header = Header::read(bytes).ok_or("index header validation failed")?;
        let (stride, ws_hint) = if legacy {
            (STRIDE, LEGACY_WS_HINT)
        } else {
            (
                u64::from(header.stride_word),
                u64::from(header.ws_hint_word),
            )
        };
        // makembindex writes neither; NCBI's search of an Nmer length of 0 takes `hkey_width - 1`
        // below 0, and `GetMinOffset` of a stride from 2^15 on shifts an `int` by 32 bits or more.
        if header.hkey_width == 0 {
            return Err("an Nmer length of 0".to_string());
        }
        if stride == 0 || stride >= MAX_STRIDE {
            return Err(format!("a stride of {stride}"));
        }
        let total_words = bytes.len() / 4;
        let mut at = HEADER_SIZE / 4;
        let word = |at: usize| -> Result<Word, String> {
            map.u32s(at * 4, 1)
                .map(|w| w[0])
                .ok_or_else(|| "a section past the end of the file".to_string())
        };
        let section = |at: usize, len: usize| -> Result<Section, String> {
            match at.checked_add(len) {
                Some(end) if end <= total_words => Ok(Section { at, len }),
                _ => Err("a section past the end of the file".to_string()),
            }
        };
        // `COffsetData_Base`: the number of words of the offset lists, the hash table
        // (`4^hkey_width + 1` words), then the lists.
        let data_len = word(at)? as usize;
        at += 1;
        let hash = section(at, (1usize << (2 * header.hkey_width)) + 1)?;
        at += hash.len;
        let data = section(at, data_len)?;
        at += data.len;
        // `CSubjectMap(map, header)`: the bytes of the lengths and the offset bits, the lengths,
        // the bytes of the local id map and the map (four words a local id), then `Load`.
        let lengths_len = word(at)? as usize / 4;
        let offset_bits = word(at + 1)?;
        at += 2;
        let lengths = section(at, lengths_len)?;
        at += lengths.len;
        let lid_map_len = word(at)? as usize / 4;
        at += 1;
        let lid_map = section(at, lid_map_len)?;
        at += lid_map.len;
        // `Load`: the bytes of the first chunks and the chunk starts; the first chunk of each
        // sequence (`stop - start + 1` words); the chunk starts and one more word (the first
        // copy of the packed sequences' size); then the size again and the packed sequences.
        let subjects_and_chunks = word(at)?;
        at += 1;
        let n_subjects = (header.stop.wrapping_sub(header.start) as usize) + 1;
        let chunk_bytes = u64::from(subjects_and_chunks)
            .checked_sub(4 * n_subjects as u64)
            .ok_or("a subject map shorter than its sequences")?;
        let subjects = section(at, n_subjects)?;
        at += subjects.len;
        let chunks = section(at, 1 + (chunk_bytes / 4) as usize)?;
        at += chunks.len;
        let seq_store_len = word(at)? as usize;
        at += 1;
        let seq_store = Section {
            at: at * 4,
            len: seq_store_len,
        };
        if seq_store
            .at
            .checked_add(seq_store.len)
            .is_none_or(|end| end > bytes.len())
        {
            return Err("packed sequences past the end of the file".to_string());
        }
        // `offset_mask_ = (1 << offset_bits_) - 1` of an `int`.
        if !(1..=30).contains(&offset_bits) {
            return Err("offset bits out of range".to_string());
        }
        let volume = IndexVolume {
            map,
            hkey_width: u64::from(header.hkey_width),
            stride,
            ws_hint,
            start: header.start,
            stop: header.stop,
            offset_bits,
            hash,
            data,
            lid_map,
            subjects,
            chunks,
            seq_store,
        };
        volume.check_subject_map()?;
        Ok(volume)
    }

    /// The subject map's words point inside their sections: each local id's chunks (a
    /// non-empty range of the chunk starts, the first at the local id's first byte, which lies
    /// in the packed sequences), and each sequence's first chunk (0 for a sequence without
    /// chunks). (A local id's last word is its first byte plus its number of bases, which the
    /// search's right extension takes as a byte bound, as NCBI's does.)
    fn check_subject_map(&self) -> Result<(), String> {
        let chunks = self.chunks();
        let seq_store_len = self.seq_store().len() as u64;
        for lid in self.lid_map().chunks_exact(4) {
            let [start, end, seq_start] = [lid[0], lid[1], lid[2]];
            if start >= end
                || end as usize > chunks.len()
                || chunks[start as usize] != seq_start
                || u64::from(seq_start) > seq_store_len
            {
                return Err("a local id outside its sections".to_string());
            }
        }
        if self
            .subjects()
            .iter()
            .any(|&first| first as usize > chunks.len())
        {
            return Err("a sequence's first chunk outside the chunks".to_string());
        }
        Ok(())
    }

    /// The Nmer hash table: for each Nmer the 1-based start of its list in [`Self::offsets`]
    /// (0 for none), and one more word.
    pub(crate) fn hash_table(&self) -> &[Word] {
        words(&self.map, self.hash)
    }

    /// The offset lists (`data_start_`, whose first word is the writer's leading 0).
    pub(crate) fn offsets(&self) -> &[Word] {
        words(&self.map, self.data)
    }

    /// The local id map: four words a local id (its first chunk, one past its last chunk, its
    /// first byte and one past its last byte in the packed sequences).
    pub(crate) fn lid_map(&self) -> &[Word] {
        words(&self.map, self.lid_map)
    }

    /// The 1-based first chunk of each sequence of the volume (and one more word).
    pub(crate) fn subjects(&self) -> &[Word] {
        words(&self.map, self.subjects)
    }

    /// The first byte of each chunk in the packed sequences, and the sentinel word after them
    /// (`NumChunks` counts it).
    pub(crate) fn chunks(&self) -> &[Word] {
        words(&self.map, self.chunks)
    }

    /// The packed sequences.
    pub(crate) fn seq_store(&self) -> &[u8] {
        &self.map.bytes()[self.seq_store.at..self.seq_store.at + self.seq_store.len]
    }

    /// The zero bytes that NCBI's mapping of the file shows after the packed sequences, which
    /// end the file: the rest of the file's last page of 4096 bytes (NCBI's search can read
    /// them where an exact match runs past the volume's last sequence; past them it crashes).
    pub(crate) fn seq_store_padding(&self) -> usize {
        let end = self.seq_store.at + self.seq_store.len;
        if end == self.map.len() {
            (PAGE_SIZE - end % PAGE_SIZE) % PAGE_SIZE
        } else {
            0
        }
    }

    /// `GetMinOffset(stride)`: the first offset value that is not a boundary code.
    pub(crate) fn min_offset(&self) -> u64 {
        min_offset(self.stride)
    }
}

#[cfg(test)]
#[path = "volume_tests.rs"]
mod tests;
