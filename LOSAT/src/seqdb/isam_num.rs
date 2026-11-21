//! The numeric ISAMs of a volume, GIs (`.pni`/`.pnd`, `.nni`/`.nnd`) and PIGs (`.ppi`/`.ppd`):
//! NCBI's `CSeqDBIsam` for the lookup of a number (`IdToOid`, `PigToOid`). The string ISAM is in
//! [`super::isam`].
//!
//! The index file holds a header and every `page_size`-th row of the data file; the data file
//! holds the rows, sorted by key: a big-endian key (4 bytes, or 8 for `eNumericLongId`) and a
//! big-endian OID of 4 bytes.

use std::path::Path;

use super::file::DbFile;
use super::volume::with_extension;
use super::{SeqDbError, SeqType};

/// `ISAM_VERSION`.
const ISAM_VERSION: i32 = 1;
/// `eNumeric`.
const NUMERIC: i32 = 0;
/// `eNumericLongId`.
const NUMERIC_LONG_ID: i32 = 5;
/// `MEMORY_ONLY_PAGE_SIZE`: an index without a data file, which makeblastdb does not write.
const MEMORY_ONLY_PAGE_SIZE: i32 = 1;
/// `10 * sizeof(Int4)`: the header `x_InitSearch` needs.
const HEADER_BYTES: u64 = 40;
/// `m_KeySampleOffset`.
const KEY_SAMPLE_OFFSET: u64 = 9 * 4;

/// An open numeric ISAM.
#[derive(Debug)]
pub struct NumericIsam {
    index: DbFile,
    data: DbFile,
    /// `m_LongId`: keys of 8 bytes.
    long_id: bool,
    num_terms: i32,
    num_samples: i32,
    page_size: i32,
    /// `m_FirstKey`, `m_LastKey`.
    first_key: i64,
    last_key: i64,
}

fn be_i32(bytes: &[u8]) -> i32 {
    i32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

impl NumericIsam {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:161-175
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:64-118
    // The constructor runs `x_FindIndexBounds` whatever `x_InitSearch` returned. A data file of
    // another length still has its first and last rows: the lookups then find nothing (probe
    // H104). An index with a bad header has no samples and NCBI crashes on it (exit status
    // 139, H98, H101): LOSAT refuses it.
    /// The GI ISAM of the volume `base`: `None` when its files do not open (NCBI's silent
    /// catch) or the data file has another length than the index says (every lookup misses).
    pub fn open_volume(base: &Path, seqtype: SeqType) -> Result<Option<Self>, SeqDbError> {
        let (Ok(index), Ok(data)) = (
            DbFile::open(&with_extension(base, seqtype, "ni")),
            DbFile::open(&with_extension(base, seqtype, "nd")),
        ) else {
            return Ok(None);
        };
        Self::from_files(index, data, "GI")
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:136-145
    /// The PIG ISAM of the volume `base` (`.ppi`/`.ppd`, [`NumericIsam::pig_file_names`]): `None`
    /// when one of its files is no regular file (`CSeqDBIsam::IndexExists`,
    /// `super::isam::open_if_exists`) or the data file has another length than the index says.
    /// A file that is one and does not open is NCBI's exception of the constructor (no catch
    /// here, unlike the GIs).
    pub fn open_pig_volume(base: &Path, seqtype: SeqType) -> Result<Option<Self>, SeqDbError> {
        let (index_path, data_path) = Self::pig_file_names(base, seqtype);
        match super::isam::open_if_exists(&index_path, &data_path)? {
            Some((index, data)) => Self::from_files(index, data, "PIG"),
            None => Ok(None),
        }
    }

    /// The index and data files of the PIG ISAM of the volume `base` (`x_MakeFilenames` with
    /// `file_ext_char` `p`).
    pub fn pig_file_names(
        base: &Path,
        seqtype: SeqType,
    ) -> (std::path::PathBuf, std::path::PathBuf) {
        (
            with_extension(base, seqtype, "pi"),
            with_extension(base, seqtype, "pd"),
        )
    }

    /// The ISAM of `index` and `data`, an index of the numbers `what` names.
    fn from_files(index: DbFile, data: DbFile, what: &str) -> Result<Option<Self>, SeqDbError> {
        let bad = |message: String| SeqDbError::Unsupported {
            path: index.path().to_path_buf(),
            message,
        };
        let crashes = |message: String| SeqDbError::NcbiCrashes {
            path: index.path().to_path_buf(),
            message,
        };
        if index.len() < HEADER_BYTES {
            return Err(crashes(format!(
                "a {what} index file of {} bytes",
                index.len()
            )));
        }
        let header = index.read_range(0, HEADER_BYTES)?;
        let word = |at: usize| be_i32(&header[at * 4..]);
        let version = word(0);
        let isam_type = word(1);
        if version != ISAM_VERSION || (isam_type != NUMERIC && isam_type != NUMERIC_LONG_ID) {
            return Err(crashes(format!(
                "a {what} index file of version {version} and type {isam_type}"
            )));
        }
        let (num_terms, num_samples, page_size) = (word(3), word(4), word(5));
        if num_terms < 1 || num_samples < 1 || page_size <= MEMORY_ONLY_PAGE_SIZE {
            return Err(bad(format!(
                "a {what} index file of {num_terms} terms, {num_samples} samples and pages of \
                 {page_size} terms"
            )));
        }
        if i64::from(word(2)) != data.len() as i64 {
            return Ok(None);
        }
        let mut isam = NumericIsam {
            index,
            data,
            long_id: isam_type == NUMERIC_LONG_ID,
            num_terms,
            num_samples,
            page_size,
            first_key: 0,
            last_key: 0,
        };
        isam.first_key = isam.row(&isam.data, 0)?.0;
        let (start, count) = isam.page(num_samples - 1);
        isam.last_key = isam
            .row(&isam.data, row_offset(isam.term_size(), start, count - 1))?
            .0;
        Ok(Some(isam))
    }

    /// `m_TermSize`.
    fn term_size(&self) -> u64 {
        if self.long_id {
            12
        } else {
            8
        }
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1253-1265
    /// The key and the OID of the row at `offset` of `file`.
    fn row(&self, file: &DbFile, offset: u64) -> Result<(i64, i32), SeqDbError> {
        let end = offset.saturating_add(self.term_size());
        let bytes = file.read_range(offset, end)?;
        if self.long_id {
            let mut key = [0u8; 8];
            key.copy_from_slice(&bytes[..8]);
            Ok((i64::from_be_bytes(key), be_i32(&bytes[8..])))
        } else {
            let key = u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
            Ok((i64::from(key), be_i32(&bytes[4..])))
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:123-137
    /// The first row and the number of rows of the page of `sample`.
    fn page(&self, sample: i32) -> (i64, i64) {
        let start = i64::from(sample) * i64::from(self.page_size);
        let count = if sample + 1 == self.num_samples {
            i64::from(self.num_terms) - start
        } else {
            i64::from(self.page_size)
        };
        (start, count.max(1))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:153-156,162-231
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:429-481
    /// `IdToOid`: the OID the file gives the key `number`.
    pub fn id_to_oid(&self, number: i64) -> Result<Option<i32>, SeqDbError> {
        if number < self.first_key || number > self.last_key {
            return Ok(None);
        }
        let term_size = self.term_size();
        let (mut start, mut stop) = (0i32, self.num_samples - 1);
        let mut sample = 0i32;
        while stop >= start {
            sample = (((stop as u32).wrapping_add(start as u32)) >> 1) as i32;
            let offset = KEY_SAMPLE_OFFSET + term_size * sample as u64;
            let (key, oid) = self.row(&self.index, offset)?;
            if key == number {
                return Ok(Some(oid));
            }
            if number < key {
                sample -= 1;
                stop = sample;
            } else {
                start = sample + 1;
            }
        }
        if sample < 0 || sample >= self.num_samples {
            return Ok(None);
        }
        let (page_start, count) = self.page(sample);
        let (mut first, mut last) = (page_start, page_start + count - 1);
        while first <= last {
            let current = (first + last) / 2;
            let (key, oid) = self.row(&self.data, row_offset(term_size, current, 0))?;
            match key.cmp(&number) {
                std::cmp::Ordering::Greater => last = current - 1,
                std::cmp::Ordering::Less => first = current + 1,
                std::cmp::Ordering::Equal => return Ok(Some(oid)),
            }
        }
        Ok(None)
    }
}

/// The byte offset of a row; past every file when a damaged header gives a row beyond 64 bits.
fn row_offset(term_size: u64, start: i64, index: i64) -> u64 {
    ((start + index) as u64).saturating_mul(term_size)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Writes the files of `rows` (sorted) in the given key width; the index samples every
    /// `page`-th row.
    fn write_isam(
        dir: &Path,
        rows: &[(i64, i32)],
        long_id: bool,
        page: usize,
    ) -> std::path::PathBuf {
        let base = dir.join("v");
        let mut data = Vec::new();
        let mut samples = Vec::new();
        for (at, &(key, oid)) in rows.iter().enumerate() {
            let mut row = Vec::new();
            if long_id {
                row.extend_from_slice(&key.to_be_bytes());
            } else {
                row.extend_from_slice(&(key as u32).to_be_bytes());
            }
            row.extend_from_slice(&oid.to_be_bytes());
            if at % page == 0 {
                samples.extend_from_slice(&row);
            }
            data.extend_from_slice(&row);
        }
        let term = if long_id { 12 } else { 8 };
        let mut index = Vec::new();
        for word in [
            ISAM_VERSION,
            if long_id { NUMERIC_LONG_ID } else { NUMERIC },
            data.len() as i32,
            rows.len() as i32,
            (samples.len() / term) as i32,
            page as i32,
            0,
            0,
            0,
        ] {
            index.extend_from_slice(&word.to_be_bytes());
        }
        index.extend_from_slice(&samples);
        index.extend_from_slice(&vec![0xff; term - 4]);
        index.extend_from_slice(&0i32.to_be_bytes());
        std::fs::write(with_extension(&base, SeqType::Protein, "ni"), index).unwrap();
        std::fs::write(with_extension(&base, SeqType::Protein, "nd"), data).unwrap();
        base
    }

    fn temp_dir(stem: &str) -> std::path::PathBuf {
        let dir = crate::utils::test_dir::temp_dir().join(format!(
            "losat_isam_num_{stem}_{}",
            crate::utils::test_dir::process_id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn every_row_is_found_in_both_key_widths() {
        for long_id in [false, true] {
            let dir = temp_dir(if long_id { "long" } else { "short" });
            let step = if long_id { 1i64 << 33 } else { 7 };
            let rows: Vec<(i64, i32)> = (1..=23).map(|i| (i * step, i as i32 + 100)).collect();
            let base = write_isam(&dir, &rows, long_id, 4);
            let isam = NumericIsam::open_volume(&base, SeqType::Protein)
                .unwrap()
                .unwrap();
            for &(key, oid) in &rows {
                assert_eq!(isam.id_to_oid(key).unwrap(), Some(oid), "key {key}");
                assert_eq!(isam.id_to_oid(key + 1).unwrap(), None, "key {}", key + 1);
                assert_eq!(isam.id_to_oid(key - 1).unwrap(), None, "key {}", key - 1);
            }
            assert_eq!(isam.id_to_oid(-5).unwrap(), None);
            assert_eq!(isam.id_to_oid(i64::MAX).unwrap(), None);
            std::fs::remove_dir_all(&dir).unwrap();
        }
    }

    #[test]
    fn missing_files_are_silent_and_a_bad_header_is_refused() {
        let dir = temp_dir("bad");
        let base = write_isam(&dir, &[(5, 1), (9, 2)], false, 256);
        let nd = with_extension(&base, SeqType::Protein, "nd");
        let ni = with_extension(&base, SeqType::Protein, "ni");
        let data = std::fs::read(&nd).unwrap();
        std::fs::write(&nd, &data[..8]).unwrap();
        assert!(NumericIsam::open_volume(&base, SeqType::Protein)
            .unwrap()
            .is_none());
        std::fs::remove_file(&nd).unwrap();
        assert!(NumericIsam::open_volume(&base, SeqType::Protein)
            .unwrap()
            .is_none());
        std::fs::write(&nd, &data).unwrap();
        let mut index = std::fs::read(&ni).unwrap();
        index[3] = 2;
        std::fs::write(&ni, &index).unwrap();
        assert!(matches!(
            NumericIsam::open_volume(&base, SeqType::Protein),
            Err(error @ SeqDbError::NcbiCrashes { .. }) if error.ncbi_crashes()
                && error.to_string().contains(": unsupported database: a GI index file of version 2 and type ")
                && error.to_string().ends_with(" (NCBI crashes on it)")
        ));
        std::fs::write(&ni, &index[..20]).unwrap();
        assert!(matches!(
            NumericIsam::open_volume(&base, SeqType::Protein),
            Err(SeqDbError::NcbiCrashes { .. })
        ));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
