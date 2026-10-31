//! The numeric ISAM of a volume (`.pni`/`.pnd`, `.nni`/`.nnd`): NCBI's `CSeqDBIsam` for the
//! lookup of a GI (`IdToOid`). The string ISAM is in [`super::isam`].
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
    // ```c++
    // CSeqDBVol::x_OpenGiFile(void) const{
    //     CFastMutexGuard mtx_gurad(m_MtxGi);
    //     if (m_IsamGi.NotEmpty()) {
    //     	m_IsamGi->AddReference();
    //     }
    //     else if ((m_Idx->GetNumOIDs() != 0) && m_HasGis) {
    //     	try {
    //         m_IsamGi = new CSeqDBIsam(m_Atlas, m_VolName, (m_IsAA?'p':'n'), 'n', eGiId);
    //     	} catch(CException &) {
    //     		m_IsamGi.Reset();
    //     		m_HasGis = false;
    //     	}
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:64-118
    // ```c++
    //     TIndx info_needed = 10 * sizeof(Int4);
    //
    //     bool found_index_file =
    //         m_Atlas.GetFileSizeL(m_IndexFname, m_IndexFileLength);
    //
    //     if ((! found_index_file) || (m_IndexFileLength < info_needed)) {
    //         return eWrongFile;
    //     }
    //     ...
    //     Int4 Version = SeqDB_GetStdOrd(& FileInfo[0]);
    //
    //     if (Version != ISAM_VERSION)
    //         return eBadVersion;
    //
    //     Int4 IsamType = SeqDB_GetStdOrd(& FileInfo[1]);
    //
    //     if (IsamType == eNumericLongId && m_Type == eNumeric) {
    //         m_LongId = true;
    //         m_TermSize = 12;
    //         IsamType = eNumeric;
    //     }
    //
    //     if (IsamType != m_Type)
    //         return eBadType;
    //
    //     m_NumTerms    = SeqDB_GetStdOrd(& FileInfo[3]);
    //     m_NumSamples  = SeqDB_GetStdOrd(& FileInfo[4]);
    //     m_PageSize    = SeqDB_GetStdOrd(& FileInfo[5]);
    //     m_MaxLineSize = SeqDB_GetStdOrd(& FileInfo[6]);
    //
    //     if(m_PageSize != MEMORY_ONLY_PAGE_SIZE) {
    //         // Special case of memory-only index
    //         m_DataFileLength = SeqDB_GetStdOrd(& FileInfo[2]);
    //
    //         TIndx disk_file_length(0);
    //         bool found_data_file =
    //             m_Atlas.GetFileSizeL(m_DataFname, disk_file_length);
    //
    //         if ((! found_data_file) || (m_DataFileLength != disk_file_length)) {
    //             return eWrongFile;
    //         }
    //     }
    // ```
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
        let bad = |message: String| SeqDbError::Unsupported {
            path: index.path().to_path_buf(),
            message,
        };
        if index.len() < HEADER_BYTES {
            return Err(bad(format!(
                "a GI index file of {} bytes (NCBI crashes on it)",
                index.len()
            )));
        }
        let header = index.read_range(0, HEADER_BYTES)?;
        let word = |at: usize| be_i32(&header[at * 4..]);
        let version = word(0);
        let isam_type = word(1);
        if version != ISAM_VERSION || (isam_type != NUMERIC && isam_type != NUMERIC_LONG_ID) {
            return Err(bad(format!(
                "a GI index file of version {version} and type {isam_type} (NCBI crashes on it)"
            )));
        }
        let (num_terms, num_samples, page_size) = (word(3), word(4), word(5));
        if num_terms < 1 || num_samples < 1 || page_size <= MEMORY_ONLY_PAGE_SIZE {
            return Err(bad(format!(
                "a GI index file of {num_terms} terms, {num_samples} samples and pages of \
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
    // ```c++
    //     Uint8 x_GetNumericKey(const void *p) {
    //         if (m_LongId)
    //             return((Uint8) SeqDB_GetStdOrd((Uint8 *)p));
    //         else
    //             return((Uint4) SeqDB_GetStdOrd((Uint4 *)p));
    //     }
    //
    //     int x_GetNumericData(const void *p) {
    //         if (m_LongId)
    //             return((int) SeqDB_GetStdOrd(((Uint4 *)p)+2));
    //         else
    //             return((int) SeqDB_GetStdOrd(((Uint4 *)p)+1));
    //     }
    // ```
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
    // ```c++
    // Int4 CSeqDBIsam::x_GetPageNumElements(Int4   sample_num,
    //                                       Int4 * start)
    // {
    //     Int4 num_elements(0);
    //
    //     *start = sample_num * m_PageSize;
    //
    //     if (sample_num + 1 == m_NumSamples) {
    //         num_elements = m_NumTerms - *start;
    //     } else {
    //         num_elements = m_PageSize;
    //     }
    //
    //     return num_elements;
    // }
    // ```
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
    // ```c++
    //     if (x_OutOfBounds(Number)) {
    //         done = true;
    //         return eNotFound;
    //     }
    //     ...
    //     Int4 Start     (0);
    //     Int4 Stop      (m_NumSamples - 1);
    //
    //     while(Stop >= Start) {
    //         SampleNum = ((Uint4)(Stop + Start)) >> 1;
    //         ...
    //         Key = x_GetNumericKey (keydatap);
    //
    //         // If this is an exact match, return the master term number.
    //
    //         if (Key == Number) {
    //             if (Data != NULL) {
    //                 *Data = x_GetNumericData(keydatap);
    //             }
    //             ...
    //             done = true;
    //             return eNoError;
    //         }
    //
    //         // Otherwise, search for the next sample.
    //
    //         if ( Number < Key )
    //             Stop = --SampleNum;
    //         else
    //             Start = SampleNum +1;
    //     }
    //
    //     // If the term is out of range altogether, report not finding it.
    //
    //     if ( (SampleNum < 0) || (SampleNum >= m_NumSamples)) {
    //     ...
    //         done = true;
    //         return eNotFound;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:429-481
    // ```c++
    //     Int4 first = Start;
    //     Int4 last  = Start + NumElements - 1;
    //     ...
    //     // Search the page for the number.
    //     while (first <= last) {
    //         current = (first+last)/2;
    //
    //         Int8 Key = x_GetNumericKey((char *)KeyDataPage + current * m_TermSize);
    //
    //         if (Key > Number) {
    //             last = --current;
    //         } else if (Key < Number) {
    //             first = ++current;
    //         } else {
    //             found = true;
    //             break;
    //         }
    //     }
    // ```
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

fn row_offset(term_size: u64, start: i64, index: i64) -> u64 {
    (start + index) as u64 * term_size
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
        let dir =
            std::env::temp_dir().join(format!("losat_isam_num_{stem}_{}", std::process::id()));
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
            Err(SeqDbError::Unsupported { .. })
        ));
        std::fs::write(&ni, &index[..20]).unwrap();
        assert!(matches!(
            NumericIsam::open_volume(&base, SeqType::Protein),
            Err(SeqDbError::Unsupported { .. })
        ));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
