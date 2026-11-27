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
/// Ten times the size of an `Int4`: the header `x_InitSearch` needs.
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

    /// The index and data files of the volume's ISAM of `kind` (`x_MakeFilenames` with
    /// `file_ext_char` `n`, `t` or `p`).
    pub(crate) fn list_file_names(
        base: &Path,
        seqtype: SeqType,
        kind: ListIsam,
    ) -> (std::path::PathBuf, std::path::PathBuf) {
        match kind {
            ListIsam::Gi => (
                with_extension(base, seqtype, "ni"),
                with_extension(base, seqtype, "nd"),
            ),
            ListIsam::Ti => (
                with_extension(base, seqtype, "ti"),
                with_extension(base, seqtype, "td"),
            ),
            ListIsam::Pig => Self::pig_file_names(base, seqtype),
        }
    }

    /// The ISAM of `index` and `data`, an index of the numbers `what` names: `None` where the
    /// data file has another length than the index says (`NumericIsam::open_files`).
    fn from_files(index: DbFile, data: DbFile, what: &str) -> Result<Option<Self>, SeqDbError> {
        Ok(match Self::open_files(index, data, what)? {
            Opened::Ready(isam) => Some(isam),
            Opened::Uninitialized => None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:136-145
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:161-175
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:214-224
    // The lookups open the volume's ISAM only for a volume with OIDs: the GIs' as the single
    // lookups do (a failure to open is caught), the TIs' and PIGs' where both files are regular
    // files (`IndexExists`), whose constructor's errors are not caught.
    /// The volume's ISAM of `kind` for the lookups of an id list (`x_OpenGiFile`,
    /// `x_OpenTiFile`, `x_OpenPigFile`), and the trace-id ISAM of a single lookup: `None` where
    /// NCBI has none ("... list specified but no ISAM file found").
    pub(crate) fn open_for_list(
        base: &Path,
        seqtype: SeqType,
        kind: ListIsam,
        num_oids: u32,
    ) -> Result<Option<Opened>, SeqDbError> {
        if num_oids == 0 {
            return Ok(None);
        }
        match kind {
            ListIsam::Gi => {
                let (Ok(index), Ok(data)) = (
                    DbFile::open(&with_extension(base, seqtype, "ni")),
                    DbFile::open(&with_extension(base, seqtype, "nd")),
                ) else {
                    return Ok(None);
                };
                Self::open_files(index, data, "GI").map(Some)
            }
            ListIsam::Ti => {
                let index_path = with_extension(base, seqtype, "ti");
                let data_path = with_extension(base, seqtype, "td");
                match super::isam::open_if_exists(&index_path, &data_path)? {
                    Some((index, data)) => Self::open_files(index, data, "TI").map(Some),
                    None => Ok(None),
                }
            }
            ListIsam::Pig => {
                let (index_path, data_path) = Self::pig_file_names(base, seqtype);
                match super::isam::open_if_exists(&index_path, &data_path)? {
                    Some((index, data)) => Self::open_files(index, data, "PIG").map(Some),
                    None => Ok(None),
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1102-1170
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:1461-1520
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1492-1517
    // The constructor runs `x_FindIndexBounds` whatever `x_InitSearch` returned, with the counts
    // of the index header that `x_InitSearch` read before it compared the data file's length: it
    // reads the first row of the data file and the row before `m_NumTerms`, through the mapping
    // of the whole file without a bound (`GetFileDataPtr`). Inside the last memory page of the
    // mapping it reads the zeros after the end of the file (`file::mapped_length`); past that
    // page it reads memory that is not the file's, which LOSAT refuses. So the bounds are read in
    // both states, and a data file of another length than the index says leaves the ISAM
    // uninitialised only where they are in the mapping.
    /// The ISAM of `index` and `data`, an index of the numbers `what` names
    /// (`CSeqDBIsam::CSeqDBIsam`): uninitialized where the data file has another length than the
    /// index says (`x_InitSearch` fails). LOSAT's refusal where NCBI's reads of the first and
    /// last rows pass the memory page after the data file.
    fn open_files(index: DbFile, data: DbFile, what: &str) -> Result<Opened, SeqDbError> {
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
        // The last page's start cancels out: its last element is the row before `m_NumTerms`,
        // also where the samples' pages would end before it or start after it.
        let last = row_offset(isam.term_size(), 0, i64::from(num_terms) - 1);
        let data_len = isam.data.len();
        let mapped = super::file::mapped_length(data_len);
        if last.saturating_add(isam.term_size()) > mapped {
            return Err(SeqDbError::NcbiCrashes {
                path: isam.data.path().to_path_buf(),
                message: format!(
                    "NCBI reads row {} of the {what} data file (its last of {num_terms} terms), \
                     past the memory page after its {data_len} bytes",
                    num_terms - 1
                ),
            });
        }
        isam.first_key = isam.bound_key(0)?;
        isam.last_key = isam.bound_key(last)?;
        if i64::from(word(2)) != data_len as i64 {
            return Ok(Opened::Uninitialized);
        }
        Ok(Opened::Ready(isam))
    }

    /// The key of the row at `offset` of the data file as NCBI's mapping holds it: the bytes past
    /// the end of the file are zeros (`open_files` checked the memory page).
    fn bound_key(&self, offset: u64) -> Result<i64, SeqDbError> {
        let len = self.data.len();
        if offset.saturating_add(self.term_size()) <= len {
            return Ok(self.row(&self.data, offset)?.0);
        }
        let mut bytes = if offset < len {
            self.data.read_range(offset, len)?
        } else {
            Vec::new()
        };
        bytes.resize(self.term_size() as usize, 0);
        Ok(self.key_of(&bytes))
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
        let oid = if self.long_id {
            &bytes[8..]
        } else {
            &bytes[4..]
        };
        Ok((self.key_of(&bytes), be_i32(oid)))
    }

    /// The key at the start of the row `bytes` (`x_GetNumericKey`).
    fn key_of(&self, bytes: &[u8]) -> i64 {
        if self.long_id {
            let mut key = [0u8; 8];
            key.copy_from_slice(&bytes[..8]);
            i64::from_be_bytes(key)
        } else {
            i64::from(u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
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

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:498-516
    /// `x_LoadIndex`: the key of every sample.
    fn sample_keys(&self) -> Result<Vec<i64>, SeqDbError> {
        let term_size = self.term_size();
        (0..self.num_samples as u64)
            .map(|sample| {
                self.row(&self.index, KEY_SAMPLE_OFFSET + term_size * sample)
                    .map(|(key, _)| key)
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:548-623
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:51-115
    // The list and the keys are walked together: the list (sorted) is advanced past the ids below
    // a key and past the ids translated already, and every id equal to a key takes that row's
    // OID. So an id that an earlier volume translated keeps its OID, and of rows with the same
    // key the walk reaches the first.
    /// `IdsToOids` of a numeric list over the volume's OIDs (`x_TranslateGiList`):
    /// sets `oid` of every entry of `ids` (sorted by key) without one whose key is a key of the
    /// volume, to the OID of its row plus `vol_start`.
    pub(crate) fn ids_to_oids(&self, vol_start: u32, ids: &mut [NumOid]) -> Result<(), SeqDbError> {
        let gilist_size = ids.len();
        if gilist_size == 0 {
            return Ok(());
        }
        debug_assert!(ids.windows(2).all(|w| w[0].key <= w[1].key), "ids sorted");
        let sample_keys = self.sample_keys()?;
        let num_samples = self.num_samples as usize;
        let term_size = self.term_size();
        let mut gilist_index = 0;
        let mut sample_index: i64 = 0;
        while gilist_index < gilist_size && (sample_index as usize) < num_samples {
            advance_list(ids, &mut gilist_index, sample_keys[sample_index as usize]);
            if gilist_index >= gilist_size {
                break;
            }
            advance_keys(
                &sample_keys,
                &mut sample_index,
                num_samples,
                ids[gilist_index].key,
            );
            if sample_index < 0 {
                return Err(SeqDbError::Format {
                    path: self.index.path().to_path_buf(),
                    message: "numeric ISAM sample keys are not in order".to_string(),
                });
            }
            let (start, count) = self.page(sample_index as i32);
            let num_keys = count as usize;
            let mut keys = Vec::with_capacity(num_keys);
            let mut vals = Vec::with_capacity(num_keys);
            for at in 0..count {
                let (key, oid) = self.row(&self.data, row_offset(term_size, start, at))?;
                keys.push(key);
                vals.push(oid);
            }
            let mut index: i64 = 0;
            while gilist_index < gilist_size && (index as usize) < num_keys {
                advance_keys(&keys, &mut index, num_keys, ids[gilist_index].key);
                if index < 0 {
                    return Err(SeqDbError::Format {
                        path: self.data.path().to_path_buf(),
                        message: "numeric ISAM keys are not in order".to_string(),
                    });
                }
                let at = index as usize;
                set_translation(ids, &mut gilist_index, keys[at], vals[at], vol_start);
                index += 1;
                if index as usize >= num_keys {
                    break;
                }
                let at = index as usize;
                advance_list(ids, &mut gilist_index, keys[at]);
                set_translation(ids, &mut gilist_index, keys[at], vals[at], vol_start);
            }
            sample_index += 1;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbisam.cpp:219-302
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1428-1457
    // Every data page is walked in order; each row's key is tested against the sorted negative
    // list with parabolic search, and each row whose OID is below `vol_end` marks `oid + vol_start`
    // visible when its key is listed and included when it is not.
    /// `x_SearchNegativeMulti`: walks every row of the numeric ISAM and records each OID
    /// (plus `vol_start`, where the row's OID is below `vol_end`) in `out` as visible when its key
    /// is in `sorted_keys` and included otherwise.
    pub(crate) fn negative_scan(
        &self,
        vol_start: u32,
        vol_end: u32,
        sorted_keys: &[i64],
        out: &mut super::isam::NegativeOids,
    ) -> Result<(), SeqDbError> {
        debug_assert!(
            sorted_keys.windows(2).all(|w| w[0] <= w[1]),
            "sorted_keys sorted"
        );
        let term_size = self.term_size();
        let term_len = term_size as usize;
        let oid_off = if self.long_id { 8 } else { 4 };
        let overflow = || SeqDbError::Format {
            path: self.data.path().to_path_buf(),
            message: "numeric ISAM OID plus the start of the volume overflows".to_string(),
        };
        let mut list_index = 0usize;
        for sample_index in 0..self.num_samples {
            let (start, count) = self.page(sample_index);
            let begin = row_offset(term_size, start, 0);
            let end = row_offset(term_size, start, count);
            let page = self.data.read_range(begin, end)?;
            for row in page.chunks_exact(term_len) {
                let key = self.key_of(row);
                let val = be_i32(&row[oid_off..]);
                let found = if list_index < sorted_keys.len() {
                    find_in_negative_list(sorted_keys, &mut list_index, key)
                } else {
                    false
                };
                if i64::from(val) < i64::from(vol_end) {
                    let oid = (val as u32).checked_add(vol_start).ok_or_else(overflow)?;
                    if i32::try_from(oid).is_err() {
                        return Err(overflow());
                    }
                    if found {
                        out.add_visible_oid(oid);
                    } else {
                        out.add_included_oid(oid);
                    }
                }
            }
        }
        Ok(())
    }
}

/// A numeric ISAM opened for an id list's lookups or a single TI's ([`NumericIsam::open_for_list`]).
#[derive(Debug)]
pub(crate) enum Opened {
    Ready(NumericIsam),
    /// `m_Initialized` false: the list's lookups throw "Unable to use ISAM index in batch mode",
    /// a single lookup finds nothing.
    Uninitialized,
}

/// The numeric ISAMs that an id list's lookups read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ListIsam {
    /// `.pni`/`.pnd`, `.nni`/`.nnd`.
    Gi,
    /// `.pti`/`.ptd`, `.nti`/`.ntd`.
    Ti,
    /// `.ppi`/`.ppd`, `.npi`/`.npd`.
    Pig,
}

/// One entry of a numeric id list (`CSeqDBGiList::SGiOid`, `STiOid`): the id and its OID, set by
/// [`NumericIsam::ids_to_oids`] (NCBI: -1 while unknown).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NumOid {
    pub key: i64,
    pub oid: Option<u32>,
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:51-75
/// `s_AdvanceGiList`: moves `gi_index` past the ids below `key`, then past those translated.
fn advance_list(gis: &[NumOid], gi_index: &mut usize, key: i64) {
    let gis_size = gis.len();
    while *gi_index < gis_size && gis[*gi_index].key < key {
        *gi_index += 1;
        let mut jump = 2;
        while *gi_index + jump < gis_size && gis[*gi_index + jump].key < key {
            *gi_index += jump;
            jump *= 2;
        }
    }
    while *gi_index < gis_size && gis[*gi_index].oid.is_some() {
        *gi_index += 1;
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:77-99
/// `s_AdvanceKeyList`: moves `index` to the last key that is `<= target`; -1 when none is.
fn advance_keys(keys: &[i64], index: &mut i64, num_keys: usize, target: i64) {
    let num_keys = num_keys as i64;
    while *index < num_keys && keys[*index as usize] <= target {
        *index += 1;
        let mut jump = 2;
        while *index + jump < num_keys && keys[(*index + jump) as usize] <= target {
            *index += jump;
            jump *= 2;
        }
    }
    *index -= 1;
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:101-115
/// `s_SetTranslation`: the OID of every entry from `gi_index` on whose key is `key`.
fn set_translation(gis: &mut [NumOid], gi_index: &mut usize, key: i64, value: i32, vol_start: u32) {
    while *gi_index < gis.len() && gis[*gi_index].key == key {
        gis[*gi_index].oid = Some((value as u32).wrapping_add(vol_start));
        *gi_index += 1;
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbisam.hpp:1428-1457
/// `x_FindInNegativeList`: advances `index` past entries below `key` and reports whether `key` is
/// at `index`.
fn find_in_negative_list(ids: &[i64], index: &mut usize, key: i64) -> bool {
    let ids_size = ids.len();
    while *index < ids_size && ids[*index] < key {
        *index += 1;
        let mut jump = 2;
        while *index + jump < ids_size && ids[*index + jump] < key {
            *index += jump;
            jump += jump;
        }
    }
    *index < ids_size && ids[*index] == key
}

/// The byte offset of a row; past every file when a damaged header gives a row beyond 64 bits.
fn row_offset(term_size: u64, start: i64, index: i64) -> u64 {
    ((start + index) as u64).saturating_mul(term_size)
}

#[cfg(test)]
#[path = "isam_num_tests.rs"]
mod tests;
