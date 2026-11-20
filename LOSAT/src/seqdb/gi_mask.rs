//! The GI-indexed masks of a database (`makeblastdb -gi_mask`): the files `<name>.gni` (index),
//! `<name>.gno` (offsets) and `<name>.gnd` or `<name>.NN.gnd` (data) of each name of the alias
//! file's `MASKLIST` (`alias::AliasTree::gi_masks`), and the GI of each sequence that keys them.
//!
//! NCBI reads them with `CSeqDBGiMask` (`seqdbgimask.cpp`, `seqdbgimask.hpp`), and the GI with
//! `CSeqDBImpl::x_GetSeqGI` and the volume's `CSeqDBGiIndex` (the `.nog` file). This module is
//! the port; `mask.rs` dispatches to it where the database has GI masks, and
//! `blastinput::db_mask` and `losat blastdbcmd` read them through `mask.rs`.
//!
//! NCBI maps the files and reads its arrays in place. LOSAT reads them in chunks of 64 KiB on
//! first use (`MappedFile`), kept for the life of the database as NCBI's mapping keeps the pages
//! it touched, so a search that reads the masks of every sequence reads each file once.
//!
//! Errors: NCBI's `CSeqDBException` texts are NCBI's ([`SeqDbError::NotFound`],
//! [`SeqDbError::Format`], [`SeqDbError::Range`], [`SeqDbError::MemoryMap`]). A failed
//! `SEQDB_FILE_ASSERT` (its text names a source line) is [`SeqDbError::Integrity`] in LOSAT's
//! words; a read outside a file, where NCBI reads past its mapping, is [`SeqDbError::Format`] in
//! LOSAT's words.

use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use super::column::Blob;
use super::file::DbFile;
use super::volume::{with_extension, Volume};
use super::{SearchPath, SeqDbError, SeqType};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:124-133
// ```c++
//     string ext_i(".gmi");
//     string ext_o(".gmo");
//     string ext_d(".gmd");
//
//     const Int2 bytetest = 0x0011;
//     const char *ptr = (const char *) &bytetest;
//     if (ptr[0] == 0x11) {
//         // Use small endian instead.
//         ext_i[2] = ext_o[2] = ext_d[2] = 'n';
//     }
// ```
/// The extensions of the index, offset and data files that a host of this byte order reads; their
/// arrays are in the host's order.
const EXTENSIONS: [&str; 3] = if cfg!(target_endian = "little") {
    ["gni", "gno", "gnd"]
} else {
    ["gmi", "gmo", "gmd"]
};

/// The bytes of the fixed fields at the start of an index file (`kFixedFieldBytes`).
const FIXED_FIELD_BYTES: u64 = 32;

/// The size of the chunks in which `MappedFile` reads a file.
const CHUNK: u64 = 1 << 16;

/// A file read in chunks on first use and kept, as NCBI's memory map of it.
#[derive(Debug)]
struct MappedFile {
    file: DbFile,
    chunks: Box<[OnceLock<Box<[u8]>>]>,
}

impl MappedFile {
    fn open(path: &Path) -> Result<Self, SeqDbError> {
        let file = DbFile::open(path)?;
        let count = file.len().div_ceil(CHUNK);
        let chunks = (0..count).map(|_| OnceLock::new()).collect();
        Ok(Self { file, chunks })
    }

    fn len(&self) -> u64 {
        self.file.len()
    }

    fn path(&self) -> &Path {
        self.file.path()
    }

    fn chunk(&self, index: usize) -> Result<&[u8], SeqDbError> {
        if let Some(bytes) = self.chunks[index].get() {
            return Ok(bytes);
        }
        let start = index as u64 * CHUNK;
        let bytes = self
            .file
            .read_range(start, (start + CHUNK).min(self.len()))?;
        Ok(self.chunks[index].get_or_init(|| bytes.into_boxed_slice()))
    }

    /// The bytes at `start..start + out.len()` into `out`. A range outside the file, which NCBI
    /// reads outside its mapping, is a format error in LOSAT's words.
    fn read_at(&self, start: i64, out: &mut [u8]) -> Result<(), SeqDbError> {
        let end = start.saturating_add(out.len() as i64);
        if start < 0 || end as u64 > self.len() {
            return Err(SeqDbError::Format {
                path: self.path().to_path_buf(),
                message: format!(
                    "a GI-mask read of bytes {start}..{end} is outside the file ({} bytes)",
                    self.len()
                ),
            });
        }
        let mut at = start as u64;
        let mut filled = 0;
        while filled < out.len() {
            let chunk = self.chunk((at / CHUNK) as usize)?;
            let from = (at % CHUNK) as usize;
            let take = (chunk.len() - from).min(out.len() - filled);
            out[filled..filled + take].copy_from_slice(&chunk[from..from + take]);
            filled += take;
            at += take as u64;
        }
        Ok(())
    }

    /// The 4 bytes at `at` as a number in the host's order (the arrays of the GI-mask files).
    fn u32_at(&self, at: i64) -> Result<u32, SeqDbError> {
        let mut bytes = [0u8; 4];
        self.read_at(at, &mut bytes)?;
        Ok(u32::from_ne_bytes(bytes))
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:123-130
    // ```c++
    //     const char * GetFileDataPtr(CSeqDBFileMemMap & lease, // commented
    //                            TIndx            start,
    //                            TIndx            end) const
    //     {
    //         _ASSERT(! m_FileName.empty());
    //         SEQDB_FILE_ASSERT(start    <  end);
    //         SEQDB_FILE_ASSERT(m_Length >= end);
    // ```
    /// The two `SEQDB_FILE_ASSERT`s of `GetFileDataPtr(lease, start, end)`, `what` naming the
    /// range in LOSAT's error.
    fn check_range(&self, start: i64, end: i64, what: &str) -> Result<(), SeqDbError> {
        if start < end && end as u64 <= self.len() && end >= 0 {
            return Ok(());
        }
        Err(SeqDbError::Integrity {
            path: self.path().to_path_buf(),
            message: format!(
                "{what} (bytes {start}..{end}) is empty or passes the end of the file ({} bytes)",
                self.len()
            ),
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:204-245
// ```c++
//     template<class T> static bool s_BinarySearch(const T *keys,
//                                const int  n,
//                                const T    key,
//                                int       &idx)
//     {
//         int lower(0), upper(n - 1);
//
//         if (key > keys[upper] || key < keys[lower]) {
//             // out of range
//             idx = -1;
//             return false;
//         }
//
//         if (key == keys[upper]) {
//             idx = upper;
//             return true;
//         }
//
//         if (key == keys[lower]) {
//             idx = lower;
//             return true;
//         }
//
//         idx = (lower + upper) / 2;
//
//         while (idx != lower) {
//             if (key > keys[idx]) {
//                 lower = idx;
//                 idx = (lower + upper) / 2;
//             }
//             else if (key < keys[idx]) {
//                 upper = idx;
//                 idx = (lower + upper) / 2;
//             }
//             else {
//                 // value found
//                 return true;
//             }
//         }
//         // value not found
//         return false;
//     }
// ```
// The port reads each key where NCBI does (`key_at`), so a list with equal GIs (a version 4
// database, survey D7) gives the entry NCBI's search lands on. The indexes are 64-bit, where
// NCBI's `int` sum could overflow.
/// `s_BinarySearch` over the `n` keys that `key_at` reads: whether `key` was found, and the index
/// (the entry found; otherwise -1 for a key outside the keys, else the entry before it).
fn binary_search(
    n: i64,
    key: u32,
    key_at: impl Fn(i64) -> Result<u32, SeqDbError>,
) -> Result<(bool, i64), SeqDbError> {
    let (mut lower, mut upper) = (0i64, n - 1);
    if key > key_at(upper)? || key < key_at(lower)? {
        return Ok((false, -1));
    }
    if key == key_at(upper)? {
        return Ok((true, upper));
    }
    if key == key_at(lower)? {
        return Ok((true, lower));
    }
    let mut idx = (lower + upper) / 2;
    while idx != lower {
        let at = key_at(idx)?;
        if key > at {
            lower = idx;
            idx = (lower + upper) / 2;
        } else if key < at {
            upper = idx;
            idx = (lower + upper) / 2;
        } else {
            return Ok((true, idx));
        }
    }
    Ok((false, idx))
}

/// `setw(width) << left << text`: `text`, then spaces up to `width` bytes.
fn push_padded(out: &mut String, text: &str, width: usize) {
    out.push_str(text);
    for _ in text.len()..width {
        out.push(' ');
    }
}

/// The opened files of one algorithm (`x_Open`, `x_ReadFields`).
#[derive(Debug)]
pub(super) struct GiMaskFiles {
    index: MappedFile,
    offsets: MappedFile,
    data: Vec<MappedFile>,
    gi_size: i32,
    offset_size: i32,
    page_size: i32,
    num_index: i32,
    num_gis: i32,
    index_start: i32,
    desc: Vec<u8>,
}

/// The fields of an index file (`x_ReadFields`) and the number of data volumes.
struct Fields {
    num_vols: i32,
    gi_size: i32,
    offset_size: i32,
    page_size: i32,
    num_index: i32,
    num_gis: i32,
    index_start: i32,
    desc: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:200-251
// ```c++
// void CSeqDBGiMask::x_ReadFields(void)//
// {
//     const int kFixedFieldBytes = 32;
//
//     //m_Atlas.Lock(locked);
//
//     // First, get the 32 bytes of fields that we know exist.
//
//     CBlastDbBlob header;
//     s_GetFileRange(0, kFixedFieldBytes, m_IndexFile, m_IndexLease, header);
//
//     int fmt_version = header.ReadInt4();
//
//     if (fmt_version != 1) {
//         NCBI_THROW(CSeqDBException,
//                    eFileErr,
//                    "Gi-mask file uses unknown format_version.");
//     }
//
//     m_NumVols    = header.ReadInt4();
//     m_GiSize     = header.ReadInt4();
//     m_OffsetSize = header.ReadInt4();
//     m_PageSize   = header.ReadInt4();
//
//     m_NumIndex   = header.ReadInt4();
//     m_NumGi      = header.ReadInt4();
//
//     m_IndexStart = header.ReadInt4();
//
//     SEQDB_FILE_ASSERT(m_IndexStart >= 0);
//     SEQDB_FILE_ASSERT(m_IndexFile.GetFileLength() >= m_IndexStart);
//
//     // Now we know how long the header actually is, so expand the blob
//     // to reference the whole thing.  (The memory lease should already
//     // hold the data, so this will just adjust a few integer fields.)
//
//     s_GetFileRange(0, m_IndexStart, m_IndexFile, m_IndexLease, header);
//
//     // Get string type header fields.
//
//     m_Desc  = header.ReadString (kStringFmt);
//     m_Date  = header.ReadString (kStringFmt);
//
//     SEQDB_FILE_ASSERT(m_Desc.size());
//     SEQDB_FILE_ASSERT(m_Date.size());
//
//     // Map the index file
//     TIndx begin = m_IndexStart;
//     TIndx end = begin + m_NumIndex * (m_GiSize + m_OffsetSize);
//     m_GiIndex = (const Uint4 *)
//               m_IndexFile.GetFileDataPtr(m_IndexLease, begin, end);
// }
// ```
// `ReferTo` keeps the blob's read offset, so the strings are read from byte 32 of the header.
// The sums are 64-bit, where NCBI's `int` product could overflow.
/// `x_ReadFields` of the index file `index`.
fn read_fields(index: &MappedFile) -> Result<Fields, SeqDbError> {
    let format = |message: String| SeqDbError::Format {
        path: index.path().to_path_buf(),
        message,
    };
    index.check_range(
        0,
        FIXED_FIELD_BYTES as i64,
        "the fixed fields of the GI-mask index",
    )?;
    let mut header = [0u8; FIXED_FIELD_BYTES as usize];
    index.read_at(0, &mut header)?;
    let mut blob = Blob::new(&header, 0);
    let mut field = || blob.read_i32().map_err(format);
    if field()? != 1 {
        return Err(format(
            "Gi-mask file uses unknown format_version.".to_string(),
        ));
    }
    let num_vols = field()?;
    let gi_size = field()?;
    let offset_size = field()?;
    let page_size = field()?;
    let num_index = field()?;
    let num_gis = field()?;
    let index_start = field()?;
    let integrity = |message: String| SeqDbError::Integrity {
        path: index.path().to_path_buf(),
        message,
    };
    if index_start < 0 {
        return Err(integrity(format!(
            "the GI-mask index states a negative start of its arrays ({index_start})"
        )));
    }
    if index.len() < index_start as u64 {
        return Err(integrity(format!(
            "the GI-mask index has {} bytes, fewer than the start of its arrays ({index_start})",
            index.len()
        )));
    }
    index.check_range(0, i64::from(index_start), "the header of the GI-mask index")?;
    let mut head = vec![0u8; index_start as usize];
    index.read_at(0, &mut head)?;
    let mut blob = Blob::new(&head, FIXED_FIELD_BYTES as i64);
    let desc = blob.read_string_var().map_err(format)?.to_vec();
    let date = blob.read_string_var().map_err(format)?;
    if desc.is_empty() || date.is_empty() {
        return Err(integrity(
            "the GI-mask index has an empty description or date".to_string(),
        ));
    }
    let begin = i64::from(index_start);
    let end = begin + i64::from(num_index) * (i64::from(gi_size) + i64::from(offset_size));
    index.check_range(begin, end, "the arrays of the GI-mask index")?;
    Ok(Fields {
        num_vols,
        gi_size,
        offset_size,
        page_size,
        num_index,
        num_gis,
        index_start,
        desc,
    })
}

impl GiMaskFiles {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:137-178
    // ```c++
    //     try {
    //         CSeqDB_Path fn_i(SeqDB_ResolveDbPath(m_MaskNames[algo_id] + ext_i));
    //         CSeqDB_Path fn_o(SeqDB_ResolveDbPath(m_MaskNames[algo_id] + ext_o));
    //
    //         bool found_i = m_IndexFile.Open(fn_i);
    //         bool found_o = m_OffsetFile.Open(fn_o);
    //
    //         if (! (found_i && found_o)) {
    //             NCBI_THROW(CSeqDBException, eFileErr,
    //                        "Could not open gi-mask index files.");
    //         }
    //         m_IndexLease.Init(fn_i.GetPathS());
    //         m_OffsetLease.Init(fn_o.GetPathS());
    //
    //         m_AlgoId = algo_id;
    //         x_ReadFields();
    //
    //         if (m_NumVols == 1) {
    //             m_DataFile.push_back(new CSeqDBRawFile(m_Atlas));
    //             m_DataLease.push_back(new CSeqDBFileMemMap(m_Atlas));
    //             CSeqDB_Path fn(SeqDB_ResolveDbPath(m_MaskNames[algo_id] + ext_d));
    //             bool found = m_DataFile[0]->Open(fn);
    //             m_DataLease[0]->Init(fn.GetPathS());
    //             if (! found) {
    //                 NCBI_THROW(CSeqDBException, eFileErr,
    //                        "Could not open gi-mask data file.");
    //             }
    //         } else {
    //             for (int vol=0; vol<m_NumVols; ++vol) {
    //                 m_DataFile.push_back(new CSeqDBRawFile(m_Atlas));
    //                 m_DataLease.push_back(new CSeqDBFileMemMap(m_Atlas));
    //                 ostringstream fnd;
    //                 fnd << m_MaskNames[algo_id] << "." << vol/10 << vol%10 << ext_d;
    //                 CSeqDB_Path fn(SeqDB_ResolveDbPath(fnd.str()));
    //                 bool found = m_DataFile[vol]->Open(fn);
    //                 m_DataLease[vol]->Init(fn.GetPathS());
    //                 if (! found) {
    //                     NCBI_THROW(CSeqDBException, eFileErr,
    //                         "Could not open gi-mask data files.");
    //                 }
    //             }
    //         }
    // ```
    // The files are looked for along the search path (the current directory, `BLASTDB`, the
    // registry), never beside the alias file (survey D5). A data file that is not found maps the
    // empty name before `found` is tested (`m_DataLease[vol]->Init("")`): the corelib's
    // "To be memory mapped the file must exist: ''" (`CSeqDBException::eMemoryMappingFailure`,
    // seqdbatlas.hpp:709-719), so "Could not open gi-mask data file(s)." is never thrown (survey
    // D8). The data files of an earlier algorithm stay in NCBI's vectors, but `[vol]` is opened
    // again, so each algorithm reads its own files.
    /// `x_Open` of the mask `name`: its index, offset and data files.
    fn open(name: &str, search: &SearchPath) -> Result<Self, SeqDbError> {
        let index_path = search.resolve_file(&format!("{name}.{}", EXTENSIONS[0]));
        let offset_path = search.resolve_file(&format!("{name}.{}", EXTENSIONS[1]));
        let (Some(index_path), Some(offset_path)) = (index_path, offset_path) else {
            return Err(SeqDbError::NotFound(
                "Could not open gi-mask index files.".to_string(),
            ));
        };
        let index = MappedFile::open(Path::new(&index_path))?;
        let offsets = MappedFile::open(Path::new(&offset_path))?;
        let fields = read_fields(&index)?;
        let mut data = Vec::new();
        for vol in 0..fields.num_vols.max(0) {
            let file = if fields.num_vols == 1 {
                format!("{name}.{}", EXTENSIONS[2])
            } else {
                format!("{name}.{}{}.{}", vol / 10, vol % 10, EXTENSIONS[2])
            };
            let Some(path) = search.resolve_file(&file) else {
                return Err(SeqDbError::MemoryMap(
                    "To be memory mapped the file must exist: ''".to_string(),
                ));
            };
            data.push(MappedFile::open(Path::new(&path))?);
        }
        Ok(Self {
            index,
            offsets,
            data,
            gi_size: fields.gi_size,
            offset_size: fields.offset_size,
            page_size: fields.page_size,
            num_index: fields.num_index,
            num_gis: fields.num_gis,
            index_start: fields.index_start,
            desc: fields.desc,
        })
    }

    /// The `(vol, off)` of the record of `key`, or `None` without one: the index, then the page
    /// of the offset file (`GetMaskData`, quoted at [`GiMaskFiles::ranges_into`]).
    fn locate(&self, key: u32) -> Result<Option<(i32, i32)>, SeqDbError> {
        let index_start = i64::from(self.index_start);
        let num_index = i64::from(self.num_index);
        let (found, page) =
            binary_search(num_index, key, |i| self.index.u32_at(index_start + 4 * i))?;
        let (file, at) = if found {
            (&self.index, index_start + 4 * (num_index + page * 2))
        } else {
            if page == -1 {
                return Ok(None);
            }
            let page_size = i64::from(self.page_size);
            let num_gis = i64::from(self.num_gis);
            let mut pagesize = page_size;
            if page * page_size + pagesize > num_gis {
                pagesize = num_gis - page * page_size;
            }
            let entry = i64::from(self.gi_size) + i64::from(self.offset_size);
            let begin = page * page_size * entry;
            let end = begin + pagesize * entry;
            self.offsets
                .check_range(begin, end, "a page of the GI-mask offsets")?;
            let (found, idx) =
                binary_search(pagesize, key, |i| self.offsets.u32_at(begin + 4 * i))?;
            if !found {
                return Ok(None);
            }
            (&self.offsets, begin + 4 * (pagesize + idx * 2))
        };
        // NCBI assigns the `Uint4` values to `int`s.
        Ok(Some((file.u32_at(at)? as i32, file.u32_at(at + 4)? as i32)))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:73-112
    // ```c++
    //     x_Open(algo_id);
    //
    //     int page, vol, off;
    //
    //     if (s_BinarySearch(m_GiIndex, m_NumIndex, GI_TO(Uint4, gi), page)) {
    //         vol = m_GiIndex[m_NumIndex + page * 2];
    //         off = m_GiIndex[m_NumIndex + page * 2 + 1];
    //     } else {
    //         if (page == -1) return; // no mask;
    //
    //         int pagesize(m_PageSize);
    //         if (page * m_PageSize + pagesize > m_NumGi) {
    //             pagesize = m_NumGi - page * m_PageSize;
    //         }
    //
    //         TIndx begin = page * m_PageSize * (m_GiSize + m_OffsetSize);
    //         TIndx end = begin + pagesize * (m_GiSize + m_OffsetSize);
    //         const Uint4 * offset = (const Uint4 *)
    //               m_OffsetFile.GetFileDataPtr(m_OffsetLease, begin, end);
    //
    //         if (!s_BinarySearch(offset, pagesize, GI_TO(Uint4, gi), page))  return;
    //
    //         vol = offset[pagesize + page * 2];
    //         off = offset[pagesize + page * 2 + 1];
    //     }
    //
    //     _ASSERT(vol >= 0);
    //     _ASSERT(vol < m_NumVols);
    //
    //     // Retrieving the mask data
    //     const Int4 * datap = (const Int4 *)
    //               m_DataFile[vol]->GetFileDataPtr(*m_DataLease[vol], off, off+4);
    //     Int4 n = *datap;
    //
    //     // Remapping the mask data
    //     datap = (const Int4 *)
    //             m_DataFile[vol]->GetFileDataPtr(*m_DataLease[vol], off+4, off + 8*n + 4);
    //
    //     ranges.append(datap, n);
    //     return;
    // }
    // ```
    // The arrays are indexed as 4-byte GIs and 4-byte (vol, off) pairs whatever the header's
    // sizes say; each read is checked against the file (NCBI reads its mapping). A `vol` outside
    // the data files (NCBI's debug-only `_ASSERT`, then an index past its vector) is a format
    // error in LOSAT's words; the sums are 64-bit, where NCBI's `int` ones could overflow.
    /// `GetMaskData` of the key `key` (the sequence's GI as a `Uint4`, `GiMasks::seq_gi_key`):
    /// its ranges appended to `out` in the order stored, none without a record.
    pub(super) fn ranges_into(
        &self,
        key: u32,
        out: &mut Vec<(u32, u32)>,
    ) -> Result<(), SeqDbError> {
        let Some((vol, off)) = self.locate(key)? else {
            return Ok(());
        };
        let Some(data) = usize::try_from(vol).ok().and_then(|vol| self.data.get(vol)) else {
            return Err(SeqDbError::Format {
                path: self.index.path().to_path_buf(),
                message: format!(
                    "the GI-mask record of GI {key} is in data volume {vol}, of {} volumes",
                    self.data.len()
                ),
            });
        };
        let off = i64::from(off);
        data.check_range(off, off + 4, "the count of a GI-mask record")?;
        let mut count = [0u8; 4];
        data.read_at(off, &mut count)?;
        let n = i64::from(i32::from_ne_bytes(count));
        let end = off + 8 * n + 4;
        data.check_range(off + 4, end, "the ranges of a GI-mask record")?;
        let mut bytes = vec![0u8; (end - off - 4) as usize];
        data.read_at(off + 4, &mut bytes)?;
        for pair in bytes.chunks_exact(8) {
            out.push((
                u32::from_ne_bytes([pair[0], pair[1], pair[2], pair[3]]),
                u32::from_ne_bytes([pair[4], pair[5], pair[6], pair[7]]),
            ));
        }
        Ok(())
    }
}

/// The `.nog` file of a volume (`CSeqDBGiIndex`): the GI of each OID.
#[derive(Debug)]
struct GiIndex {
    file: MappedFile,
    size: i32,
    num_oids: i32,
}

/// The GI masks of a database whose alias file names them (`CSeqDBGiMask`): one algorithm per
/// name of `MASKLIST`, with ids 0 to N-1 in its order.
#[derive(Debug)]
pub struct GiMasks {
    names: Vec<String>,
    search: SearchPath,
    seqtype: SeqType,
    /// The files of each algorithm, opened on first use (a failed opening is tried again, as
    /// NCBI's `x_Open` is).
    opened: Box<[OnceLock<Arc<GiMaskFiles>>]>,
    /// The `.nog` file of each volume of the database, opened on first use.
    gi_indexes: Box<[OnceLock<Option<GiIndex>>]>,
}

impl GiMasks {
    /// The masks `names` (`MASKLIST`) of a database of `volumes` volumes of `seqtype`, whose
    /// files are looked for along `search`. Nothing is read here.
    pub(super) fn new(
        names: Vec<String>,
        search: SearchPath,
        seqtype: SeqType,
        volumes: usize,
    ) -> Self {
        let opened = names.iter().map(|_| OnceLock::new()).collect();
        let gi_indexes = (0..volumes).map(|_| OnceLock::new()).collect();
        Self {
            names,
            search,
            seqtype,
            opened,
            gi_indexes,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:99-105
    // ```c++
    //     void GetAvailableMaskAlgorithms(vector <int> & algo) const {
    //         algo.clear();
    //         for (unsigned int i=0; i<m_MaskNames.size(); ++i) {
    //             algo.push_back(i);
    //         }
    //         return;
    //     }
    // ```
    /// `GetAvailableMaskAlgorithms`: 0 to N-1. No file is read.
    pub fn ids(&self) -> Vec<i32> {
        (0..self.names.len() as i32).collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:130-143
    // ```c++
    //     string GetAvailableAlgorithmNames() const {
    //         CNcbiOstrstream retval;
    //         retval << endl
    //                << "Available filtering algorithm(s):"
    //                << endl << endl;
    //         retval << setw(14) << left << "Algorithm ID"
    //                << setw(40) << left << "Algorithm name" << endl;
    //         for (unsigned int id=0; id < m_MaskNames.size(); ++id) {
    //             retval << "    " << setw(10) << left << id
    //                    << setw(40) << left << m_MaskNames[id] << endl;
    //         }
    //         return CNcbiOstrstreamToString(retval);
    //     }
    // ```
    fn available_names(&self) -> String {
        let mut out = String::from("\nAvailable filtering algorithm(s):\n\n");
        push_padded(&mut out, "Algorithm ID", 14);
        push_padded(&mut out, "Algorithm name", 40);
        out.push('\n');
        for (id, name) in self.names.iter().enumerate() {
            out.push_str("    ");
            push_padded(&mut out, &id.to_string(), 10);
            push_padded(&mut out, name, 40);
            out.push('\n');
        }
        out
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:110-120
    // ```c++
    //     int GetAlgorithmId(const string & algo_name) const {
    //         for (unsigned int i=0; i<m_MaskNames.size(); ++i) {
    //             if (m_MaskNames[i] == algo_name) return i;
    //         }
    //         CNcbiOstrstream oss;
    //         oss << "Filtering algorithm " << algo_name
    //             << " does not exist." << endl;
    //         oss << GetAvailableAlgorithmNames();
    //         NCBI_THROW(CSeqDBException, eArgErr,
    //                        CNcbiOstrstreamToString(oss));
    //     }
    // ```
    /// `GetAlgorithmId`: the id of the mask `name`, or the text of NCBI's exception.
    pub fn algorithm_id(&self, name: &str) -> Result<i32, String> {
        match self.names.iter().position(|known| known == name) {
            Some(id) => Ok(id as i32),
            None => Err(format!(
                "Filtering algorithm {name} does not exist.\n{}",
                self.available_names()
            )),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.hpp:169-178
    // ```c++
    //     void x_VerifyAlgorithmId(int algo_id) const {
    //         if (algo_id < 0 || algo_id >= (int)m_MaskNames.size()) {
    //             CNcbiOstrstream oss;
    //             oss << "Filtering algorithm ID " << algo_id
    //                 << " does not exist." << endl;
    //             oss << GetAvailableAlgorithmNames();
    //             NCBI_THROW(CSeqDBException, eArgErr,
    //                        CNcbiOstrstreamToString(oss));
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:115-122
    // ```c++
    // void CSeqDBGiMask::x_Open(Int4              algo_id)
    //
    // {
    //     if (algo_id == m_AlgoId) {
    //         return;
    //     }
    //
    //     x_VerifyAlgorithmId(algo_id);
    // ```
    /// `x_Open(id)`: the files of the algorithm `id`. An id outside 0 to N-1 is NCBI's text in
    /// [`SeqDbError::Range`] (a search reports it as `BLAST Database error`).
    pub(super) fn open(&self, id: i32) -> Result<Arc<GiMaskFiles>, SeqDbError> {
        let Some(at) = usize::try_from(id).ok().filter(|&at| at < self.names.len()) else {
            return Err(SeqDbError::Range(format!(
                "Filtering algorithm ID {id} does not exist.\n{}",
                self.available_names()
            )));
        };
        if let Some(files) = self.opened[at].get() {
            return Ok(Arc::clone(files));
        }
        let files = Arc::new(GiMaskFiles::open(&self.names[at], &self.search)?);
        Ok(Arc::clone(self.opened[at].get_or_init(|| files)))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgimask.cpp:56-63
    // ```c++
    // const string &
    // CSeqDBGiMask::GetDesc(int               algo_id,
    //                       CSeqDBLockHold  & locked)
    // {
    //     //m_Atlas.Lock(locked);
    //     x_Open(algo_id);
    //     return m_Desc;
    // }
    // ```
    /// `GetDesc(id)`: the stored description of the algorithm `id` and the index file that holds
    /// it. A description that is not UTF-8 is LOSAT's refusal.
    pub(super) fn description(&self, id: i32) -> Result<(String, PathBuf), SeqDbError> {
        let files = self.open(id)?;
        let path = files.index.path().to_path_buf();
        match String::from_utf8(files.desc.clone()) {
            Ok(desc) => Ok((desc, path)),
            Err(_) => Err(SeqDbError::Unsupported {
                path,
                message: "the description of the GI mask is not UTF-8 text".to_string(),
            }),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:257-268
    // ```c++
    //     if (!m_OidFileOpened &&
    //          CSeqDBGiIndex::IndexExists(m_VolName, (m_IsAA?'p':'n')) &&
    //          m_Idx->GetNumOIDs() != 0) {
    //         m_GiIndex =
    //             new CSeqDBGiIndex(m_Atlas,
    //                               m_VolName,
    //                               (m_IsAA?'p':'n'));
    //     }
    //     m_OidFileOpened = true;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:57-64
    // ```c++
    //     if (m_NumOIDs == 0) {
    //
    //         data = m_Lease.GetFileDataPtr(0);
    //
    //         // TODO we may want to check the version number and file type
    //         m_Size = (Int4) SeqDB_GetStdOrd((Int4 *) (data+8));
    //         m_NumOIDs = (Int4) SeqDB_GetStdOrd((Int4 *) (data+12));
    //     }
    // ```
    /// The `.nog` file of the volume `at` (`volume`), opened on first use; `None` where NCBI has
    /// none (no such file, or a volume without sequences).
    fn gi_index(&self, at: usize, volume: &Volume) -> Result<Option<&GiIndex>, SeqDbError> {
        if let Some(index) = self.gi_indexes[at].get() {
            return Ok(index.as_ref());
        }
        let path = with_extension(volume.base(), self.seqtype, "og");
        let index = if super::alias::exists(&path) && volume.num_oids() != 0 {
            let file = MappedFile::open(&path)?;
            let mut header = [0u8; 8];
            file.read_at(8, &mut header)?;
            Some(GiIndex {
                file,
                size: i32::from_be_bytes([header[0], header[1], header[2], header[3]]),
                num_oids: i32::from_be_bytes([header[4], header[5], header[6], header[7]]),
            })
        } else {
            None
        };
        Ok(self.gi_indexes[at].get_or_init(|| index).as_ref())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:901-915
    // ```c++
    //     int vol_oid = 0;
    //     if (const CSeqDBVol * vol = m_VolSet.FindVol(oid, vol_oid)) {
    //         // Try lookup *.nxg first
    //         TGi gi = vol->GetSeqGI(vol_oid, locked);
    //         if (gi>=ZERO_GI) return gi;
    //         // Fall back to parsing deflines
    //         list< CRef<CSeq_id> > ids =
    //             vol->GetSeqIDs(vol_oid);
    //         ITERATE(list< CRef<CSeq_id> >, id, ids) {
    //             if ((**id).IsGi()) {
    //                 return (**id).GetGi();
    //             }
    //         }
    //         // No GI found
    //         return INVALID_GI;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:66-70
    // ```c++
    //     if (oid >= m_NumOIDs || oid < 0) return INVALID_GI;
    //
    //     TIndx offset = oid * m_Size + 32;
    //     data = m_Lease.GetFileDataPtr(offset);
    //     return GI_FROM(Uint4, SeqDB_GetStdOrd((Uint4 *) data));
    // ```
    // A `.nog` value comes from a `Uint4`, so it is never below 0 and is returned as it is (the
    // -1 that makeblastdb writes for a sequence without a GI reads as 4294967295). The mask
    // lookup takes `GI_TO(Uint4, gi)`, so `INVALID_GI` is 4294967295 too (survey D3).
    /// The key of the sequence `vol_oid` of the volume `at` (`volume`) in the GI masks:
    /// `(Uint4) x_GetSeqGI(oid)`.
    pub(super) fn seq_gi_key(
        &self,
        at: usize,
        volume: &Volume,
        vol_oid: u32,
    ) -> Result<u32, SeqDbError> {
        if let Some(index) = self.gi_index(at, volume)? {
            if i64::from(vol_oid) < i64::from(index.num_oids) {
                let offset = i64::from(vol_oid) * i64::from(index.size) + 32;
                let mut bytes = [0u8; 4];
                index.file.read_at(offset, &mut bytes)?;
                return Ok(u32::from_be_bytes(bytes));
            }
        }
        let set = volume.defline_set(vol_oid)?;
        let gi = set.ids().find_map(|id| id.gi());
        Ok(gi.map_or(u32::MAX, |gi| gi as u32))
    }
}

#[cfg(test)]
#[path = "gi_mask_tests.rs"]
mod tests;
