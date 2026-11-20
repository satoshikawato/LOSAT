//! Column files of a volume (`.?xa` index, `.?xb` / `.?xc` data): per-sequence
//! blobs that `makeblastdb` writes next to the sequence files, among them the mask data
//! of `-mask_data` (the column titled `BlastDb/MaskData`).
//!
//! NCBI reads them with `CSeqDBColumn` (`seqdbcol.cpp`) and `CBlastDbBlob`
//! (`seqdbblob.cpp`); a volume opens all its columns the first time one is asked for
//! (`CSeqDBVol::x_OpenAllColumns`). This module does the same with positioned reads.
//!
//! The integers of the index file are big-endian on every host; the data file holds the
//! blobs as they were written. The two data files of a column differ only inside the
//! blobs of the mask column, whose ranges NCBI writes big-endian to `.?xb` and
//! little-endian to `.?xc`; NCBI writes a user column (`CWriteDB::CreateUserColumn`) to
//! `.?xb` alone.
//!
//! Errors: where NCBI throws a `CSeqDBException` with a fixed text, the text is NCBI's, in
//! [`SeqDbError::Format`]. Where NCBI fails a `SEQDB_FILE_ASSERT` (its text names a source
//! line), the error is [`SeqDbError::Integrity`] in LOSAT's words (a search reports NCBI's
//! text for it); where NCBI reads past its mapping, it is a [`SeqDbError::Format`] in LOSAT's
//! words.

use std::collections::BTreeMap;
use std::path::Path;

use super::file::DbFile;
use super::volume::with_extension;
use super::{SeqDbError, SeqType};

/// The title of the column that holds the mask data (`CSeqDBImpl::x_GetMaskDataColumn`).
pub const MASK_DATA_TITLE: &[u8] = b"BlastDb/MaskData";

/// Size of the fixed fields at the start of a column index file.
const FIXED_FIELD_BYTES: u64 = 32;

/// NCBI's `CSeqDBException` text for a column file that does not parse.
const FORMAT_ERROR: &str = "CSeqDBColumn: File format error.";

/// A read cursor over bytes, as `CBlastDbBlob` reads them. Offsets are `int` in NCBI; here
/// they are `i64`, and a read whose result NCBI leaves undefined is an error.
pub(super) struct Blob<'a> {
    data: &'a [u8],
    offset: i64,
}

impl<'a> Blob<'a> {
    pub(super) fn new(data: &'a [u8], offset: i64) -> Self {
        Self { data, offset }
    }

    pub(super) fn read_offset(&self) -> i64 {
        self.offset
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:547-550
    // ```c++
    // void CBlastDbBlob::SeekRead(int offset)
    // {
    //     m_ReadOffset = offset;
    // }
    // ```
    pub(super) fn seek_read(&mut self, offset: i64) {
        self.offset = offset;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:209-226
    // ```c++
    // const char * CBlastDbBlob::x_ReadRaw(int size, int * offsetp) const
    // {
    //     _ASSERT(offsetp);
    //     _ASSERT(size >= 0);
    //
    //     CTempString s = Str();
    //
    //     int begin = *offsetp;
    //     int end = begin + size;
    //
    //     if (begin > end || end > (int)s.size()) {
    //         NCBI_THROW(CSeqDBException,
    //                    eFileErr,
    //                    "CBlastDbBlob::x_ReadRaw: hit end of data");
    //     }
    //
    //     *offsetp = end;
    //     return s.data() + begin;
    // ```
    /// `x_ReadRaw`. A read that starts before the data, which NCBI would perform on the
    /// bytes in front of its buffer, is LOSAT's own error.
    pub(super) fn read_raw(&mut self, size: i64) -> Result<&'a [u8], String> {
        let begin = self.offset;
        if begin < 0 {
            return Err(format!(
                "a read at offset {begin} of the column data starts before the data"
            ));
        }
        let end = begin.saturating_add(size);
        if begin > end || end > self.data.len() as i64 {
            return Err("CBlastDbBlob::x_ReadRaw: hit end of data".to_string());
        }
        self.offset = end;
        Ok(&self.data[begin as usize..end as usize])
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbblob.hpp:386-423
    // ```c++
    //     /// @param offsetp The offset to read at.
    //     /// @return The value that was read.
    //     template<typename TValue, int TBytes>
    //     TValue x_ReadIntFixed(int * offsetp) const
    //     {
    //         // Check that the value fits in the specified range.
    //         unsigned char * buf = (unsigned char*) x_ReadRaw(TBytes, offsetp);
    //         unsigned char * vbuf = buf - 8 + TBytes;
    //
    //         TValue x = vbuf[7];
    //
    //         switch(TBytes) {
    //         case 8: x |= Uint8(vbuf[0]) << 56;
    //         case 7: x |= Uint8(vbuf[1]) << 48;
    //         case 6: x |= Uint8(vbuf[2]) << 40;
    //         case 5: x |= Uint8(vbuf[3]) << 32;
    //         case 4: x |= Uint4(vbuf[4]) << 24;
    //         case 3: x |= Uint4(vbuf[5]) << 16;
    //         case 2: x |= Uint4(vbuf[6]) << 8;
    //         case 1:
    //             break;
    //         default:
    //             _ASSERT(0);
    //         }
    //
    //         if (TBytes < sizeof(TValue)) {
    //             // This only applies to 'shortened' formats, such as an
    //             // Int8 packed in 5 bytes or an Int4 packed in 3 bytes.
    //             // It only affects signed values, its purposes is to fix
    //             // the numeric sign.  The current design does not use
    //             // shortened formats anywhere.
    //
    //             int bits = (sizeof(TValue)-TBytes) * 8;
    //             x = (x << bits) >> bits;
    //         }
    //
    //         return x;
    //     }
    // ```
    /// `ReadInt4`: four bytes, most significant first.
    pub(super) fn read_i32(&mut self) -> Result<i32, String> {
        let bytes = self.read_raw(4)?;
        Ok(i32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    /// `ReadInt8`: eight bytes, most significant first.
    pub(super) fn read_i64(&mut self) -> Result<i64, String> {
        let bytes = self.read_raw(8)?;
        let mut array = [0u8; 8];
        array.copy_from_slice(bytes);
        Ok(i64::from_be_bytes(array))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:92-115
    // ```c++
    // Int8 CBlastDbBlob::x_ReadVarInt(int * offsetp) const
    // {
    //     CTempString all = Str();
    //     Int8 rv(0);
    //
    //     for(size_t i = *offsetp; i < all.size(); i++) {
    //         int ch = all[i];
    //
    //         if (ch & 0x80) {
    //             // middle
    //             rv = (rv << 7) | (ch & 0x7F);
    //         } else {
    //             // end
    //             rv = (rv << 6) | (ch & 0x3F);
    //             *offsetp = static_cast<int>(i+1);
    //
    //             return (ch & 0x40) ? -rv : rv;
    //         }
    //     }
    //
    //     NCBI_THROW(CSeqDBException,
    //                eFileErr,
    //                "CBlastDbBlob::ReadVarInt: eof while reading integer.");
    // }
    // ```
    /// `ReadVarInt`. The shifts wrap where NCBI's `Int8` arithmetic overflows.
    pub(super) fn read_varint(&mut self) -> Result<i64, String> {
        let mut value: i64 = 0;
        let mut at = usize::try_from(self.offset).unwrap_or(usize::MAX);
        while at < self.data.len() {
            let byte = self.data[at];
            if byte & 0x80 != 0 {
                value = (value << 7) | i64::from(byte & 0x7F);
            } else {
                value = (value << 6) | i64::from(byte & 0x3F);
                self.offset = at as i64 + 1;
                return Ok(if byte & 0x40 != 0 {
                    value.wrapping_neg()
                } else {
                    value
                });
            }
            at += 1;
        }
        Err("CBlastDbBlob::ReadVarInt: eof while reading integer.".to_string())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:169-207
    // ```c++
    // CTempString CBlastDbBlob::x_ReadString(EStringFormat fmt, int * offsetp) const
    // {
    //     int sz = 0;
    //
    //     if (fmt == eSize4) {
    //         sz = x_ReadIntFixed<int,4>(offsetp);
    //     } else if (fmt == eSizeVar) {
    //         sz = static_cast<int>( x_ReadVarInt(offsetp));
    //     }
    //
    //     const char * datap = "";
    //
    //     if (fmt == eNUL) {
    //         CTempString ts = Str();
    //         int zoffset = -1;
    //
    //         for(size_t i = *offsetp; i < ts.size(); i++) {
    //             if (ts[i] == (char)0) {
    //                 zoffset = static_cast<int>(i);
    //                 break;
    //             }
    //         }
    //
    //         if (zoffset == -1) {
    //             NCBI_THROW(CSeqDBException,
    //                        eFileErr,
    //                        "CBlastDbBlob::ReadString: Unterminated string.");
    //         }
    //
    //         datap = ts.data() + *offsetp;
    //         sz = zoffset - *offsetp;
    //         *offsetp = zoffset+1;
    //     } else {
    //         datap = x_ReadRaw(sz, offsetp);
    //     }
    //
    //     return CTempString(datap, sz);
    // }
    // #endif
    // ```
    /// `ReadString(eSizeVar)`: a variable-length size, then that many bytes. The size is
    /// truncated to an `int` as NCBI does.
    pub(super) fn read_string_var(&mut self) -> Result<&'a [u8], String> {
        let size = self.read_varint()? as i32;
        self.read_raw(i64::from(size))
    }

    /// `ReadString(eNUL)`: the bytes up to the next NUL, which is consumed.
    pub(super) fn read_string_nul(&mut self) -> Result<&'a [u8], String> {
        let from = usize::try_from(self.offset).unwrap_or(usize::MAX);
        let nul = self
            .data
            .get(from..)
            .and_then(|rest| rest.iter().position(|&byte| byte == 0));
        match nul {
            Some(at) => {
                self.offset = (from + at + 1) as i64;
                Ok(&self.data[from..from + at])
            }
            None => Err("CBlastDbBlob::ReadString: Unterminated string.".to_string()),
        }
    }
}

/// One column of a volume: its title, metadata and per-sequence blobs.
#[derive(Debug)]
pub struct Column {
    index: DbFile,
    data: DbFile,
    title: Vec<u8>,
    date: Vec<u8>,
    num_oids: i32,
    data_length: i64,
    offset_array_start: i32,
    /// `std::map<string,string>`: ordered by the bytes of the key.
    metadata: BTreeMap<Vec<u8>, Vec<u8>>,
}

impl Column {
    /// The title of the column (not empty).
    pub fn title(&self) -> &[u8] {
        &self.title
    }

    /// The date the column was written.
    #[cfg(test)]
    pub fn date(&self) -> &[u8] {
        &self.date
    }

    /// The number of sequences the column has a blob for.
    #[cfg(test)]
    pub fn num_oids(&self) -> i32 {
        self.num_oids
    }

    /// The length of the data section, as the index file states it.
    #[cfg(test)]
    pub fn data_length(&self) -> i64 {
        self.data_length
    }

    /// The metadata of the column, ordered by key.
    pub fn metadata(&self) -> &BTreeMap<Vec<u8>, Vec<u8>> {
        &self.metadata
    }

    /// Path of the index file.
    pub fn index_path(&self) -> &Path {
        self.index.path()
    }

    /// Path of the data file.
    pub fn data_path(&self) -> &Path {
        self.data.path()
    }

    /// Whether the numbers of the data file read byte-reversed on this host: a little-endian
    /// host that reads the big-endian `.?xb` file because the column has no `.?xc` (`open_all`).
    pub fn reads_byte_reversed(&self) -> bool {
        cfg!(target_endian = "little")
            && !self
                .data
                .path()
                .extension()
                .and_then(|extension| extension.to_str())
                .is_some_and(|extension| extension.ends_with('c'))
    }

    fn format_error(&self, message: impl Into<String>) -> SeqDbError {
        SeqDbError::Format {
            path: self.index.path().to_path_buf(),
            message: message.into(),
        }
    }

    /// The error of a failed `SEQDB_FILE_ASSERT` of the index file.
    fn integrity_error(&self, message: impl Into<String>) -> SeqDbError {
        SeqDbError::Integrity {
            path: self.index.path().to_path_buf(),
            message: message.into(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcol.cpp:102-145
    // ```c++
    // CSeqDBColumn::CSeqDBColumn(const string   & basename,
    //                            const string   & index_extn,
    //                            const string   & data_extn,
    //                            CSeqDBLockHold * lockedp)
    //     : m_AtlasHolder      (lockedp, true),
    //     ...
    //     try {
    //         CSeqDB_Path fn1(basename + "." + index_extn);
    //         CSeqDB_Path fn2(basename + "." + data_extn);
    //
    //         bool found1 = m_IndexFile.Open(fn1);
    //         bool found2 = m_DataFile.Open(fn2);
    //
    //         if (! (found1 && found2)) {
    //             NCBI_THROW(CSeqDBException, eFileErr,
    //                        "Could not open database column files.");
    //         }
    //
    //         x_ReadFields(*lockedp);
    //         x_ReadMetaData(*lockedp);
    //     }
    // ```
    /// Opens the column whose index file is `index_path` and data file `data_path`.
    pub(super) fn open(index_path: &Path, data_path: &Path) -> Result<Self, SeqDbError> {
        let index = DbFile::open(index_path)?;
        let data = DbFile::open(data_path)?;
        let mut column = Self {
            index,
            data,
            title: Vec::new(),
            date: Vec::new(),
            num_oids: 0,
            data_length: 0,
            offset_array_start: 0,
            metadata: BTreeMap::new(),
        };
        let meta_data_start = column.read_fields()?;
        column.read_metadata(meta_data_start)?;
        Ok(column)
    }

    /// `x_ReadFields`. Returns `m_MetaDataStart`.
    fn read_fields(&mut self) -> Result<i32, SeqDbError> {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcol.cpp:206-271
        // ```c++
        // void CSeqDBColumn::x_ReadFields(CSeqDBLockHold & locked)
        // {
        //     const int kFixedFieldBytes = 32;
        //
        //     m_Atlas.Lock(locked);
        //
        //     // First, get the 32 bytes of fields that we know exist.
        //
        //     CBlastDbBlob header;
        //     x_GetFileRange(0, kFixedFieldBytes, e_Index, false, header, locked);
        //
        //     int fmt_version = header.ReadInt4();
        //
        //     if (fmt_version != 1) {
        //         NCBI_THROW(CSeqDBException,
        //                    eFileErr,
        //                    "Column file uses unknown format_version.");
        //     }
        //
        //     int column_type = header.ReadInt4();
        //
        //     if (column_type != 1) {
        //         NCBI_THROW(CSeqDBException,
        //                    eFileErr,
        //                    "Column file uses unknown data type.");
        //     }
        //
        //     int offset_size = header.ReadInt4();
        //
        //     if (offset_size != 4) {
        //         NCBI_THROW(CSeqDBException,
        //                    eFileErr,
        //                    "Column file uses unsupported offset size.");
        //     }
        //
        //     m_NumOIDs = header.ReadInt4();
        //     m_DataLength = header.ReadInt8();
        //     m_MetaDataStart = header.ReadInt4();
        //     m_OffsetArrayStart = header.ReadInt4();
        //
        //     SEQDB_FILE_ASSERT(m_NumOIDs || (! m_DataLength));
        //
        //     SEQDB_FILE_ASSERT(m_MetaDataStart >= 0);
        //     SEQDB_FILE_ASSERT(m_OffsetArrayStart >= m_MetaDataStart);
        //     SEQDB_FILE_ASSERT(m_IndexFile.GetFileLength() >= m_OffsetArrayStart);
        //
        //     // Now we know how long the header actually is, so expand the blob
        //     // to reference the whole thing.  (The memory lease should already
        //     // hold the data, so this will just adjust a few integer fields.)
        //
        //     x_GetFileRange(0, m_MetaDataStart, e_Index, false, header, locked);
        //
        //     // Get string type header fields.
        //
        //     m_Title = header.ReadString (kStringFmt);
        //     m_Date  = header.ReadString (kStringFmt);
        //
        //     SEQDB_FILE_ASSERT(m_Title.size());
        //     SEQDB_FILE_ASSERT(m_Date.size());
        //
        //     if (header.GetReadOffset() != m_MetaDataStart) {
        //         NCBI_THROW(CSeqDBException,
        //                    eFileErr,
        //                    "CSeqDBColumn: File format error.");
        //     }
        // }
        // ```
        // Reading 32 bytes of a shorter file fails `SEQDB_FILE_ASSERT(m_Length >= end)`.
        if self.index.len() < FIXED_FIELD_BYTES {
            return Err(self.integrity_error(format!(
                "the column index file has {} bytes, fewer than its 32 fixed bytes",
                self.index.len()
            )));
        }
        let header = self.index.read_range(0, FIXED_FIELD_BYTES)?;
        let mut blob = Blob::new(&header, 0);
        let fmt_version = blob.read_i32().map_err(|e| self.format_error(e))?;
        if fmt_version != 1 {
            return Err(self.format_error("Column file uses unknown format_version."));
        }
        let column_type = blob.read_i32().map_err(|e| self.format_error(e))?;
        if column_type != 1 {
            return Err(self.format_error("Column file uses unknown data type."));
        }
        let offset_size = blob.read_i32().map_err(|e| self.format_error(e))?;
        if offset_size != 4 {
            return Err(self.format_error("Column file uses unsupported offset size."));
        }
        self.num_oids = blob.read_i32().map_err(|e| self.format_error(e))?;
        self.data_length = blob.read_i64().map_err(|e| self.format_error(e))?;
        let meta_data_start = blob.read_i32().map_err(|e| self.format_error(e))?;
        self.offset_array_start = blob.read_i32().map_err(|e| self.format_error(e))?;

        if self.num_oids == 0 && self.data_length != 0 {
            return Err(self.integrity_error(format!(
                "the column index file states {} data bytes for no sequences",
                self.data_length
            )));
        }
        if meta_data_start < 0 {
            return Err(self.integrity_error(format!(
                "the column index file states a negative metadata start ({meta_data_start})"
            )));
        }
        if self.offset_array_start < meta_data_start {
            return Err(self.integrity_error(format!(
                "the column index file starts its offset array ({}) before its metadata \
                 ({meta_data_start})",
                self.offset_array_start
            )));
        }
        if self.index.len() < self.offset_array_start as u64 {
            return Err(self.integrity_error(format!(
                "the column index file has {} bytes, fewer than the start of its offset array \
                 ({})",
                self.index.len(),
                self.offset_array_start
            )));
        }
        // Reading the range 0..meta_data_start fails `SEQDB_FILE_ASSERT(start < end)` when it
        // is empty.
        if meta_data_start == 0 {
            return Err(self.integrity_error("the column index file states a metadata start of 0"));
        }
        let head = self.index.read_range(0, meta_data_start as u64)?;
        // The blob keeps its read offset (32) when it is pointed at the longer range.
        let mut blob = Blob::new(&head, FIXED_FIELD_BYTES as i64);
        self.title = blob
            .read_string_var()
            .map_err(|e| self.format_error(e))?
            .to_vec();
        self.date = blob
            .read_string_var()
            .map_err(|e| self.format_error(e))?
            .to_vec();
        if self.title.is_empty() {
            return Err(self.integrity_error("the column index file has an empty title"));
        }
        if self.date.is_empty() {
            return Err(self.integrity_error("the column index file has an empty date"));
        }
        if blob.read_offset() != i64::from(meta_data_start) {
            return Err(self.format_error(FORMAT_ERROR));
        }
        Ok(meta_data_start)
    }

    /// `x_ReadMetaData`.
    fn read_metadata(&mut self, meta_data_start: i32) -> Result<(), SeqDbError> {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcol.cpp:273-320
        // ```c++
        // void CSeqDBColumn::x_ReadMetaData(CSeqDBLockHold & locked)
        // {
        //     m_Atlas.Lock(locked);
        //
        //     int begin = m_MetaDataStart;
        //     int end = m_OffsetArrayStart;
        //
        //     _ASSERT(begin > 0 && end > begin);
        //
        //     CBlastDbBlob metadata;
        //     x_GetFileRange(begin, end, e_Index, false, metadata, locked);
        //
        //     Int8 count8 = metadata.ReadVarInt();
        //
        //     if (count8 >> 31) {
        //         NCBI_THROW(CSeqDBException,
        //                    eFileErr,
        //                    "CSeqDBColumn: File format error.");
        //     }
        //
        //     int count = (int) count8;
        //
        //     for(int j = 0; j < count; j++) {
        //         string key = metadata.ReadString(kStringFmt);
        //         string value = metadata.ReadString(kStringFmt);
        //
        //         if (m_MetaData.find(key) != m_MetaData.end()) {
        //             NCBI_THROW(CSeqDBException,
        //                        eFileErr,
        //                        "CSeqDBColumn: Error; duplicate metadata key.");
        //         }
        //
        //         m_MetaData[key] = value;
        //     }
        //
        //     // Align to an 8 byte multiple; eString means that we can change
        //     // the alignment of this field without losing compatibility.
        //
        //     metadata.SkipPadBytes(8, CBlastDbBlob::eString);
        //
        //     int header_bytes = m_OffsetArrayStart - m_MetaDataStart;
        //
        //     if (metadata.GetReadOffset() != header_bytes) {
        //         NCBI_THROW(CSeqDBException,
        //                    eFileErr,
        //                    "CSeqDBColumn: File format error.");
        //     }
        // }
        // ```
        // Reading the range `begin..end` fails `SEQDB_FILE_ASSERT(start < end)` when it is
        // empty.
        if self.offset_array_start <= meta_data_start {
            return Err(self.integrity_error("the column index file has no metadata"));
        }
        let bytes = self
            .index
            .read_range(meta_data_start as u64, self.offset_array_start as u64)?;
        let mut blob = Blob::new(&bytes, 0);
        let count = blob.read_varint().map_err(|e| self.format_error(e))?;
        if count >> 31 != 0 {
            return Err(self.format_error(FORMAT_ERROR));
        }
        for _ in 0..count {
            let key = blob
                .read_string_var()
                .map_err(|e| self.format_error(e))?
                .to_vec();
            let value = blob
                .read_string_var()
                .map_err(|e| self.format_error(e))?
                .to_vec();
            if self.metadata.contains_key(&key) {
                return Err(self.format_error("CSeqDBColumn: Error; duplicate metadata key."));
            }
            self.metadata.insert(key, value);
        }
        // Align to an 8 byte multiple; eString means that we can change the alignment of
        // this field without losing compatibility.
        blob.read_string_nul().map_err(|e| self.format_error(e))?;
        let header_bytes = i64::from(self.offset_array_start) - i64::from(meta_data_start);
        if blob.read_offset() != header_bytes {
            return Err(self.format_error(FORMAT_ERROR));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcol.cpp:322-352
    // ```c++
    // void CSeqDBColumn::GetBlob(int              oid,
    //                            CBlastDbBlob   & blob,
    //                            bool             keep,
    //                            CSeqDBLockHold * lockedp)
    // {
    //     _ASSERT(0 == blob.Size());
    //
    //     CSeqDBLockHold locked2(m_Atlas);
    //
    //     if (lockedp == NULL) {
    //         lockedp = & locked2;
    //     }
    //
    //     int item_size = 4;
    //     int istart = m_OffsetArrayStart + item_size*oid;
    //     int iend = istart + (2 * item_size);
    //
    //     CBlastDbBlob offsets;
    //     x_GetFileRange(istart, iend, e_Index, false, offsets, *lockedp);
    //
    //     int dstart = offsets.ReadInt4();
    //     int dend = offsets.ReadInt4();
    //
    //     SEQDB_FILE_ASSERT(dend >= dstart);
    //
    //     if (dend > dstart) {
    //         x_GetFileRange(dstart, dend, e_Data, keep, blob, *lockedp);
    //     } else {
    //         _ASSERT(! blob.Size());
    //     }
    // }
    // ```
    /// The blob of the sequence `oid` (volume-relative); empty when the sequence has none
    /// ([`Column::read_blobs`] of the one sequence).
    pub fn blob(&self, oid: u32) -> Result<Vec<u8>, SeqDbError> {
        let mut blobs = ColumnBlobs::default();
        self.read_blobs(oid, oid.saturating_add(1), &mut blobs);
        blobs.take_blob(oid)
    }

    /// The end of the sequences of `first..end` (volume-relative) whose two entries of the offset
    /// array can be read: those of the column whose entries the index file holds (`GetBlob`'s
    /// `SEQDB_FILE_ASSERT(m_Length >= end)` of the index file); `first` when there are none.
    fn entries_end(&self, first: u32, end: u32) -> u32 {
        let in_file = (self.index.len() as i64 - i64::from(self.offset_array_start) - 4) / 4;
        let in_file = u32::try_from(in_file.max(0)).unwrap_or(u32::MAX);
        let num_oids = u32::try_from(self.num_oids.max(0)).unwrap_or(0);
        end.min(num_oids).min(in_file).max(first)
    }

    /// The error of the sequence `oid` at the end of the entries that can be read
    /// ([`Column::entries_end`]).
    fn entries_end_error(&self, oid: u32) -> SeqDbError {
        if oid >= u32::try_from(self.num_oids.max(0)).unwrap_or(0) {
            SeqDbError::Range(format!(
                "sequence {oid} is outside the column ({} sequences)",
                self.num_oids
            ))
        } else {
            self.integrity_error(format!(
                "the offset array of the column ends before the entry of sequence {oid}"
            ))
        }
    }

    /// The offset array's entries of the sequences `first..=last` into `entries` (cleared first),
    /// with one read (`last` below [`Column::entries_end`]'s end, or at it).
    fn read_entries(
        &self,
        first: u32,
        last: u32,
        entries: &mut Vec<i32>,
    ) -> Result<(), SeqDbError> {
        entries.clear();
        let istart = i64::from(self.offset_array_start) + 4 * i64::from(first);
        let mut bytes = vec![0u8; 4 * (last - first + 1) as usize];
        self.index.read_exact_at(&mut bytes, istart as u64)?;
        entries.extend(
            bytes
                .chunks_exact(4)
                .map(|entry| i32::from_be_bytes([entry[0], entry[1], entry[2], entry[3]])),
        );
        Ok(())
    }

    /// `GetBlob`'s checks of the blob of the sequence `oid`, whose data range is `dstart..dend`,
    /// in NCBI's order.
    fn check_blob(&self, oid: u32, dstart: i32, dend: i32) -> Result<(), SeqDbError> {
        if dend < dstart {
            Err(self.integrity_error(format!(
                "the data range of sequence {oid} ends ({dend}) before it starts ({dstart})"
            )))
        } else if dend == dstart {
            Ok(())
        } else if dstart < 0 {
            // A negative start is undefined in NCBI (a pointer before the mapping).
            Err(self.format_error(format!(
                "the data range of sequence {oid} starts at {dstart}"
            )))
        } else {
            // A range past the end of the data file fails `SEQDB_FILE_ASSERT(m_Length >= end)`:
            // the checks of `DbFile::read_range` before its buffer is sized, as an integrity error.
            self.data
                .x_check_range(dstart as u64, dend as u64)
                .map_err(|error| match error {
                    SeqDbError::Format { path, message } => SeqDbError::Integrity { path, message },
                    error => error,
                })
        }
    }

    /// Whether a sequence of `first..end` (volume-relative) has a blob, by the offset array
    /// alone, read a window at a time: `Some(true)` at the first that has one whose data range
    /// reads, `Some(false)` when a sequence before it cannot be read ([`Column::read_blobs`]'s
    /// checks: a reader of the blobs in order meets that error first), `None` when every one has
    /// an empty blob.
    pub fn first_blob_in(&self, first: u32, end: u32) -> Option<bool> {
        const WINDOW: u32 = 1 << 14;
        let mut entries = Vec::new();
        let mut window = first;
        while window < end {
            let window_end = end.min(window.saturating_add(WINDOW));
            let read_end = self.entries_end(window, window_end);
            if read_end > window {
                if self.read_entries(window, read_end, &mut entries).is_err() {
                    return Some(false);
                }
                for (oid, pair) in (window..).zip(entries.windows(2)) {
                    match self.check_blob(oid, pair[0], pair[1]) {
                        Err(_) => return Some(false),
                        Ok(()) if pair[1] > pair[0] => return Some(true),
                        Ok(()) => {}
                    }
                }
            }
            if read_end < window_end {
                return Some(false);
            }
            window = window_end;
        }
        None
    }

    /// The blobs of the sequences `first..end` (volume-relative, `first < end`) into `blobs`,
    /// whose buffers it reuses: one read of the offset array's entries of the sequences and one
    /// of the data range of their blobs, as `Volume::read_block` reads a block's sequences (NCBI
    /// maps the files and points into them, `GetBlob` per sequence). Each sequence has
    /// [`Column::blob`]'s checks in the order of the sequences; the first that fails keeps its
    /// error for itself and the sequences after it ([`ColumnBlobs::blob`]). A read of the data
    /// range that fails is the error of its first sequence with a blob.
    pub fn read_blobs(&self, first: u32, end: u32, blobs: &mut ColumnBlobs) {
        blobs.first = first;
        blobs.offsets.clear();
        blobs.data.clear();
        blobs.data_start = 0;
        blobs.error = None;
        let read_end = self.entries_end(first, end);
        if read_end < end {
            blobs.error = Some((read_end, self.entries_end_error(read_end)));
        }
        if read_end == first {
            return;
        }
        if let Err(error) = self.read_entries(first, read_end, &mut blobs.offsets) {
            blobs.offsets.clear();
            blobs.error = Some((first, error));
            return;
        }
        // The checks of each sequence in order; the data range of those before the first that
        // fails.
        let mut data: Option<(i32, i32, u32)> = None;
        for oid in first..read_end {
            let i = (oid - first) as usize;
            let (dstart, dend) = (blobs.offsets[i], blobs.offsets[i + 1]);
            if let Err(error) = self.check_blob(oid, dstart, dend) {
                blobs.error = Some((oid, error));
                break;
            }
            if dend > dstart {
                data = Some(match data {
                    None => (dstart, dend, oid),
                    Some((low, high, at)) => (low.min(dstart), high.max(dend), at),
                });
            }
        }
        if let Some((low, high, at)) = data {
            blobs.data.resize((high - low) as usize, 0);
            if let Err(error) = self.data.read_exact_at(&mut blobs.data, low as u64) {
                blobs.data.clear();
                blobs.error = Some((at, error));
                return;
            }
            blobs.data_start = low;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3325-3392
// ```c++
// void CSeqDBVol::x_OpenAllColumns(CSeqDBLockHold & locked)
// {
//     //m_Atlas.Lock(locked);
//
//     if (m_HaveColumns) {
//         return;
//     }
//
//     string alpha("abcdefghijklmnopqrstuvwxyz");
//     string ei("??a"), ed("??b"), ed2("??c");
//
//     ei[0] = ed[0] = ed2[0] = (m_IsAA ? 'p' : 'n');
//
//     map<string,int> unique_titles;
//
//     for(size_t i = 0; i < alpha.size(); i++) {
//         ei[1] = ed[1] = ed2[1] = alpha[i];
//
//         if (CSeqDBColumn::ColumnExists(m_VolName, ei, m_Atlas)) {
//
//             bool big   = CSeqDBColumn::ColumnExists(m_VolName, ed, m_Atlas);
//             bool small = CSeqDBColumn::ColumnExists(m_VolName, ed2, m_Atlas);
//
//             if ( ! (big || small)) continue;
//
//             CRef<CSeqDBColumn> col;
//
//             const Int2 bytetest = 0x0011;
//             const char * ptr = (const char *) &bytetest;
//             if (ptr[0] == 0x11 && small) {
//                 col.Reset(new CSeqDBColumn(m_VolName, ei, ed2, & locked));
//             } else {
//                 col.Reset(new CSeqDBColumn(m_VolName, ei, ed, & locked));
//             }
//
//             string errmsg, errarg;
//
//             string title = col->GetTitle();
//
//             if (unique_titles[title]) {
//                 errmsg = "duplicate column title";
//                 errarg = title;
//             } else {
//                 unique_titles[title] = 1;
//             }
//
//             int noidc(col->GetNumOIDs()), noidv(m_Idx->GetNumOIDs());
//
//             if (noidc != noidv) {
//                 errmsg = "column has wrong #oids";
//                 errarg = NStr::IntToString(noidc) + " vs "
//                     + NStr::IntToString(noidv);
//             }
//
//             if (errmsg.size()) {
//                 if (errarg.size()) {
//                     errmsg += string(" [") + errarg + "].";
//                 }
//                 NCBI_THROW(CSeqDBException, eFileErr,
//                            string("Error: ") + errmsg);
//             }
//
//             m_Columns.push_back(col);
//         }
//     }
//
//     m_HaveColumns = true;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:311-318
// ```c++
//     for(int col_i = 0; col_i < (int)m_Columns.size(); col_i++) {
//         _ASSERT(col_i * 2 < (int) blobs.size());
//         if (col_i == maskcol_id) {
//              m_Columns[col_i]->AddBlob(*blobs[col_i * 2], *blobs[col_i * 2 + 1]);
//         } else {
//              m_Columns[col_i]->AddBlob(*blobs[col_i * 2]);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:93-107
// ```c++
// void CWriteDB_Column::AddBlob(const CBlastDbBlob & blob)
// {
//     // Note that data size is the size *after* the blob has been
//     // written.  The initial (zero) offset is written during file
//     // creation.
//
//     Int8 data_size = m_DFile->WriteBlob(blob);
//     m_IFile->WriteBlobIndex(data_size);
// }
//
// void CWriteDB_Column::AddBlob(const CBlastDbBlob & blob, const CBlastDbBlob & blob2)
// {
//     AddBlob(blob);
//     if (m_UseBothByteOrder) m_DFile2->WriteBlob(blob2);
// }
// ```
// The index file holds the offsets of the `.?xb` blobs (`WriteInt4`, big-endian), and the
// `.?xc` blob of a sequence has the size of its `.?xb` blob, so one index serves either
// data file. Only the mask column writes a `.?xc` (`blob2`, its ranges `WriteInt4_LE`); a
// user column has `.?xb` alone. A blob is the file's bytes either way (`GetBlob`);
// `mask::read_ranges` reads the ranges in host order, as NCBI's copy does, so a column read
// from its `.?xb` alone on a little-endian host gives each number byte-reversed, a
// deterministic result that LOSAT reproduces (line SW, SW8 DECISIONS 40; a refusal remains
// where a consumer would index letters with those ranges).
/// Opens every column of the volume `base` (`CSeqDBVol::x_OpenAllColumns`): for the letters
/// `a` to `z`, the columns whose index file exists, in letter order. `volume_oids` is the
/// number of sequences of the volume.
pub(super) fn open_all(
    base: &Path,
    seqtype: SeqType,
    volume_oids: u32,
) -> Result<Vec<Column>, SeqDbError> {
    let mut columns: Vec<Column> = Vec::new();
    for letter in b'a'..=b'z' {
        let letter = char::from(letter);
        let index_path = with_extension(base, seqtype, &format!("{letter}a"));
        if !super::alias::exists(&index_path) {
            continue;
        }
        let big_path = with_extension(base, seqtype, &format!("{letter}b"));
        let small_path = with_extension(base, seqtype, &format!("{letter}c"));
        let big = super::alias::exists(&big_path);
        let small = super::alias::exists(&small_path);
        if !(big || small) {
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3345-3358
        // ```c++
        //             bool big   = CSeqDBColumn::ColumnExists(m_VolName, ed, m_Atlas);
        //             bool small = CSeqDBColumn::ColumnExists(m_VolName, ed2, m_Atlas);
        //
        //             if ( ! (big || small)) continue;
        //
        //             CRef<CSeqDBColumn> col;
        //
        //             const Int2 bytetest = 0x0011;
        //             const char * ptr = (const char *) &bytetest;
        //             if (ptr[0] == 0x11 && small) {
        //                 col.Reset(new CSeqDBColumn(m_VolName, ei, ed2, & locked));
        //             } else {
        //                 col.Reset(new CSeqDBColumn(m_VolName, ei, ed, & locked));
        //             }
        // ```
        // `ptr[0] == 0x11` is a little-endian host.
        let data_path = if cfg!(target_endian = "little") && small {
            small_path
        } else {
            big_path
        };
        let column = Column::open(&index_path, &data_path)?;

        let mut message = String::new();
        let mut argument = String::new();
        if columns.iter().any(|other| other.title == column.title) {
            message = "duplicate column title".to_string();
            argument = String::from_utf8_lossy(&column.title).into_owned();
        }
        if column.num_oids != volume_oids as i32 {
            message = "column has wrong #oids".to_string();
            argument = format!("{} vs {}", column.num_oids, volume_oids);
        }
        if !message.is_empty() {
            if !argument.is_empty() {
                message = format!("{message} [{argument}].");
            }
            return Err(SeqDbError::Format {
                path: index_path,
                message: format!("Error: {message}"),
            });
        }
        columns.push(column);
    }
    Ok(columns)
}

/// The blobs of a run of sequences of one column ([`Column::read_blobs`]): the offset array's
/// entries of the run and the data range of their blobs, kept to be read again (a search's
/// worker reuses one per block of its preliminary stage). No NCBI counterpart: NCBI points into
/// its mapping of the files.
#[derive(Debug, Default)]
pub struct ColumnBlobs {
    /// The first sequence of the run (volume-relative).
    first: u32,
    /// The offset array's entries from `first` on: sequence `first + i` has the data
    /// `offsets[i]..offsets[i + 1]`.
    offsets: Vec<i32>,
    /// The data file's bytes from `data_start` on: the blobs of the sequences before `error`.
    data: Vec<u8>,
    data_start: i32,
    /// The first sequence whose blob cannot be read, and its error, which the sequences after it
    /// share.
    error: Option<(u32, SeqDbError)>,
}

impl ColumnBlobs {
    /// The blob of the sequence `oid` of the run (empty when it has none), or the error of the
    /// first sequence of the run that cannot be read, at it or before it.
    pub fn blob(&self, oid: u32) -> Result<&[u8], &SeqDbError> {
        if let Some((failed, error)) = &self.error {
            if oid >= *failed {
                return Err(error);
            }
        }
        let i = (oid - self.first) as usize;
        let (dstart, dend) = (self.offsets[i], self.offsets[i + 1]);
        if dend == dstart {
            return Ok(&[]);
        }
        Ok(&self.data[(dstart - self.data_start) as usize..(dend - self.data_start) as usize])
    }

    /// [`ColumnBlobs::blob`] of `oid`, owned.
    fn take_blob(self, oid: u32) -> Result<Vec<u8>, SeqDbError> {
        if let Ok(blob) = self.blob(oid) {
            return Ok(blob.to_vec());
        }
        Err(self.into_error(oid))
    }

    /// The error that [`ColumnBlobs::blob`] gives `oid`, owned (`oid` is at or after the first
    /// sequence that cannot be read).
    pub(super) fn into_error(self, oid: u32) -> SeqDbError {
        match self.error {
            Some((failed, error)) if oid >= failed => error,
            _ => SeqDbError::Range(format!(
                "LOSAT internal error: the blob of sequence {oid} was read"
            )),
        }
    }
}
