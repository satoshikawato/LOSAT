//! Column files of a volume (`.?xa` index, `.?xb` / `.?xc` data): per-sequence
//! blobs that `makeblastdb` writes next to the sequence files, among them the mask data
//! of `-mask_data` (the column titled `BlastDb/MaskData`).
//!
//! NCBI reads them with `CSeqDBColumn` (`c++/src/objtools/blast/seqdb_reader/seqdbcol.cpp`) and `CBlastDbBlob`
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
    pub(super) fn seek_read(&mut self, offset: i64) {
        self.offset = offset;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:209-226
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
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:311-318
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:93-107
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1955-2052
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2082-2115
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcol.cpp:368-375
// NCBI numbers a column of the database by its title (`GetColumnId`: one entry per title,
// holding each volume's own column of that title or none); LOSAT names it by the title. The
// metadata of the database's column is made once from the volumes in their order, and a key
// keeps the value of the first volume that has it. A sequence of a volume without a column of
// the title has an empty blob.
impl super::SeqDb {
    /// `ListColumns`: the titles of the columns of all the volumes, each once, in the byte
    /// order of the titles. Every volume opens its columns (`Volume::columns`), in the order of
    /// the volumes, so the first volume whose columns cannot be opened gives the error.
    pub fn column_titles(&self) -> Result<Vec<Vec<u8>>, SeqDbError> {
        let mut all = std::collections::BTreeSet::new();
        for volume in self.volumes() {
            for column in volume.columns()? {
                all.insert(column.title().to_vec());
            }
        }
        Ok(all.into_iter().collect())
    }

    /// `GetColumnMetaData(GetColumnId(title))`: the metadata of the columns titled `title` of
    /// all the volumes, a key with the value of the first volume (in the order of the volumes)
    /// that has it.
    pub fn column_meta_data(&self, title: &[u8]) -> Result<BTreeMap<Vec<u8>, Vec<u8>>, SeqDbError> {
        let mut merged = BTreeMap::new();
        for volume in self.volumes() {
            let Some(column) = volume
                .columns()?
                .iter()
                .find(|column| column.title() == title)
            else {
                continue;
            };
            for (key, value) in column.metadata() {
                merged.entry(key.clone()).or_insert_with(|| value.clone());
            }
        }
        Ok(merged)
    }

    /// The blobs of the column titled `title` (`GetColumnBlob`) for the sequences of `range`,
    /// read into `block` with one read of the column's offset array and one of its data
    /// ([`Column::read_blobs`]); none when the volume of `range` has no column of the title.
    /// [`ColumnBlock::blob`] gives each sequence's.
    pub fn read_column_block(
        &self,
        title: &[u8],
        range: super::BlockRange,
        block: &mut ColumnBlock,
    ) {
        block.error = None;
        block.present = false;
        let Some(volume) = self.volumes().get(range.volume) else {
            block.error = Some(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes().len()
            )));
            return;
        };
        let columns = match volume.columns() {
            Ok(columns) => columns,
            Err(error) => {
                block.error = Some(error);
                return;
            }
        };
        let Some(column) = columns.iter().find(|column| column.title() == title) else {
            return;
        };
        block.present = true;
        column.read_blobs(range.first, range.end, &mut block.blobs);
    }
}

/// The blobs of one column of the database for the sequences of a block
/// ([`super::SeqDb::read_column_block`]), kept with their buffers from one block to the next.
/// No NCBI counterpart: NCBI points into its mapping of the column for each sequence.
#[derive(Debug, Default)]
pub struct ColumnBlock {
    /// The error of every sequence of the block: the volume's columns could not be opened.
    error: Option<SeqDbError>,
    /// Whether the volume has a column of the title (none: every blob is empty).
    present: bool,
    blobs: ColumnBlobs,
}

impl ColumnBlock {
    /// `GetColumnBlob` of the sequence `vol_oid` (volume-relative) of the block read: its blob
    /// as the data file holds it, empty when it has none. An error is given once, owned: the
    /// caller stops at it.
    pub fn blob(&mut self, vol_oid: u32) -> Result<&[u8], SeqDbError> {
        if let Some(error) = self.error.take() {
            return Err(error);
        }
        if !self.present {
            return Ok(&[]);
        }
        if self.blobs.blob(vol_oid).is_err() {
            return Err(std::mem::take(&mut self.blobs).into_error(vol_oid));
        }
        Ok(self.blobs.blob(vol_oid).unwrap_or(&[]))
    }
}
