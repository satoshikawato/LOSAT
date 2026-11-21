//! The files of a column of a volume (`CWriteDB_Column`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_column.cpp`): the index file
//! (`.?aa`), the data file (`.?ab`) and the second data file (`.?ac`). Every column of a
//! volume has both byte orders; only the mask column writes its second blobs, so the columns
//! that `makeblastdb -input_type blastdb` copies (`CreateUserColumn`) never create their `.?ac`.

use std::collections::BTreeMap;

use crate::writedb::files::WriteDbFile;
use crate::writedb::WriteDbError;

/// `TColumnMeta`: the metadata of a column, by key in the order of
/// `std::map<string, string>`.
pub(crate) type ColumnMeta = BTreeMap<Vec<u8>, Vec<u8>>;

/// The two blobs of one column for one sequence (`m_Blobs[2 * col_id]` and
/// `m_Blobs[2 * col_id + 1]`): the mask column's offsets in big-endian and in little-endian byte
/// order, or a copied column's blob and its second one (`SetBlobData` called twice), which the
/// volume checks against the column's size limit and does not write. Empty for a sequence
/// without data in the column.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BlobPair<'a> {
    pub(crate) blob: &'a [u8],
    pub(crate) blob2: &'a [u8],
}

impl BlobPair<'_> {
    /// No data for the sequence.
    pub(crate) const EMPTY: BlobPair<'static> = BlobPair {
        blob: &[],
        blob2: &[],
    };
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1493-1535
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:547-552
// NCBI makes a column (the mask column too) in the open volume as well, with an empty blob for
// each sequence already written there, and adds metadata to the open volume's column.
// makeblastdb makes its columns and their metadata before the first sequence.
/// LOSAT's refusal of a column, of the mask column or of a column's metadata made after the
/// first volume was created.
pub(crate) fn unsupported_late_column() -> WriteDbError {
    WriteDbError::Unsupported(
        "the columns of a database (the mask column too) and their metadata are made before \
         its first volume is written; LOSAT does not add a column, a masking algorithm or \
         column metadata to an open volume"
            .to_string(),
    )
}

/// The creation dates of the columns' index files.
#[derive(Debug, Default)]
pub(crate) struct ColumnDates {
    /// The dates by volume, set by a test (`set`); a volume without one reads the clock.
    dates: BTreeMap<usize, Vec<u8>>,
    /// The date when the clock cannot be read.
    fallback: Option<Vec<u8>>,
}

impl ColumnDates {
    /// Test hook: the creation date of the columns of volume `index`.
    #[cfg(test)]
    pub(crate) fn set(&mut self, index: usize, date: &[u8]) {
        self.dates.insert(index, date.to_vec());
    }

    /// The date of the build, for a clock that cannot be read.
    pub(crate) fn set_fallback(&mut self, date: &[u8]) {
        self.fallback = Some(date.to_vec());
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:199-199
    /// The date stored in the index file of a column of volume `index`: the one a test set, or
    /// the current local time in NCBI's `CTime` text form.
    pub(crate) fn date(&self, index: usize) -> Result<Vec<u8>, WriteDbError> {
        if let Some(date) = self.dates.get(&index) {
            return Ok(date.clone());
        }
        self.now()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_gimask.cpp:197-197
    /// The current local time in NCBI's `CTime` text form, the date of the build when the clock
    /// cannot be read: the date of a column made now, and of the index files of a GI mask
    /// (`super::gi_mask`), which no test sets.
    pub(crate) fn now(&self) -> Result<Vec<u8>, WriteDbError> {
        if let Some(now) = super::time::LocalTime::now() {
            return Ok(now.log_time().into_bytes());
        }
        self.fallback.clone().ok_or_else(|| {
            WriteDbError::Other(
                "the system has no clock for the dates of the database files".to_string(),
            )
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:253-280
// Only non-negative values are written here.
fn write_var_int(out: &mut Vec<u8>, x: u64) {
    let mut buf = [0u8; 16];
    let mut ptr = buf.len();
    let mut ux = x;
    ptr -= 1;
    buf[ptr] = (ux & 0x3F) as u8;
    ux >>= 6;
    while ux != 0 {
        ptr -= 1;
        buf[ptr] = (ux & 0x7F) as u8 | 0x80;
        ux >>= 7;
    }
    out.extend_from_slice(&buf[ptr..]);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:393-403
/// `WriteString(str, eSizeVar)`.
pub(crate) fn write_string_var(out: &mut Vec<u8>, text: &[u8]) {
    write_var_int(out, text.len() as u64);
    out.extend_from_slice(text);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:562-585
/// `WritePadBytes` of the string kind to a multiple of `align`: 1 to `align` bytes, `#` and a
/// NUL at the end.
pub(crate) fn write_pad_bytes_string(out: &mut Vec<u8>, align: usize) {
    let pads = align - out.len() % align;
    out.extend(std::iter::repeat_n(b'#', pads - 1));
    out.push(0);
}

/// `CWriteDB_ColumnData`: a data file and the length of what it holds.
#[derive(Debug)]
struct ColumnData {
    file: WriteDbFile,
    /// `m_DataLength`.
    data_length: u64,
}

impl ColumnData {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:342-353
    /// `WriteBlob`: the length of the file after the blob.
    fn write_blob(&mut self, blob: &[u8]) -> Result<u64, WriteDbError> {
        if blob.is_empty() {
            return Ok(self.data_length);
        }
        if !self.file.created() {
            self.file.create_now();
        }
        self.data_length = u64::from(self.file.write(blob)?);
        Ok(self.data_length)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:362-365
    fn can_fit(&self, size: usize) -> bool {
        self.data_length.wrapping_add(size as u64) < self.file.max_file_size()
    }
}

/// `CWriteDB_ColumnIndex`: the index file, kept in memory until the column is
/// closed.
#[derive(Debug)]
struct ColumnIndex {
    file: WriteDbFile,
    title: Vec<u8>,
    date: Vec<u8>,
    meta: ColumnMeta,
    /// `m_Offsets` without its first (zero) element: the end of each blob in
    /// the data file, as the 4-byte integers of the file.
    offsets: Vec<u32>,
    /// `m_OIDs`.
    oids: u32,
    /// `m_DataLength`: set when the first offset is written.
    data_length: u64,
}

impl ColumnIndex {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:231-288
    // NCBI builds the header when the first offset is written and again when
    // the file is flushed (for the final count and data length); the second
    // build overwrites the first byte for byte except in those two fields,
    // so the header is built once here, with the final values.
    fn build_header(&self, data_length: u64) -> Vec<u8> {
        let mut header = Vec::with_capacity(96);
        for value in [1u32, 1, 4, self.oids] {
            header.extend_from_slice(&value.to_be_bytes());
        }
        header.extend_from_slice(&data_length.to_be_bytes());
        let meta_data_p = header.len();
        header.extend_from_slice(&0u32.to_be_bytes());
        let array_offset_p = header.len();
        header.extend_from_slice(&0u32.to_be_bytes());
        write_string_var(&mut header, &self.title);
        write_string_var(&mut header, &self.date);
        let meta_off = header.len() as u32;
        header[meta_data_p..meta_data_p + 4].copy_from_slice(&meta_off.to_be_bytes());
        write_var_int(&mut header, self.meta.len() as u64);
        for (key, value) in &self.meta {
            write_string_var(&mut header, key);
            write_string_var(&mut header, value);
        }
        write_pad_bytes_string(&mut header, 8);
        let array_off = header.len() as u32;
        header[array_offset_p..array_offset_p + 4].copy_from_slice(&array_off.to_be_bytes());
        header
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:206-229
    fn write_blob_index(&mut self, offset: u64) {
        if self.oids == 0 {
            self.data_length = self.build_header(0).len() as u64 + 4;
        }
        self.offsets.push(offset as u32);
        self.oids += 1;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:313-316
    fn can_fit(&self) -> bool {
        const ENTRY_SIZE: u64 = 4;
        self.data_length + ENTRY_SIZE < self.file.max_file_size()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:290-311
    /// `x_Flush` followed by `CWriteDB_File::Close`; `data_length` is that
    /// of the data file.
    fn close(&mut self, data_length: u64) -> Result<(), WriteDbError> {
        if data_length != 0 {
            if !self.file.created() {
                self.file.create_now();
            }
            let mut out = self.build_header(data_length);
            out.reserve(4 * (self.offsets.len() + 1));
            out.extend_from_slice(&0u32.to_be_bytes());
            for offset in &self.offsets {
                out.extend_from_slice(&offset.to_be_bytes());
            }
            self.file.write(&out)?;
        }
        self.offsets = Vec::new();
        self.file.close();
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.hpp:258-261
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1206-1210
// A new volume makes each column with the default of `CreateColumn`, both byte orders: every
// column has `m_DFile2`. A column with one byte order is made only in an open volume, which
// LOSAT refuses (`unsupported_late_column`).
/// `CWriteDB_Column`: the files of a column with both byte orders (`AddByteOrder`).
#[derive(Debug)]
pub(crate) struct Column {
    index: ColumnIndex,
    data: ColumnData,
    /// `m_DFile2`.
    data2: ColumnData,
    closed: bool,
}

/// The letters of the column ids (`CWriteDB_Volume::CreateColumn`): a volume has at most 36
/// columns.
pub(crate) const COLUMN_LETTERS: &[u8; 36] = b"abcdefghijklmnopqrstuvwxyz0123456789";

impl Column {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:45-78
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:505-543
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:199-199
    // The files are `.<p|n><letter of col_id><a|b|c>`. No file is created here. `date` is the
    // current local time in NCBI's `CTime` text form. The caller has checked `col_id < 36`
    // (`COLUMN_LETTERS`).
    /// The column `col_id` of volume `index`, titled `title`, with both byte orders.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        dbname: &str,
        protein: bool,
        index: u32,
        col_id: usize,
        title: &[u8],
        meta: &ColumnMeta,
        date: &[u8],
        max_file_size: u64,
    ) -> Self {
        let extension = |last: u8| -> String {
            let bytes = [
                if protein { b'p' } else { b'n' },
                COLUMN_LETTERS[col_id],
                last,
            ];
            String::from_utf8_lossy(&bytes).into_owned()
        };
        Self {
            index: ColumnIndex {
                file: WriteDbFile::lazy(dbname, &extension(b'a'), index, max_file_size),
                title: title.to_vec(),
                date: date.to_vec(),
                meta: meta.clone(),
                offsets: Vec::new(),
                oids: 0,
                data_length: 0,
            },
            data: ColumnData {
                file: WriteDbFile::lazy(dbname, &extension(b'b'), index, max_file_size),
                data_length: 0,
            },
            data2: ColumnData {
                file: WriteDbFile::lazy(dbname, &extension(b'c'), index, max_file_size),
                data_length: 0,
            },
            closed: false,
        }
    }

    pub(crate) fn can_fit(&self, size: usize) -> bool {
        self.index.can_fit() && self.data.can_fit(size)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:93-107
    /// `AddBlob(blob, blob2)` (`AddBlob(blob)` is the call with an empty `blob2`); an empty blob
    /// adds an entry without data, and an empty `blob2` writes nothing.
    pub(crate) fn add_blob(&mut self, blob: &[u8], blob2: &[u8]) -> Result<(), WriteDbError> {
        let data_size = self.data.write_blob(blob)?;
        self.index.write_blob_index(data_size);
        self.data2.write_blob(blob2)?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:109-114
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.index.close(self.data.data_length)?;
        self.data.file.close();
        self.data2.file.close();
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:121-133
    pub(crate) fn rename_single(&mut self) {
        self.index.file.rename_single();
        self.data.file.rename_single();
        self.data2.file.rename_single();
    }

    pub(crate) fn rename_file_index(&mut self, num_digits: u32) {
        self.index.file.rename_file_index(num_digits);
        self.data.file.rename_file_index(num_digits);
        self.data2.file.rename_file_index(num_digits);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:84-91
    /// `ListFiles(files, true)`: the files of a column that holds data.
    pub(crate) fn list_files(&self, files: &mut Vec<String>) {
        if self.data.data_length != 0 {
            files.push(self.index.file.filename().to_string());
            files.push(self.data.file.filename().to_string());
            files.push(self.data2.file.filename().to_string());
        }
    }
}
