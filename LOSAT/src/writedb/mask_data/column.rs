//! The files of a column of a volume (`CWriteDB_Column`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_column.cpp`): the index file
//! (`.?aa`), the data file (`.?ab`) and, for a column with two byte orders,
//! the second data file (`.?ac`). The mask column is the only one LOSAT
//! writes, and it has both byte orders.

use std::collections::BTreeMap;

use crate::writedb::files::WriteDbFile;
use crate::writedb::WriteDbError;

/// `TColumnMeta`: the metadata of a column, by key in the order of
/// `std::map<string, string>`.
pub(crate) type ColumnMeta = BTreeMap<Vec<u8>, Vec<u8>>;

/// The title of the mask column (`x_GetMaskDataColumnId`).
pub(crate) const MASK_COLUMN_TITLE: &[u8] = b"BlastDb/MaskData";

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
/// `WritePadBytes(align, eString)`: 1 to `align` bytes, `#` and a NUL at the end.
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

/// `CWriteDB_Column` with `AddByteOrder`: the three files of a column.
#[derive(Debug)]
pub(crate) struct Column {
    index: ColumnIndex,
    data: ColumnData,
    data2: ColumnData,
    closed: bool,
}

impl Column {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:45-78
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:512-527
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:199-199
    // The column is the first of the volume (`col_id` 0, letter `a`). No file
    // is created here. `date` is `CTime(CTime::eCurrent).AsString()`.
    pub(crate) fn new(
        dbname: &str,
        protein: bool,
        index: u32,
        meta: &ColumnMeta,
        date: &[u8],
        max_file_size: u64,
    ) -> Self {
        let extension = |letter: char| -> &'static str {
            match (protein, letter) {
                (true, 'a') => "paa",
                (true, 'b') => "pab",
                (true, _) => "pac",
                (false, 'a') => "naa",
                (false, 'b') => "nab",
                (false, _) => "nac",
            }
        };
        Self {
            index: ColumnIndex {
                file: WriteDbFile::lazy(dbname, extension('a'), index, max_file_size),
                title: MASK_COLUMN_TITLE.to_vec(),
                date: date.to_vec(),
                meta: meta.clone(),
                offsets: Vec::new(),
                oids: 0,
                data_length: 0,
            },
            data: ColumnData {
                file: WriteDbFile::lazy(dbname, extension('b'), index, max_file_size),
                data_length: 0,
            },
            data2: ColumnData {
                file: WriteDbFile::lazy(dbname, extension('c'), index, max_file_size),
                data_length: 0,
            },
            closed: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:116-119
    pub(crate) fn can_fit(&self, size: usize) -> bool {
        self.index.can_fit() && self.data.can_fit(size)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:93-107
    /// `AddBlob(blob, blob2)`; an empty blob adds an entry without data.
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
