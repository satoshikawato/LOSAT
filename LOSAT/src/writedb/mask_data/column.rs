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
pub type ColumnMeta = BTreeMap<Vec<u8>, Vec<u8>>;

/// The title of the mask column (`x_GetMaskDataColumnId`).
pub const MASK_COLUMN_TITLE: &[u8] = b"BlastDb/MaskData";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:253-280
// ```c++
// int CBlastDbBlob::x_WriteVarInt(Int8 x, int * offsetp)
// {
//     // The variable length integer is written into the end of the 16
//     // byte array shown below.
//     ...
//     char buf[16];
//     int end_ptr = sizeof(buf);
//     int ptr = end_ptr;
//
//     Uint8 ux = (Uint8)((x >= 0) ? x : -x);
//
//     buf[--ptr] = (ux & 0x3F);
//     ux >>= 6;
//     ...
//     while(ux) {
//         buf[--ptr] = (ux & 0x7F) | 0x80;
//         ux >>= 7;
//     }
//
//     int bytes = end_ptr - ptr;
//
//     x_WriteRaw(buf + ptr, bytes, offsetp);
// ```
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
// ```c++
// int CBlastDbBlob::x_WriteString(CTempString str, EStringFormat fmt, int * offsetp)
// {
//     int start_off = offsetp ? *offsetp : m_WriteOffset;
//
//     if (fmt == eSize4) {
//         x_WriteIntFixed<int,4>(static_cast<int>(str.size()), offsetp);
//     } else if (fmt == eSizeVar) {
//         x_WriteVarInt(str.size(), offsetp);
//     }
//
//     x_WriteRaw(str.data(),static_cast<int>( str.size()), offsetp);
// ```
/// `WriteString(str, eSizeVar)`.
fn write_string_var(out: &mut Vec<u8>, text: &[u8]) {
    write_var_int(out, text.len() as u64);
    out.extend_from_slice(text);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbblob.cpp:562-585
// ```c++
// void CBlastDbBlob::WritePadBytes(int align, EPadding fmt)
// {
//     vector<char> pad;
//     CTempString tmp;
//
//     int pads = align ? (m_WriteOffset % align) : 0;
//
//     if (fmt == eSimple) {
//         pads = pads ? (align - pads) : 0;
//     } else {
//         pads = align - pads;
//     }
//
//     if (fmt == eSimple) {
//         for(int i = 0; i < pads; i++) {
//             x_WriteRaw("#", 1, NULL);
//         }
//     } else {
//         for(int i = 1; i < pads; i++) {
//             x_WriteRaw("#", 1, NULL);
//         }
//         char ch = (char)0;
//         x_WriteRaw(& ch, 1, NULL);
//     }
// ```
/// `WritePadBytes(8, eString)`: 1 to 8 bytes, `#` and a NUL at the end.
fn write_pad_bytes_string(out: &mut Vec<u8>, align: usize) {
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
    // ```c++
    // Int8 CWriteDB_ColumnData::WriteBlob(const CBlastDbBlob & blob)
    // {
    //     if (! blob.Size()) {
    //         return m_DataLength;
    //     }
    //
    //     if (! m_Created) {
    //         Create();
    //     }
    //
    //     return m_DataLength = Write(blob.Str());
    // }
    // ```
    /// `WriteBlob`: the length of the file after the blob.
    fn write_blob(&mut self, blob: &[u8]) -> Result<u64, WriteDbError> {
        if blob.is_empty() {
            return Ok(self.data_length);
        }
        self.file.create_file()?;
        self.data_length = u64::from(self.file.write(blob)?);
        Ok(self.data_length)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:362-365
    // ```c++
    // bool CWriteDB_ColumnData::CanFit(int size) const
    // {
    //     return Uint8(m_DataLength + size) < m_MaxFileSize;
    // }
    // ```
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
    // ```c++
    // void CWriteDB_ColumnIndex::x_BuildHeaderFields()
    // {
    //     // The Blob type makes a great binary data stream type.
    //
    //     const int kFormatVersion = 1; // SeqDB has one of these.
    //     const int kColumnType    = 1; // Blob (only choice right now)
    //     const int kOffsetSize    = 4; // Data file offset size (always 4)
    //
    //     m_Header->SeekWrite(0);
    //     m_Header->WriteInt4(kFormatVersion);
    //     m_Header->WriteInt4(kColumnType);
    //     m_Header->WriteInt4(kOffsetSize);
    //     m_Header->WriteInt4(m_OIDs);
    //     m_Header->WriteInt8(m_DataFile->GetDataLength());
    // }
    //
    // void CWriteDB_ColumnIndex::x_BuildHeaderStrings()
    // {
    //     int meta_data_p = m_Header->GetWriteOffset();
    //     m_Header->WriteInt4(0); // metadata start
    //
    //     int array_offset_p = m_Header->GetWriteOffset();
    //     m_Header->WriteInt4(0); // offset array start
    //
    //     m_Header->WriteString(m_Title, kStringFmt);
    //     m_Header->WriteString(m_Date, kStringFmt);
    //
    //     int meta_off = m_Header->GetWriteOffset();
    //     m_Header->WriteInt4(meta_off, meta_data_p);
    //
    //     x_BuildMetaData();
    //
    //     // Align to an 8 byte multiple; eString means that we can change
    //     // the alignment of this field without losing compatibility.
    //     m_Header->WritePadBytes(8, CBlastDbBlob::eString);
    //
    //     int array_off = m_Header->GetWriteOffset();
    //     m_Header->WriteInt4(array_off, array_offset_p);
    // }
    //
    // void CWriteDB_ColumnIndex::x_BuildMetaData()
    // {
    //     m_Header->WriteVarInt(m_MetaData.size());
    //
    //     ITERATE(TColumnMeta, iter, m_MetaData) {
    //         CTempString key = iter->first, value = iter->second;
    //         m_Header->WriteString(key, kStringFmt);
    //         m_Header->WriteString(value, kStringFmt);
    //     }
    // }
    // ```
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
    // ```c++
    // void CWriteDB_ColumnIndex::WriteBlobIndex(Int8 offset)
    // {
    //     _ASSERT(0 == (offset >> 32));
    //
    //     if (m_Header.Empty()) {
    //         m_Header.Reset(new CBlastDbBlob(256));
    //         m_Offsets.Reset(new CBlastDbBlob(4096));
    //
    //         // We build these now so that m_DataLength is accurate.  They
    //         // will be rebuilt just before they are written, when the file
    //         // is closed.
    //
    //         x_BuildHeaderFields();
    //         x_BuildHeaderStrings();
    //
    //         // Offset of first data element (always zero).
    //         m_Offsets->WriteInt4(0);
    //
    //         m_DataLength = m_Header->Size() + m_Offsets->Size();
    //     }
    //
    //     m_Offsets->WriteInt4((Int4) offset);
    //     m_OIDs ++;
    // }
    // ```
    fn write_blob_index(&mut self, offset: u64) {
        if self.oids == 0 {
            self.data_length = self.build_header(0).len() as u64 + 4;
        }
        self.offsets.push(offset as u32);
        self.oids += 1;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:313-316
    // ```c++
    // bool CWriteDB_ColumnIndex::CanFit() const
    // {
    //     return (m_DataLength + kEntrySize) < m_MaxFileSize;
    // }
    // ```
    fn can_fit(&self) -> bool {
        const ENTRY_SIZE: u64 = 4;
        self.data_length + ENTRY_SIZE < self.file.max_file_size()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:290-311
    // ```c++
    // void CWriteDB_ColumnIndex::x_Flush()
    // {
    //     if (! m_DataFile->Empty()) {
    //         if (! m_Created) {
    //             Create();
    //         }
    //
    //         // These need to be rebuilt to write the correct values for
    //         // OID count, total length, and possibly meta data.
    //
    //         x_BuildHeaderFields();
    //         x_BuildHeaderStrings();
    //
    //         Write(m_Header->Str());
    //         Write(m_Offsets->Str());
    //
    //         // We're done with these now, so free up the memory.
    //
    //         m_Header.Reset();
    //         m_Offsets.Reset();
    //     }
    // }
    // ```
    /// `x_Flush` followed by `CWriteDB_File::Close`; `data_length` is that
    /// of the data file.
    fn close(&mut self, data_length: u64) -> Result<(), WriteDbError> {
        if data_length != 0 {
            self.file.create_file()?;
            let mut out = self.build_header(data_length);
            out.reserve(4 * (self.offsets.len() + 1));
            out.extend_from_slice(&0u32.to_be_bytes());
            for offset in &self.offsets {
                out.extend_from_slice(&offset.to_be_bytes());
            }
            self.file.write(&out)?;
        }
        self.offsets = Vec::new();
        self.file.close()
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
    // ```c++
    // CWriteDB_Column::CWriteDB_Column(const string      & dbname,
    //                                  const string      & extn1,
    //                                  const string      & extn2,
    //                                  int                 index,
    //                                  const string      & title,
    //                                  const TColumnMeta & meta,
    //                                  Uint8               max_file_size)
    //                 : m_UseBothByteOrder(false)
    // {
    //     m_DFile.Reset(new CWriteDB_ColumnData(dbname,
    //                                           extn2,
    //                                           index,
    //                                           max_file_size));
    //
    //     m_IFile.Reset(new CWriteDB_ColumnIndex(dbname,
    //                                            extn1,
    //                                            index,
    //                                            *m_DFile,
    //                                            title,
    //                                            meta,
    //                                            max_file_size));
    // }
    //
    // void CWriteDB_Column::AddByteOrder(const string      & dbname,
    //                               const string      & extn,
    //                               int                 index,
    //                               Uint8               max_file_size)
    // {
    //     m_UseBothByteOrder = true;
    //     m_DFile2.Reset(new CWriteDB_ColumnData(dbname,
    //                                           extn,
    //                                           index,
    //                                           max_file_size));
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:512-527
    // ```c++
    //     string extn(m_Protein ? "p??" : "n??");
    //     ...
    //     extn[1] = "abcdefghijklmnopqrstuvwxyz0123456789"[col_id];
    //
    //     string extn2 = extn;
    //     string extn3 = extn;
    //
    //     extn[2] = 'a';
    //     extn2[2] = 'b';
    //     extn3[2] = 'c';
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:199-199
    // ```c++
    //     m_Date = CTime(CTime::eCurrent).AsString();
    // ```
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
                file: WriteDbFile::deferred(dbname, extension('a'), index, max_file_size),
                title: MASK_COLUMN_TITLE.to_vec(),
                date: date.to_vec(),
                meta: meta.clone(),
                offsets: Vec::new(),
                oids: 0,
                data_length: 0,
            },
            data: ColumnData {
                file: WriteDbFile::deferred(dbname, extension('b'), index, max_file_size),
                data_length: 0,
            },
            data2: ColumnData {
                file: WriteDbFile::deferred(dbname, extension('c'), index, max_file_size),
                data_length: 0,
            },
            closed: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:116-119
    // ```c++
    // bool CWriteDB_Column::CanFit(int size) const
    // {
    //     return m_IFile->CanFit() && m_DFile->CanFit(size);
    // }
    // ```
    pub(crate) fn can_fit(&self, size: usize) -> bool {
        self.index.can_fit() && self.data.can_fit(size)
    }

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
    /// `AddBlob(blob, blob2)`; an empty blob adds an entry without data.
    pub(crate) fn add_blob(&mut self, blob: &[u8], blob2: &[u8]) -> Result<(), WriteDbError> {
        let data_size = self.data.write_blob(blob)?;
        self.index.write_blob_index(data_size);
        self.data2.write_blob(blob2)?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:109-114
    // ```c++
    // void CWriteDB_Column::Close()
    // {
    //     m_IFile->Close();
    //     m_DFile->Close();
    //     if (m_UseBothByteOrder) m_DFile2->Close();
    // }
    // ```
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.index.close(self.data.data_length)?;
        self.data.file.close()?;
        self.data2.file.close()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:121-133
    // ```c++
    // void CWriteDB_Column::RenameSingle()
    // {
    //     m_IFile->RenameSingle();
    //     m_DFile->RenameSingle();
    //     if (m_UseBothByteOrder) m_DFile2->RenameSingle();
    // }
    //
    // void CWriteDB_Column::RenameFileIndex(unsigned int num_digits)
    // {
    //     m_IFile->RenameFileIndex(num_digits);
    //     m_DFile->RenameFileIndex(num_digits);
    //     if (m_UseBothByteOrder) m_DFile2->RenameFileIndex(num_digits);
    // }
    // ```
    pub(crate) fn rename_single(&mut self) -> Result<(), WriteDbError> {
        self.index.file.rename_single()?;
        self.data.file.rename_single()?;
        self.data2.file.rename_single()
    }

    pub(crate) fn rename_file_index(&mut self, num_digits: u32) -> Result<(), WriteDbError> {
        self.index.file.rename_file_index(num_digits)?;
        self.data.file.rename_file_index(num_digits)?;
        self.data2.file.rename_file_index(num_digits)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:84-91
    // ```c++
    // void CWriteDB_Column::ListFiles(vector<string> & files, bool skip_empty) const
    // {
    //     if (! (skip_empty && m_DFile->Empty())) {
    //         files.push_back(m_IFile->GetFilename());
    //         files.push_back(m_DFile->GetFilename());
    //         if (m_UseBothByteOrder) files.push_back(m_DFile2->GetFilename());
    //     }
    // }
    // ```
    /// `ListFiles(files, true)`: the files of a column that holds data.
    pub(crate) fn list_files(&self, files: &mut Vec<String>) {
        if self.data.data_length != 0 {
            files.push(self.index.file.filename().to_string());
            files.push(self.data.file.filename().to_string());
            files.push(self.data2.file.filename().to_string());
        }
    }

    /// Drops the open files without writing what is still buffered.
    pub(crate) fn abandon(&mut self) {
        self.closed = true;
        self.index.file.abandon();
        self.data.file.abandon();
        self.data2.file.abandon();
    }
}
