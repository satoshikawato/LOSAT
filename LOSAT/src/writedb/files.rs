//! The three files of a volume as CWriteDB writes them
//! (`c++/src/objtools/blast/seqdb_writer/writedb_files.cpp` and
//! `c++/include/objtools/blast/seqdb_writer/writedb_files.hpp`): the index
//! file, the header file and the sequence file, each with NCBI's rule for
//! "does another sequence fit".

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use super::WriteDbError;

/// `CWriteDB_File`: one output file of a volume.
#[derive(Debug)]
pub(crate) struct WriteDbFile {
    writer: Option<BufWriter<File>>,
    base_name: String,
    extension: &'static str,
    index: u32,
    use_index: bool,
    fname: String,
    /// `m_Offset`: the stream position, an `unsigned int` in NCBI.
    offset: u32,
    max_file_size: u64,
    /// `m_Created`: whether the file exists. The index, header and sequence
    /// files of a volume are created with it (`always_create`); the files of
    /// a column only when they get data (`deferred`).
    created: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:235-257
// ```c++
// string CWriteDB_File::MakeShortName(const string & base, int index)
// {
//     ostringstream fns;
//
//     fns << base;
//     fns << ".";
//     fns << (index / 10);
//     fns << (index % 10);
//
//     return fns.str();
// }
//
// void CWriteDB_File::x_MakeFileName()
// {
//     if (m_UseIndex) {
//         m_Fname = MakeShortName(m_BaseName, m_Index);
//     } else {
//         m_Fname = m_BaseName;
//     }
//
//     m_Fname += ".";
//     m_Fname += m_Extension;
// }
// ```
/// `CWriteDB_File::MakeShortName`: `base.NN` (two digits, more from volume
/// 100 on).
pub(crate) fn make_short_name(base: &str, index: u32) -> String {
    format!("{base}.{}{}", index / 10, index % 10)
}

impl WriteDbFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:178-202
    // ```c++
    //     if (m_MaxFileSize == 0) {
    //         m_MaxFileSize = x_DefaultByteLimit();
    //     } else {
    // #ifdef _DEBUG
    //         _ASSERT(max_file_size <= MAX_FILE_SIZE);
    // #endif
    //     }
    //
    //     m_Nul.resize(1);
    //     m_Nul[0] = (char) 0;
    //
    //     m_UseIndex = (index >= 0);
    //     x_MakeFileName();
    //
    //     if (always_create) {
    //         Create();
    //     }
    // }
    //
    // void CWriteDB_File::Create()
    // {
    //     _ASSERT(! m_Created);
    //     m_Created = true;
    //     m_RealFile.open(m_Fname.c_str(), ios::out | ios::binary);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:174-178
    // ```c++
    //     Uint8 x_DefaultByteLimit()
    //     {
    //         // 1 gb (marketing version) - 1; about a billion
    //         return 1000*1000*1000 - 1;
    //     }
    // ```
    // NCBI does not check that the file opened or that a write succeeded;
    // LOSAT reports the failure instead of leaving a damaged database.
    fn create(
        base_name: &str,
        extension: &'static str,
        index: u32,
        max_file_size: u64,
    ) -> Result<Self, WriteDbError> {
        let fname = format!("{}.{extension}", make_short_name(base_name, index));
        let file = File::create(&fname).map_err(|source| WriteDbError::Io {
            path: PathBuf::from(&fname),
            source,
        })?;
        let mut created = Self::deferred(base_name, extension, index, max_file_size);
        created.writer = Some(BufWriter::with_capacity(1 << 16, file));
        created.created = true;
        Ok(created)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:178-194
    // ```c++
    //     if (m_MaxFileSize == 0) {
    //         m_MaxFileSize = x_DefaultByteLimit();
    //     }
    //     ...
    //     m_UseIndex = (index >= 0);
    //     x_MakeFileName();
    //
    //     if (always_create) {
    //         Create();
    //     }
    // ```
    /// `CWriteDB_File` with `always_create` false: the name is made, the file
    /// is not created (`create_file` does that when the first data arrive).
    pub(crate) fn deferred(
        base_name: &str,
        extension: &'static str,
        index: u32,
        max_file_size: u64,
    ) -> Self {
        Self {
            writer: None,
            base_name: base_name.to_string(),
            extension,
            index,
            use_index: true,
            fname: format!("{}.{extension}", make_short_name(base_name, index)),
            offset: 0,
            max_file_size: if max_file_size == 0 {
                1000 * 1000 * 1000 - 1
            } else {
                max_file_size
            },
            created: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:197-202
    // ```c++
    // void CWriteDB_File::Create()
    // {
    //     _ASSERT(! m_Created);
    //     m_Created = true;
    //     m_RealFile.open(m_Fname.c_str(), ios::out | ios::binary);
    // }
    // ```
    /// `Create` of a file made with `deferred`.
    pub(crate) fn create_file(&mut self) -> Result<(), WriteDbError> {
        if self.created {
            return Ok(());
        }
        let file = File::create(&self.fname).map_err(|source| WriteDbError::Io {
            path: PathBuf::from(&self.fname),
            source,
        })?;
        self.created = true;
        self.writer = Some(BufWriter::with_capacity(1 << 16, file));
        Ok(())
    }

    /// `m_Created`.
    pub(crate) fn is_created(&self) -> bool {
        self.created
    }

    /// `m_MaxFileSize`.
    pub(crate) fn max_file_size(&self) -> u64 {
        self.max_file_size
    }

    fn io_error(&self, source: std::io::Error) -> WriteDbError {
        WriteDbError::Io {
            path: PathBuf::from(&self.fname),
            source,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:218-221
    // ```c++
    //     m_RealFile.write(data.data(), data.length());
    //
    //     m_Offset += data.length();
    //     return m_Offset;
    // ```
    // `m_Offset` is an `unsigned int`: past 4 GiB it wraps and the offsets of
    // the index no longer name the data. makeblastdb limits `-max_file_sz`
    // below 4 GiB, but a volume takes its first sequence whatever its size;
    // LOSAT stops instead of writing offsets that wrapped.
    pub(crate) fn write(&mut self, data: &[u8]) -> Result<u32, WriteDbError> {
        let end = u64::from(self.offset) + data.len() as u64;
        if end > u64::from(u32::MAX) {
            return Err(WriteDbError::Creation(format!(
                "{} would pass 4 GiB, which the 32-bit offsets of a BLAST database volume \
                 cannot address",
                self.fname
            )));
        }
        let Some(writer) = self.writer.as_mut() else {
            return Err(self.io_error(std::io::Error::other("the file is closed")));
        };
        if let Err(source) = writer.write_all(data) {
            return Err(self.io_error(source));
        }
        self.offset = end as u32;
        Ok(self.offset)
    }

    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if let Some(mut writer) = self.writer.take() {
            if let Err(source) = writer.flush() {
                return Err(self.io_error(source));
            }
        }
        Ok(())
    }

    pub(crate) fn filename(&self) -> &str {
        &self.fname
    }

    /// Drops the file handle without writing what is still buffered.
    pub(crate) fn abandon(&mut self) {
        if let Some(writer) = self.writer.take() {
            let (_file, _unwritten) = writer.into_parts();
        }
    }

    /// Drops the file handle and removes the file: what a volume does with
    /// the files it created when one of its files cannot be created.
    pub(crate) fn discard(&mut self) {
        self.abandon();
        let _ = std::fs::remove_file(&self.fname);
    }

    // The error names both paths: the system's message does not say which
    // of the two it is about (a destination that is a directory, a source
    // that is gone).
    //
    // NCBI renames a file that was never created as well: `CDirEntry::Rename`
    // finds nothing to rename and the name changes all the same. LOSAT does
    // the same (the name of a file of a column without data is still the one
    // `list_files` would give).
    fn rename_to(&mut self, new_name: String) -> Result<(), WriteDbError> {
        if self.created {
            std::fs::rename(&self.fname, &new_name).map_err(|source| WriteDbError::Rename {
                from: PathBuf::from(&self.fname),
                to: PathBuf::from(&new_name),
                source,
            })?;
        }
        self.fname = new_name;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:267-277
    // ```c++
    // void CWriteDB_File::RenameSingle()
    // {
    //     _ASSERT(m_UseIndex == true);
    //
    //     string nm1 = m_Fname;
    //     m_UseIndex = false;
    //     x_MakeFileName();
    //
    //     CDirEntry fn1(nm1);
    //     fn1.Rename(m_Fname, CDirEntry::fRF_Overwrite);
    // }
    // ```
    pub(crate) fn rename_single(&mut self) -> Result<(), WriteDbError> {
        self.use_index = false;
        self.rename_to(format!("{}.{}", self.base_name, self.extension))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:279-304
    // ```c++
    // void CWriteDB_File::RenameFileIndex(unsigned int num_digits)
    // {
    // 	_ASSERT(num_digits > 2);
    // 	unsigned int orig_num_digits = log10(m_Index) +1;
    // 	if(orig_num_digits == num_digits) {
    // 		return;
    // 	}
    //
    //     string orig_fname = m_Fname;
    //     ostringstream fns;
    //
    //     fns << m_BaseName;
    //     fns << ".";
    //     for (unsigned int i=2; i< num_digits; i++){
    //     	fns << "0";
    //     }
    //     fns << (m_Index / 10);
    //     fns << (m_Index % 10);
    //     fns << ".";
    //     fns << m_Extension;
    //
    //     m_Fname = fns.str();
    //
    //     CDirEntry fn(orig_fname);
    //     fn.Rename(m_Fname, CDirEntry::fRF_Overwrite);
    // }
    // ```
    // `num_digits - 2` zeros go in front of the name `MakeShortName` gives,
    // whatever the number of digits of the index: with four digits volume 5
    // becomes `0005` and volume 123 becomes `00123`. For index 0 the
    // conversion of `log10(0) + 1` (minus infinity) to `unsigned int` gives 0
    // on the reference platform (x86-64), so volume 0 is renamed as well.
    pub(crate) fn rename_file_index(&mut self, num_digits: u32) -> Result<(), WriteDbError> {
        if orig_num_digits(self.index) == num_digits {
            return Ok(());
        }
        let zeros = "0".repeat(num_digits.saturating_sub(2) as usize);
        self.rename_to(format!(
            "{}.{zeros}{}{}.{}",
            self.base_name,
            self.index / 10,
            self.index % 10,
            self.extension
        ))
    }
}

/// `log10(m_Index) + 1` as an `unsigned int`: the number of decimal digits
/// of a positive index, and 0 for index 0 (see `rename_file_index`).
pub(crate) fn orig_num_digits(index: u32) -> u32 {
    if index == 0 {
        0
    } else {
        index.ilog10() + 1
    }
}

/// `CWriteDB_IndexFile`: the index file, kept in memory until the volume is
/// closed.
#[derive(Debug)]
pub(crate) struct IndexFile {
    file: WriteDbFile,
    protein: bool,
    title: Vec<u8>,
    date: Vec<u8>,
    oids: u32,
    data_size: u64,
    letters: u64,
    max_length: u32,
    hdr: Vec<u32>,
    seq: Vec<u32>,
    amb: Vec<u32>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:355-360
// ```c++
// int CWriteDB_IndexFile::x_Overhead(const string & T,
//                                    const string & D)
// {
//     return 4 * sizeof(int) + sizeof(long)
//             + 2 * sizeof(int) + T.size() + D.size();
// }
// ```
// `sizeof(long)` is 8 on the reference platform (Linux x86-64).
fn overhead(title: &[u8], date: &[u8]) -> u64 {
    (4 * 4 + 8 + 2 * 4 + title.len() + date.len()) as u64
}

impl IndexFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:329-344
    // ```c++
    //     m_Overhead = x_Overhead(title, date);
    //     if (dbver == eBDB_Version5) {
    //         m_Overhead = x_Overhead(title, x_MakeLmdbName(), date);
    //     } else {
    //         m_Overhead = x_Overhead(title, date);
    //     }
    //     m_Overhead = s_RoundUp(m_Overhead, 8);
    //     m_DataSize = m_Overhead;
    //
    //     // The '1' added to the sequence offset array refers to the fact
    //     // that sequence files contain an initial NUL byte.  This seems to
    //     // be for the benefit of the protein database scanning code, but
    //     // it is also done for nucleotide databases.
    //
    //     m_Hdr.push_back(0);
    //     m_Seq.push_back(1);
    // ```
    pub(crate) fn create(
        dbname: &str,
        protein: bool,
        title: &[u8],
        date: &[u8],
        index: u32,
        max_file_size: u64,
    ) -> Result<Self, WriteDbError> {
        let file = WriteDbFile::create(
            dbname,
            if protein { "pin" } else { "nin" },
            index,
            max_file_size,
        )?;
        Ok(Self {
            file,
            protein,
            title: title.to_vec(),
            date: date.to_vec(),
            oids: 0,
            data_size: overhead(title, date).div_ceil(8) * 8,
            letters: 0,
            max_length: 0,
            hdr: vec![0],
            seq: vec![1],
            amb: Vec::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:228-236
    // ```c++
    //     bool CanFit()
    //     {
    //         _ASSERT(m_MaxFileSize > 1024UL);
    //
    //         if (m_OIDs == 0)
    //             return true;
    //
    //         return m_DataSize < (m_MaxFileSize - 12UL);
    //     }
    // ```
    // The subtraction is unsigned: a limit under 12 bytes wraps and always
    // fits.
    pub(crate) fn can_fit(&self) -> bool {
        if self.oids == 0 {
            return true;
        }
        self.data_size < self.file.max_file_size.wrapping_sub(12)
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:246-258
    // ```c++
    //     void AddSequence(int length, unsigned int hdr, unsigned int seq)
    //     {
    //         if (length > m_MaxLength) {
    //             m_MaxLength = length;
    //         }
    //
    //         m_OIDs++;
    //         m_Letters += length;
    //         m_DataSize += 8;
    //
    //         m_Hdr.push_back(hdr);
    //         m_Seq.push_back(seq);
    //     }
    // ```
    /// Adds a protein: its length and the ends of its header and sequence.
    pub(crate) fn add_protein(&mut self, length: u32, hdr: u32, seq: u32) {
        self.max_length = self.max_length.max(length);
        self.oids += 1;
        self.letters += u64::from(length);
        self.data_size += 8;
        self.hdr.push(hdr);
        self.seq.push(seq);
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:269-282
    // ```c++
    //     void AddSequence(int length, unsigned int hdr, unsigned int seq, unsigned int amb)
    //     {
    //         if (length > m_MaxLength) {
    //             m_MaxLength = length;
    //         }
    //
    //         m_OIDs++;
    //         m_Letters += length;
    //
    //         m_DataSize += 12;
    //         m_Hdr.push_back(hdr);
    //         m_Seq.push_back(amb); // Not a bug.
    //         m_Amb.push_back(seq); // Also not a bug.
    //     }
    // ```
    /// Adds a nucleotide sequence: `seq` is the end of its 2-bit data (where
    /// its ambiguity data start) and `amb` the end of its ambiguity data
    /// (where the next sequence starts).
    pub(crate) fn add_nucleotide(&mut self, length: u32, hdr: u32, seq: u32, amb: u32) {
        self.max_length = self.max_length.max(length);
        self.oids += 1;
        self.letters += u64::from(length);
        self.data_size += 12;
        self.hdr.push(hdr);
        self.seq.push(amb);
        self.amb.push(seq);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:373-430
    // ```c++
    //     string pad_date = m_Date;
    //     int count = 0;
    //     const string lmdb_name = use_lmdb ? x_MakeLmdbName() : "";
    //     int overhead = use_lmdb
    //             ? x_Overhead(m_Title, lmdb_name, pad_date)
    //             : x_Overhead(m_Title, pad_date);
    //     while (overhead & 0x7) {
    //         pad_date.append(m_Nul);
    //         if (count != -1) {
    //             _ASSERT(count++ < 8);
    //         }
    //         overhead = use_lmdb
    //                 ? x_Overhead(m_Title, lmdb_name, pad_date)
    //                 : x_Overhead(m_Title, pad_date);
    //     }
    //
    //     // Write header
    //
    //     ostream & F = m_RealFile;
    //
    //     s_WriteInt4  (F, format_version);
    //     s_WriteInt4  (F, seq_type);
    //     if (!lmdb_name.empty()) {
    //         s_WriteInt4  (F, m_Index);
    //         s_WriteString(F, m_Title);
    //         s_WriteString(F, lmdb_name);
    //     } else {
    //         s_WriteString(F, m_Title);
    //     }
    //     s_WriteString(F, pad_date);
    //     s_WriteInt4  (F, m_OIDs);
    //     s_WriteInt8LE(F, m_Letters);
    //     s_WriteInt4  (F, m_MaxLength);
    //
    //     for(unsigned i = 0; i < m_Hdr.size(); i++) {
    //         s_WriteInt4(F, m_Hdr[i]);
    //     }
    //
    //     for(unsigned i = 0; i < m_Seq.size(); i++) {
    //         s_WriteInt4(F, m_Seq[i]);
    //     }
    //
    //     // Should loop m_OID times, or not at all.
    //     for(unsigned i = 0; i < m_Amb.size(); i++) {
    //         s_WriteInt4(F, m_Amb[i]);
    //     }
    //     ...
    //     if (m_Amb.size()) {
    //         s_WriteInt4(F, m_Seq.back());
    //     }
    // ```
    /// `x_Flush` (format version 4) and `Close`.
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        let mut pad_date = self.date.clone();
        while overhead(&self.title, &pad_date) & 0x7 != 0 {
            pad_date.push(0);
        }
        let mut out = Vec::with_capacity(
            32 + self.title.len() + pad_date.len() + 4 * (self.hdr.len() * 3 + 1),
        );
        let write_int4 =
            |out: &mut Vec<u8>, value: u32| out.extend_from_slice(&value.to_be_bytes());
        write_int4(&mut out, 4);
        write_int4(&mut out, u32::from(self.protein));
        write_int4(&mut out, self.title.len() as u32);
        out.extend_from_slice(&self.title);
        write_int4(&mut out, pad_date.len() as u32);
        out.extend_from_slice(&pad_date);
        write_int4(&mut out, self.oids);
        out.extend_from_slice(&self.letters.to_le_bytes());
        write_int4(&mut out, self.max_length);
        for array in [&self.hdr, &self.seq, &self.amb] {
            for &value in array {
                write_int4(&mut out, value);
            }
        }
        if !self.amb.is_empty() {
            let last = *self.seq.last().unwrap_or(&1);
            write_int4(&mut out, last);
        }
        self.file.write(&out)?;
        self.hdr = Vec::new();
        self.seq = Vec::new();
        self.amb = Vec::new();
        self.file.close()
    }

    pub(crate) fn file(&self) -> &WriteDbFile {
        &self.file
    }

    pub(crate) fn file_mut(&mut self) -> &mut WriteDbFile {
        &mut self.file
    }
}

/// `CWriteDB_HeaderFile`.
#[derive(Debug)]
pub(crate) struct HeaderFile {
    file: WriteDbFile,
    data_size: u64,
}

impl HeaderFile {
    pub(crate) fn create(
        dbname: &str,
        protein: bool,
        index: u32,
        max_file_size: u64,
    ) -> Result<Self, WriteDbError> {
        let file = WriteDbFile::create(
            dbname,
            if protein { "phr" } else { "nhr" },
            index,
            max_file_size,
        )?;
        Ok(Self { file, data_size: 0 })
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:365-382
    // ```c++
    //     bool CanFit(int size)
    //     {
    //         _ASSERT(size >= 0);
    //
    //         if (m_DataSize == 0UL) {
    //             return true;
    //         }
    //
    //         return (m_DataSize + (Uint8) size) < m_MaxFileSize;
    //     }
    //
    //     /// Add binary header data to this file.
    //     /// @param binhdr Binary ASN.1 version of header data. [in]
    //     /// @param offset Offset of end of header data. [out]
    //     void AddSequence(const string & binhdr, unsigned int & offset)
    //     {
    //         m_DataSize = offset = Write(binhdr);
    //     }
    // ```
    pub(crate) fn can_fit(&self, size: usize) -> bool {
        if self.data_size == 0 {
            return true;
        }
        self.data_size + (size as u64) < self.file.max_file_size
    }

    /// Writes one header record and returns the offset of its end.
    pub(crate) fn add_sequence(&mut self, binhdr: &[u8]) -> Result<u32, WriteDbError> {
        let offset = self.file.write(binhdr)?;
        self.data_size = u64::from(offset);
        Ok(offset)
    }

    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        self.file.close()
    }

    pub(crate) fn file(&self) -> &WriteDbFile {
        &self.file
    }

    pub(crate) fn file_mut(&mut self) -> &mut WriteDbFile {
        &mut self.file
    }
}

/// `CWriteDB_SequenceFile`.
#[derive(Debug)]
pub(crate) struct SequenceFile {
    file: WriteDbFile,
}

impl SequenceFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:481-487
    // ```c++
    //     // Only protein sequences need the inter-sequence NUL bytes.
    //     // The first null written here is for nucleotide sequences.
    //     // It doesn't seem necessary, but formatdb provides it, so I
    //     // will too.
    //
    //     WriteWithNull(string());
    // }
    // ```
    pub(crate) fn create(
        dbname: &str,
        protein: bool,
        index: u32,
        max_file_size: u64,
    ) -> Result<Self, WriteDbError> {
        let mut file = WriteDbFile::create(
            dbname,
            if protein { "psq" } else { "nsq" },
            index,
            max_file_size,
        )?;
        file.write(&[0])?;
        Ok(Self { file })
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:419-434
    // ```c++
    //     bool CanFit(int size, int letters)
    //     {
    //         _ASSERT(size >= 0);
    //         _ASSERT(letters >= 0);
    //
    //         if (m_Offset <= 1) {
    //             return true;
    //         }
    //
    //         if ((m_BaseLimit != 0)  &&
    //             ((m_Letters + (Uint8) letters) > m_BaseLimit)) {
    //             return false;
    //         }
    //
    //         return ((Uint8)(m_Offset + size) < m_MaxFileSize);
    //     }
    // ```
    // makeblastdb sets no limit on the letters of a volume (`m_BaseLimit` is
    // 0). `m_Offset + size` is a 32-bit unsigned sum, as in NCBI.
    pub(crate) fn can_fit(&self, size: usize) -> bool {
        if self.file.offset <= 1 {
            return true;
        }
        u64::from(self.file.offset.wrapping_add(size as u32)) < self.file.max_file_size
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:443-452
    // ```c++
    //     void AddSequence(const string & sequence,
    //                      unsigned int & offset,
    //                      int            length)
    //     {
    // #ifdef _DEBUG
    //         _ASSERT(m_Protein);
    // #endif
    //         offset = WriteWithNull(sequence);
    //         m_Letters += length;
    //     }
    // ```
    /// Writes a protein and the NUL after it; returns the offset after the
    /// NUL.
    pub(crate) fn add_protein(&mut self, sequence: &[u8]) -> Result<u32, WriteDbError> {
        self.file.write(sequence)?;
        self.file.write(&[0])
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_files.hpp:463-475
    // ```c++
    //     void AddSequence(const string & sequence,
    //                      const string & ambig,
    //                      unsigned int & off_seq,
    //                      unsigned int & off_amb,
    //                      int          length)
    //     {
    // #ifdef _DEBUG
    //         _ASSERT(! m_Protein);
    // #endif
    //         off_seq = Write(sequence);
    //         off_amb = Write(ambig);
    //         m_Letters += length;
    //     }
    // ```
    /// Writes a nucleotide sequence and its ambiguity data; returns the
    /// offsets after each.
    pub(crate) fn add_nucleotide(
        &mut self,
        sequence: &[u8],
        ambig: &[u8],
    ) -> Result<(u32, u32), WriteDbError> {
        let off_seq = self.file.write(sequence)?;
        let off_amb = self.file.write(ambig)?;
        Ok((off_seq, off_amb))
    }

    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        self.file.close()
    }

    pub(crate) fn file(&self) -> &WriteDbFile {
        &self.file
    }

    pub(crate) fn file_mut(&mut self) -> &mut WriteDbFile {
        &mut self.file
    }
}
