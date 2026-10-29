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
        Ok(Self {
            writer: Some(BufWriter::with_capacity(1 << 16, file)),
            base_name: base_name.to_string(),
            extension,
            index,
            use_index: true,
            fname,
            offset: 0,
            max_file_size: if max_file_size == 0 {
                1000 * 1000 * 1000 - 1
            } else {
                max_file_size
            },
        })
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
    fn write(&mut self, data: &[u8]) -> Result<u32, WriteDbError> {
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

    fn close(&mut self) -> Result<(), WriteDbError> {
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
    fn rename_to(&mut self, new_name: String) -> Result<(), WriteDbError> {
        std::fs::rename(&self.fname, &new_name).map_err(|source| WriteDbError::Rename {
            from: PathBuf::from(&self.fname),
            to: PathBuf::from(&new_name),
            source,
        })?;
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

/// `EBlastDbVersion`: the format version of the database being written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbVersion {
    /// `eBDB_Version4`: the volumes only.
    V4,
    /// `eBDB_Version5`: the volumes, whose index files name the LMDB file, and the LMDB
    /// files of the database (`super::lmdb`).
    V5,
}

impl DbVersion {
    /// The version of `-blastdb_version` (makeblastdb allows 4 and 5 only).
    pub fn from_number(version: i32) -> Option<Self> {
        match version {
            4 => Some(DbVersion::V4),
            5 => Some(DbVersion::V5),
            _ => None,
        }
    }

    /// `(int) m_Version`: the number written as the format version of an index file.
    pub fn number(self) -> u32 {
        match self {
            DbVersion::V4 => 4,
            DbVersion::V5 => 5,
        }
    }
}

/// The offset arrays of a closed index file, emptied, for the index file of
/// the next volume.
///
/// NCBI frees the arrays of a volume when it writes them and the next volume
/// grows its own from nothing. With glibc that puts the allocator's history
/// into the peak: the free of a small first volume's arrays (mapped blocks
/// of 1-32 MiB) raises the mmap threshold to their size, the next volumes
/// grow theirs in the heap up to that size, and the heap copies they leave
/// when they move to a mapping stay resident unless they exceed the trim
/// threshold (UniRef50, 5 volumes: 43 MiB stayed for the rest of the build).
/// Arrays that are handed on are grown in place and freed only when the
/// database is closed, so no free moves the thresholds, and their pages are
/// not faulted in again for every volume.
#[derive(Debug, Default)]
pub(crate) struct OffsetArrays {
    hdr: Vec<u32>,
    seq: Vec<u32>,
    amb: Vec<u32>,
}

/// `CWriteDB_IndexFile`: the index file, kept in memory until the volume is
/// closed.
#[derive(Debug)]
pub(crate) struct IndexFile {
    file: WriteDbFile,
    protein: bool,
    title: Vec<u8>,
    date: Vec<u8>,
    /// `m_Version`.
    version: DbVersion,
    /// `m_Index`: the number of the volume, written into a version 5 index file.
    index: u32,
    /// `x_MakeLmdbName()` for version 5, empty for version 4.
    lmdb_name: Vec<u8>,
    /// The platform's `sizeof(long)` in `x_Overhead`.
    sizeof_long: u64,
    oids: u32,
    data_size: u64,
    letters: u64,
    max_length: u32,
    hdr: Vec<u32>,
    seq: Vec<u32>,
    amb: Vec<u32>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:346-360
// ```c++
// int CWriteDB_IndexFile::x_Overhead(const string & T,
//                                    const string & lmdbName,
//                                    const string & D)
// {
//     return 5 * sizeof(int) + sizeof(long)
//             + 3 * sizeof(int) + T.size() + lmdbName.size() + D.size();
// }
//
// int CWriteDB_IndexFile::x_Overhead(const string & T,
//                                    const string & D)
// {
//     return 4 * sizeof(int) + sizeof(long)
//             + 2 * sizeof(int) + T.size() + D.size();
// }
// ```
// `sizeof(int)` is 4 and `sizeof(long)` the platform's (`Platform::sizeof_long`):
// 8 on Linux and macOS (LP64), 4 on Windows (LLP64), so the padding of the
// date and the room the header takes from `-max_file_sz` follow the platform
// whose makeblastdb is reproduced. The version 5 form counts the volume number
// and the length of the LMDB name as two more `int`s.
fn overhead(
    version: DbVersion,
    sizeof_long: u64,
    title: &[u8],
    lmdb_name: &[u8],
    date: &[u8],
) -> u64 {
    match version {
        DbVersion::V4 => 4 * 4 + sizeof_long + 2 * 4 + (title.len() + date.len()) as u64,
        DbVersion::V5 => {
            5 * 4 + sizeof_long + 3 * 4 + (title.len() + lmdb_name.len() + date.len()) as u64
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:438-448
// ```c++
// /// Form name of lmdb database file.
// const string CWriteDB_IndexFile::x_MakeLmdbName()
// {
//     string suffix = (m_Protein ? ".pdb" : ".ndb");
//     size_t last_slash = m_BaseName.find_last_of(CFile::GetPathSeparator());
//     if (last_slash == m_BaseName.npos) {
//         return m_BaseName + suffix;
//     } else {
//         return m_BaseName.substr(last_slash + 1) + suffix;
//     }
// }
// ```
// `m_BaseName` is the database name without the number of the volume, and
// `CFile::GetPathSeparator()` is the separator of the platform ('/', or '\\'
// on Windows).
/// `x_MakeLmdbName`: the file name of the database's LMDB file, without its
/// directory.
pub(crate) fn make_lmdb_name(base_name: &str, protein: bool) -> String {
    let suffix = if protein { ".pdb" } else { ".ndb" };
    let name = match base_name.rfind(std::path::MAIN_SEPARATOR) {
        Some(last_slash) => &base_name[last_slash + 1..],
        None => base_name,
    };
    format!("{name}{suffix}")
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
    // `m_DataSize` starts at the header's size, so the LMDB name of version 5
    // takes room from the limit of `-max_file_sz` (`can_fit`).
    /// `arrays` are those of the previous volume (`take_offset_arrays`), or
    /// new ones for the first volume.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn create(
        dbname: &str,
        protein: bool,
        title: &[u8],
        date: &[u8],
        index: u32,
        max_file_size: u64,
        version: DbVersion,
        sizeof_long: u64,
        arrays: OffsetArrays,
    ) -> Result<Self, WriteDbError> {
        let file = WriteDbFile::create(
            dbname,
            if protein { "pin" } else { "nin" },
            index,
            max_file_size,
        )?;
        let lmdb_name = match version {
            DbVersion::V4 => Vec::new(),
            DbVersion::V5 => make_lmdb_name(dbname, protein).into_bytes(),
        };
        let OffsetArrays {
            mut hdr,
            mut seq,
            mut amb,
        } = arrays;
        hdr.clear();
        seq.clear();
        amb.clear();
        hdr.push(0);
        seq.push(1);
        Ok(Self {
            file,
            protein,
            title: title.to_vec(),
            date: date.to_vec(),
            version,
            index,
            data_size: overhead(version, sizeof_long, title, &lmdb_name, date).div_ceil(8) * 8,
            lmdb_name,
            sizeof_long,
            oids: 0,
            letters: 0,
            max_length: 0,
            hdr,
            seq,
            amb,
        })
    }

    /// `m_Letters`: the residues or bases of the volume.
    pub(crate) fn letters(&self) -> u64 {
        self.letters
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
    /// `x_Flush` and `Close`.
    // NCBI writes the offsets one by one to the file's stream. They go out
    // here through a buffer of a fixed size: a copy of the whole index next
    // to the offset arrays would double what a volume of many sequences
    // holds at its close (UniRef50: 16.6 million sequences in one volume,
    // 127 MiB of arrays and 127 MiB of copy; makedb-v5 survey/MEM_RCA.md).
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        let mut pad_date = self.date.clone();
        while overhead(
            self.version,
            self.sizeof_long,
            &self.title,
            &self.lmdb_name,
            &pad_date,
        ) & 0x7
            != 0
        {
            pad_date.push(0);
        }
        let mut out =
            Vec::with_capacity(40 + self.title.len() + self.lmdb_name.len() + pad_date.len());
        let write_int4 =
            |out: &mut Vec<u8>, value: u32| out.extend_from_slice(&value.to_be_bytes());
        write_int4(&mut out, self.version.number());
        write_int4(&mut out, u32::from(self.protein));
        if !self.lmdb_name.is_empty() {
            write_int4(&mut out, self.index);
            write_int4(&mut out, self.title.len() as u32);
            out.extend_from_slice(&self.title);
            write_int4(&mut out, self.lmdb_name.len() as u32);
            out.extend_from_slice(&self.lmdb_name);
        } else {
            write_int4(&mut out, self.title.len() as u32);
            out.extend_from_slice(&self.title);
        }
        write_int4(&mut out, pad_date.len() as u32);
        out.extend_from_slice(&pad_date);
        write_int4(&mut out, self.oids);
        out.extend_from_slice(&self.letters.to_le_bytes());
        write_int4(&mut out, self.max_length);
        self.file.write(&out)?;
        let mut chunk = [0u8; 1 << 14];
        for array in [&self.hdr, &self.seq, &self.amb] {
            for values in array.chunks(chunk.len() / 4) {
                let bytes = &mut chunk[..values.len() * 4];
                for (out, value) in bytes.chunks_exact_mut(4).zip(values) {
                    out.copy_from_slice(&value.to_be_bytes());
                }
                self.file.write(bytes)?;
            }
        }
        if !self.amb.is_empty() {
            let last = *self.seq.last().unwrap_or(&1);
            self.file.write(&last.to_be_bytes())?;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:432-435
        // ```c++
        //     vector<unsigned int> tmp1, tmp2, tmp3;
        //     m_Hdr.swap(tmp1);
        //     m_Seq.swap(tmp2);
        //     m_Amb.swap(tmp3);
        // ```
        // NCBI frees the arrays here; LOSAT empties them and keeps their
        // memory for the next volume's index file (`OffsetArrays`).
        self.hdr.clear();
        self.seq.clear();
        self.amb.clear();
        self.file.close()
    }

    /// The offset arrays, emptied by `close`, for the next volume.
    pub(crate) fn take_offset_arrays(&mut self) -> OffsetArrays {
        OffsetArrays {
            hdr: std::mem::take(&mut self.hdr),
            seq: std::mem::take(&mut self.seq),
            amb: std::mem::take(&mut self.amb),
        }
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
