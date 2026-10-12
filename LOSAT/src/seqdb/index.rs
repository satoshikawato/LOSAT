//! Index file of a database volume (`.pin` / `.nin`), format versions 4 and 5.

use std::path::Path;

use super::{SeqDbError, SeqType};

/// Parsed index file: the volume's summary values and its offset arrays.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexFile {
    /// Format version, 4 or 5.
    pub version: u32,
    /// Sequence type stored in the file.
    pub seqtype: SeqType,
    /// Volume number (version 5 only).
    pub volume_number: Option<u32>,
    /// Volume title, as stored.
    pub title: Vec<u8>,
    /// Name of the LMDB file (version 5 only). LOSAT does not read that file.
    pub lmdb_file: Option<Vec<u8>>,
    /// Creation date, as stored (may end in NUL padding).
    pub date: Vec<u8>,
    /// Number of sequences (ordinal ids) in the volume.
    pub num_oids: u32,
    /// Total number of residues or bases in the volume.
    pub volume_length: u64,
    /// Length of the longest sequence.
    pub max_length: u32,
    hdr_offsets: Vec<u32>,
    seq_offsets: Vec<u32>,
    amb_offsets: Vec<u32>,
}

struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    path: &'a Path,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, len: usize) -> Result<&'a [u8], SeqDbError> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|&end| end <= self.bytes.len());
        let Some(end) = end else {
            return Err(SeqDbError::Format {
                path: self.path.to_path_buf(),
                message: format!(
                    "index file ends inside a field ({} bytes wanted at offset {}, file has {})",
                    len,
                    self.pos,
                    self.bytes.len()
                ),
            });
        };
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:71-79
    // ```c++
    // TIndx CSeqDBRawFile::ReadSwapped(CSeqDBFileMemMap & lease,
    //                                  TIndx            offset,
    //                                  Uint4          * value) const
    //
    // {
    //     *value = SeqDB_GetStdOrd((Uint4 *) lease.GetFileDataPtr(m_FileName,offset));
    //
    //     return offset + sizeof(*value);
    // }
    // ```
    fn be_u32(&mut self) -> Result<u32, SeqDbError> {
        let bytes = self.take(4)?;
        Ok(u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:81-89
    // ```c++
    // TIndx CSeqDBRawFile::ReadSwapped(CSeqDBFileMemMap & lease,
    //                                  TIndx            offset,
    //                                  Uint8          * value) const
    //
    // {
    //     *value = SeqDB_GetBroken((Int8 *) lease.GetFileDataPtr(m_FileName,offset));
    //
    //     return offset + sizeof(*value);
    // }
    // ```
    // `SeqDB_GetBroken` reads the one little-endian field of the format.
    fn broken_u64(&mut self) -> Result<u64, SeqDbError> {
        let bytes = self.take(8)?;
        let mut raw = [0u8; 8];
        raw.copy_from_slice(bytes);
        Ok(u64::from_le_bytes(raw))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:91-104
    // ```c++
    // TIndx CSeqDBRawFile::ReadSwapped(CSeqDBFileMemMap & lease,
    //                                  TIndx            offset,
    //                                  string         * value) const
    // {
    //     Uint4 len = 0;
    //
    //     len = SeqDB_GetStdOrd((Int4 *) lease.GetFileDataPtr(m_FileName,offset));
    //
    //     offset += sizeof(len);
    //
    //     value->assign(lease.GetFileDataPtr(offset), (int) len);
    //
    //     return offset + len;
    // }
    // ```
    fn string(&mut self) -> Result<Vec<u8>, SeqDbError> {
        let len = self.be_u32()? as usize;
        Ok(self.take(len)?.to_vec())
    }

    fn offsets(&mut self, count: usize) -> Result<Vec<u32>, SeqDbError> {
        let bytes = self.take(count.checked_mul(4).ok_or_else(|| SeqDbError::Format {
            path: self.path.to_path_buf(),
            message: "index file offset array size overflows".to_string(),
        })?)?;
        Ok(bytes
            .chunks_exact(4)
            .map(|word| u32::from_be_bytes([word[0], word[1], word[2], word[3]]))
            .collect())
    }
}

impl IndexFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:175-235
    // ```c++
    //     offset = x_ReadSwapped(m_Lease, offset, & f_format_version);
    //     ...
    //     	 bool dbv5 = ( 5 == f_format_version );
    //         if (!dbv5 && (f_format_version != 4)) {
    //             NCBI_THROW(CSeqDBException,
    //                        eFileErr,
    //                        "Error: Not a valid version 4 or 5 database.");
    //         }
    //
    //         offset = x_ReadSwapped(m_Lease, offset, & f_db_seqtype);
    //         if (dbv5) {
    //         	offset = x_ReadSwapped(m_Lease, offset, & m_Volume);
    //         }
    //
    //         offset = x_ReadSwapped(m_Lease, offset, & m_Title);
    //         if (dbv5) {
    //             offset = x_ReadSwapped(m_Lease, offset, & m_LMDBFile);
    //        }
    //
    //         offset = x_ReadSwapped(m_Lease, offset, & m_Date);
    //         ...
    //         offset = x_ReadSwapped(m_Lease, offset, & m_NumOIDs);
    //         ...
    //         offset = x_ReadSwapped(m_Lease, offset, & m_VolLen);
    //         ...
    //         offset = x_ReadSwapped(m_Lease, offset, & m_MaxLen);
    //
    //         TIndx region_bytes = 4 * (m_NumOIDs + 1);
    //
    //         off1   = offset;
    //         off2   = off1 + region_bytes;
    //         off3   = off2 + region_bytes;
    //         offend = off3 + region_bytes;
    //     ...
    //     char db_seqtype = ((f_db_seqtype == 1) ? 'p' : 'n');
    //
    //     if (db_seqtype != x_GetSeqType()) {
    //         NCBI_THROW(CSeqDBException,
    //                    eFileErr,
    //                    "Error: requested sequence type does not match DB.");
    //     }
    //
    //     m_OffHdr = off1; m_EndHdr = off2;
    //     m_OffSeq = off2; m_EndSeq = off3;
    //
    //     if (db_seqtype == 'n') {
    //         m_OffAmb = off3; m_EndAmb = offend;
    //     } else {
    //         m_OffAmb = m_EndAmb = 0;
    //     }
    // ```
    // NCBI keeps the arrays in the mapped file; LOSAT copies them (4 or 8, for
    // nucleotide 12, bytes per sequence). Arrays that are cut short or not in
    // ascending order, which NCBI would follow out of its mapping, are
    // explicit format errors.
    /// Parses the bytes of an index file opened as `seqtype`.
    pub fn parse(bytes: &[u8], seqtype: SeqType, path: &Path) -> Result<Self, SeqDbError> {
        let format = |message: &str| SeqDbError::Format {
            path: path.to_path_buf(),
            message: message.to_string(),
        };
        let mut cursor = Cursor {
            bytes,
            pos: 0,
            path,
        };
        let version = cursor.be_u32()?;
        let dbv5 = version == 5;
        if !dbv5 && version != 4 {
            return Err(format("Error: Not a valid version 4 or 5 database."));
        }
        let db_seqtype = cursor.be_u32()?;
        let volume_number = if dbv5 { Some(cursor.be_u32()?) } else { None };
        let title = cursor.string()?;
        let lmdb_file = if dbv5 { Some(cursor.string()?) } else { None };
        let date = cursor.string()?;
        let num_oids = cursor.be_u32()?;
        let volume_length = cursor.broken_u64()?;
        let max_length = cursor.be_u32()?;

        let file_seqtype = if db_seqtype == 1 {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        if file_seqtype != seqtype {
            return Err(format("Error: requested sequence type does not match DB."));
        }
        // NCBI holds the count in a signed 32-bit `int`.
        if num_oids > i32::MAX as u32 {
            return Err(format(
                "index file holds more sequences than a signed 32-bit count",
            ));
        }

        let region = num_oids as usize + 1;
        let hdr_offsets = cursor.offsets(region)?;
        let seq_offsets = cursor.offsets(region)?;
        let amb_offsets = match seqtype {
            SeqType::Nucleotide => cursor.offsets(region)?,
            SeqType::Protein => Vec::new(),
        };

        let index = Self {
            version,
            seqtype,
            volume_number,
            title,
            lmdb_file,
            date,
            num_oids,
            volume_length,
            max_length,
            hdr_offsets,
            seq_offsets,
            amb_offsets,
        };
        index.check_offsets(path)?;
        Ok(index)
    }

    fn check_offsets(&self, path: &Path) -> Result<(), SeqDbError> {
        let format = |message: String| SeqDbError::Format {
            path: path.to_path_buf(),
            message,
        };
        if let Some(oid) = self.hdr_offsets.windows(2).position(|w| w[0] > w[1]) {
            return Err(format(format!(
                "header offsets of sequence {oid} are not in ascending order"
            )));
        }
        if let Some(oid) = self.seq_offsets.windows(2).position(|w| w[0] > w[1]) {
            return Err(format(format!(
                "sequence offsets of sequence {oid} are not in ascending order"
            )));
        }
        match self.seqtype {
            SeqType::Protein => {
                // One NUL byte precedes the first sequence and one follows each.
                if self.num_oids > 0 && self.seq_offsets[0] == 0 {
                    return Err(format(
                        "the first protein sequence has no sentinel byte before it".to_string(),
                    ));
                }
                if let Some(oid) = self.seq_offsets.windows(2).position(|w| w[0] == w[1]) {
                    return Err(format(format!(
                        "protein sequence {oid} has no sentinel byte after it"
                    )));
                }
            }
            SeqType::Nucleotide => {
                for oid in 0..self.num_oids as usize {
                    let start = self.seq_offsets[oid];
                    let amb = self.amb_offsets[oid];
                    let next = self.seq_offsets[oid + 1];
                    // The 2-bit data hold at least the byte with the remainder count.
                    if !(start < amb && amb <= next) {
                        return Err(format(format!(
                            "ambiguity offset of sequence {oid} lies outside its sequence data"
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:664-669
    // ```c++
    // void
    // CSeqDBIdxFile::GetHdrStartEnd(int oid, TIndx & start, TIndx & end) const
    // {
    //     start = SeqDB_GetStdOrd(& x_GetHdr()[oid]);
    //     end   = SeqDB_GetStdOrd(& x_GetHdr()[oid+1]);
    // }
    // ```
    /// Byte range of the header (`Blast-def-line-set`) of `oid` in the header file.
    pub fn header_range(&self, oid: u32) -> (u64, u64) {
        let oid = oid as usize;
        (
            u64::from(self.hdr_offsets[oid]),
            u64::from(self.hdr_offsets[oid + 1]),
        )
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:671-681
    // ```c++
    // void
    // CSeqDBIdxFile::GetSeqStartEnd(int oid, TIndx & start, TIndx & end) const
    // {
    //     start = SeqDB_GetStdOrd(& x_GetSeq()[oid]);
    //
    //     if ('p' == x_GetSeqType()) {
    //         end = SeqDB_GetStdOrd(& x_GetSeq()[oid+1]);
    //     } else {
    //         end = SeqDB_GetStdOrd(& x_GetAmb()[oid]);
    //     }
    // }
    // ```
    /// Byte range of the sequence data of `oid` in the sequence file. For a
    /// protein the range includes the NUL byte that follows the residues; for
    /// a nucleotide it ends where the ambiguity data start.
    pub fn sequence_range(&self, oid: u32) -> (u64, u64) {
        let oid = oid as usize;
        let start = u64::from(self.seq_offsets[oid]);
        let end = match self.seqtype {
            SeqType::Protein => u64::from(self.seq_offsets[oid + 1]),
            SeqType::Nucleotide => u64::from(self.amb_offsets[oid]),
        };
        (start, end)
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:651-662
    // ```c++
    // bool
    // CSeqDBIdxFile::GetAmbStartEnd(int oid, TIndx & start, TIndx & end) const
    // {
    //     if ('n' == x_GetSeqType()) {
    //         start = SeqDB_GetStdOrd(& x_GetAmb()[oid]);
    //         end   = SeqDB_GetStdOrd(& x_GetSeq()[oid+1]);
    //
    //         return (start <= end);
    //     }
    //
    //     return false;
    // }
    // ```
    /// Byte range of the ambiguity data of `oid` (nucleotide volumes only).
    pub fn ambiguity_range(&self, oid: u32) -> Option<(u64, u64)> {
        match self.seqtype {
            SeqType::Protein => None,
            SeqType::Nucleotide => {
                let oid = oid as usize;
                Some((
                    u64::from(self.amb_offsets[oid]),
                    u64::from(self.seq_offsets[oid + 1]),
                ))
            }
        }
    }

    /// Offset in the sequence file where the data of `oid` start; `oid` may be
    /// `num_oids`, which gives the end of the last sequence's data.
    pub fn sequence_start(&self, oid: u32) -> u64 {
        u64::from(self.seq_offsets[oid as usize])
    }

    /// The sequence start offsets, `num_oids + 1` entries.
    pub fn sequence_offsets(&self) -> &[u32] {
        &self.seq_offsets
    }

    /// The ambiguity start offsets, `num_oids + 1` entries (empty for protein).
    pub fn ambiguity_offsets(&self) -> &[u32] {
        &self.amb_offsets
    }
}
