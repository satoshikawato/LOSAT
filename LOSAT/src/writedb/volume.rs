//! One volume of a database being written (`CWriteDB_Volume`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp`).

use super::convert::find_sequence_length;
use super::files::{make_short_name, orig_num_digits, HeaderFile, IndexFile, SequenceFile};
use super::WriteDbError;

/// `CWriteDB_Volume` without ISAM indices, columns and OID masks.
#[derive(Debug)]
pub(crate) struct Volume {
    dbname: String,
    protein: bool,
    vol_name: String,
    index: u32,
    idx: IndexFile,
    hdr: HeaderFile,
    seq: SequenceFile,
    /// `m_OID`: the number of sequences written to this volume.
    oid: u32,
    open: bool,
}

impl Volume {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:65-84
    // ```c++
    //     m_VolName = CWriteDB_File::MakeShortName(m_DbName, m_Index);
    //
    //     m_Idx.Reset(new CWriteDB_IndexFile(dbname,
    //                                        protein,
    //                                        title,
    //                                        date,
    //                                        index,
    //                                        max_file_size,
    //                                        dbver));
    //
    //     m_Hdr.Reset(new CWriteDB_HeaderFile(dbname,
    //                                         protein,
    //                                         index,
    //                                         max_file_size));
    //
    //     m_Seq.Reset(new CWriteDB_SequenceFile(dbname,
    //                                           protein,
    //                                           index,
    //                                           max_file_size,
    //                                           max_letters));
    // ```
    pub(crate) fn create(
        dbname: &str,
        protein: bool,
        title: &[u8],
        date: &[u8],
        index: u32,
        max_file_size: u64,
    ) -> Result<Self, WriteDbError> {
        Ok(Self {
            dbname: dbname.to_string(),
            protein,
            vol_name: make_short_name(dbname, index),
            index,
            idx: IndexFile::create(dbname, protein, title, date, index, max_file_size)?,
            hdr: HeaderFile::create(dbname, protein, index, max_file_size)?,
            seq: SequenceFile::create(dbname, protein, index, max_file_size)?,
            oid: 0,
            open: true,
        })
    }

    /// `GetOID`: the number the next sequence gets in this volume.
    pub(crate) fn oid(&self) -> u32 {
        self.oid
    }

    /// `GetVolumeName`: the path of the volume without extension.
    pub(crate) fn volume_name(&self) -> &str {
        &self.vol_name
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:167-185
    // ```c++
    //     if (! (seq.size() && binhdr.size())) {
    //             NCBI_THROW(CWriteDBException,
    //                        eArgErr,
    //                        "Error: Cannot find CBioseq or deflines.");
    //     }
    //
    //     _ASSERT(m_Open);
    //
    //     int length = (m_Protein
    //                   ? (int) seq.size()
    //                   : x_FindNuclLength(seq));
    //
    //     bool overfull = false;
    //
    //     if (! (m_Idx->CanFit() &&
    //            m_Hdr->CanFit((int)binhdr.size()) &&
    //            m_Seq->CanFit((int)(seq.size() + ambig.size()), length))) {
    //         overfull = true;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:219-226
    // ```c++
    //     // Exception - if volume has no data, ignore the file size limits;
    //     // otherwise there would be either a hard failure or an infinite
    //     // recursion of building empty volumes.  Building a volume that's
    //     // too big is considered preferable to either of these outcomes.
    //
    //     if (m_OID && overfull) {
    //         return false;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:254-264
    // ```c++
    //     unsigned int off_hdr(0), off_seq(0), off_amb(0);
    //
    //     m_Hdr->AddSequence(binhdr, off_hdr);
    //
    //     if (m_Protein) {
    //         m_Seq->AddSequence(seq, off_seq, length);
    //         m_Idx->AddSequence((int) seq.size(), off_hdr, off_seq);
    //     } else {
    //         m_Seq->AddSequence(seq, ambig, off_seq, off_amb, length);
    //         m_Idx->AddSequence(length, off_hdr, off_seq, off_amb);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:321-323
    // ```c++
    //     m_OID ++;
    //
    //     return true;
    // ```
    /// `WriteSequence`: adds one sequence, or returns `false` when the
    /// volume already holds a sequence and this one does not fit a file.
    pub(crate) fn write_sequence(
        &mut self,
        seq: &[u8],
        ambig: &[u8],
        binhdr: &[u8],
    ) -> Result<bool, WriteDbError> {
        if seq.is_empty() || binhdr.is_empty() {
            return Err(WriteDbError::Creation(
                "Error: Cannot find CBioseq or deflines.".to_string(),
            ));
        }
        let length = find_sequence_length(self.protein, seq);
        let overfull = !(self.idx.can_fit()
            && self.hdr.can_fit(binhdr.len())
            && self.seq.can_fit(seq.len() + ambig.len()));
        if self.oid != 0 && overfull {
            return Ok(false);
        }
        let off_hdr = self.hdr.add_sequence(binhdr)?;
        if self.protein {
            let off_seq = self.seq.add_protein(seq)?;
            self.idx.add_protein(length, off_hdr, off_seq);
        } else {
            let (off_seq, off_amb) = self.seq.add_nucleotide(seq, ambig)?;
            self.idx.add_nucleotide(length, off_hdr, off_seq, off_amb);
        }
        self.oid += 1;
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:336-342
    // ```c++
    //     if (m_Open) {
    //         m_Open = false;
    //
    //         // close each file.
    //         m_Idx->Close();
    //         m_Hdr->Close();
    //         m_Seq->Close();
    // ```
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if self.open {
            self.open = false;
            self.idx.close()?;
            self.hdr.close()?;
            self.seq.close()?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:378-384
    // ```c++
    //     _ASSERT(! m_Open);
    //     m_VolName = m_DbName;
    //
    //     // rename all files to 'single volume' notation.
    //     m_Idx->RenameSingle();
    //     m_Hdr->RenameSingle();
    //     m_Seq->RenameSingle();
    // ```
    pub(crate) fn rename_single(&mut self) -> Result<(), WriteDbError> {
        self.vol_name = self.dbname.clone();
        self.idx.file_mut().rename_single()?;
        self.hdr.file_mut().rename_single()?;
        self.seq.file_mut().rename_single()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:419-428
    // ```c++
    //     _ASSERT(! m_Open);
    //     m_Idx->RenameFileIndex(num_digits);
    //     m_Hdr->RenameFileIndex(num_digits);
    //     m_Seq->RenameFileIndex(num_digits);
    //
    //     if (log10(m_Index) +1 < num_digits) {
    //     	string index_filename = m_Idx->GetFilename();
    //     	size_t t = index_filename.find_last_of(".");
    //     	m_VolName = index_filename.substr(0, t);
    //     }
    // ```
    // The comparison is in floating point: it holds for index 0
    // (`log10(0)` is minus infinity) and for every index with fewer decimal
    // digits than `num_digits`.
    pub(crate) fn rename_file_index(&mut self, num_digits: u32) -> Result<(), WriteDbError> {
        self.idx.file_mut().rename_file_index(num_digits)?;
        self.hdr.file_mut().rename_file_index(num_digits)?;
        self.seq.file_mut().rename_file_index(num_digits)?;
        if self.index == 0 || orig_num_digits(self.index) < num_digits {
            let index_filename = self.idx.file().filename();
            if let Some(t) = index_filename.rfind('.') {
                self.vol_name = index_filename[..t].to_string();
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:461-465
    // ```c++
    // void CWriteDB_Volume::ListFiles(vector<string> & files) const
    // {
    //     files.push_back(m_Idx->GetFilename());
    //     files.push_back(m_Hdr->GetFilename());
    //     files.push_back(m_Seq->GetFilename());
    // ```
    pub(crate) fn list_files(&self, files: &mut Vec<String>) {
        files.push(self.idx.file().filename().to_string());
        files.push(self.hdr.file().filename().to_string());
        files.push(self.seq.file().filename().to_string());
    }

    /// Drops the open files without writing what is still buffered, so that
    /// a failed build can remove them on every platform.
    pub(crate) fn abandon(&mut self) {
        self.open = false;
        self.idx.file_mut().abandon();
        self.hdr.file_mut().abandon();
        self.seq.file_mut().abandon();
    }
}
