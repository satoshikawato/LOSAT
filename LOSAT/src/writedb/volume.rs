//! One volume of a database being written (`CWriteDB_Volume`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp`).

use super::convert::find_sequence_length;
use super::files::{make_short_name, orig_num_digits, HeaderFile, IndexFile, SequenceFile};
use super::mask_data::column::{Column, ColumnMeta};
use super::WriteDbError;

/// `CWriteDB_Volume` without ISAM indices and OID masks, and with the mask
/// column as its only column.
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
    /// `m_Columns`: the mask column (`BlastDb/MaskData`) when the database
    /// has mask data.
    mask_column: Option<Column>,
}

/// The two blobs of the mask column for one sequence (the offsets in
/// big-endian and in little-endian byte order); empty for a sequence
/// without mask data.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MaskBlobs<'a> {
    pub(crate) blob: &'a [u8],
    pub(crate) blob2: &'a [u8],
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
        // NCBI does not test whether the three files opened
        // (`WriteDbFile::create`). LOSAT stops at the first that cannot be
        // created, and this volume is then not one of the build's volumes,
        // whose files the failed build removes; so the files created here
        // before the failure are removed here.
        let mut idx = IndexFile::create(dbname, protein, title, date, index, max_file_size)?;
        let mut hdr = match HeaderFile::create(dbname, protein, index, max_file_size) {
            Ok(hdr) => hdr,
            Err(error) => {
                idx.file_mut().discard();
                return Err(error);
            }
        };
        let seq = match SequenceFile::create(dbname, protein, index, max_file_size) {
            Ok(seq) => seq,
            Err(error) => {
                idx.file_mut().discard();
                hdr.file_mut().discard();
                return Err(error);
            }
        };
        Ok(Self {
            dbname: dbname.to_string(),
            protein,
            vol_name: make_short_name(dbname, index),
            index,
            idx,
            hdr,
            seq,
            oid: 0,
            open: true,
            mask_column: None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:505-554
    // ```c++
    // int CWriteDB_Volume::CreateColumn(const string      & title,
    //                                   const TColumnMeta & meta,
    //                                   Uint8               max_sz,
    //                                   bool                mbo)
    // {
    //     int col_id = m_Columns.size();
    //     ...
    //     CRef<CWriteDB_Column> new_col
    //         (new CWriteDB_Column(m_DbName,
    //                              extn,
    //                              extn2,
    //                              m_Index,
    //                              title,
    //                              meta,
    //                              max_sz));
    //
    //     /* For support of multiple byte orders */
    //     if (mbo) new_col->AddByteOrder(m_DbName,
    //                              extn3,
    //                              m_Index,
    //                              max_sz);
    //
    //     // If the OID is not zero, then add all the blank records for the
    //     // prior OIDs to the new column.
    //
    //     CBlastDbBlob blank;
    //
    //     for(int j = 0; j < m_OID; j++) {
    //         if (mbo) new_col->AddBlob(blank, blank);
    //         else     new_col->AddBlob(blank);
    //     }
    //
    //     m_Columns.push_back(new_col);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1206-1210
    // ```c++
    //             for(size_t i = 0; i < m_ColumnTitles.size(); i++) {
    //                 m_Volume->CreateColumn(m_ColumnTitles[i],
    //                                        m_ColumnMetas[i],
    //                                        m_MaxFileSize);
    //             }
    // ```
    // The mask column is the volume's only column (column 0, letter `a`) and
    // `CWriteDB_Volume::CreateColumn`'s `mbo` is true by default. `date` is
    // the creation date of the column's index file.
    pub(crate) fn create_mask_column(
        &mut self,
        meta: &ColumnMeta,
        date: &[u8],
        max_file_size: u64,
    ) -> Result<(), WriteDbError> {
        let mut column = Column::new(
            &self.dbname,
            self.protein,
            self.index,
            meta,
            date,
            max_file_size,
        );
        for _ in 0..self.oid {
            column.add_blob(&[], &[])?;
        }
        self.mask_column = Some(column);
        Ok(())
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
        mask: MaskBlobs<'_>,
    ) -> Result<bool, WriteDbError> {
        if seq.is_empty() || binhdr.is_empty() {
            return Err(WriteDbError::Creation(
                "Error: Cannot find CBioseq or deflines.".to_string(),
            ));
        }
        let length = find_sequence_length(self.protein, seq);
        let mut overfull = !(self.idx.can_fit()
            && self.hdr.can_fit(binhdr.len())
            && self.seq.can_fit(seq.len() + ambig.len()));
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:209-216
        // ```c++
        //     for(int blob_i = 0; blob_i < (int) blobs.size(); blob_i++) {
        //         _ASSERT(blob_i / 2 < (int) m_Columns.size());
        //
        //         if (! m_Columns[blob_i / 2]->CanFit(blobs[blob_i]->Size())) {
        //             overfull = true;
        //             break;
        //         }
        //     }
        // ```
        // Both blobs of the mask column are checked, the first one first.
        if let Some(column) = &self.mask_column {
            if !(column.can_fit(mask.blob.len()) && column.can_fit(mask.blob2.len())) {
                overfull = true;
            }
        }
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
        // Every sequence adds its blobs to the column, empty ones too.
        if let Some(column) = &mut self.mask_column {
            column.add_blob(mask.blob, mask.blob2)?;
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
        // The columns are closed after the other files, each time `Close`
        // is called (NCBI: writedb_volume.cpp:362-367).
        if let Some(column) = &mut self.mask_column {
            column.close()?;
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
        self.seq.file_mut().rename_single()?;
        if let Some(column) = &mut self.mask_column {
            column.rename_single()?;
        }
        Ok(())
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
        if let Some(column) = &mut self.mask_column {
            column.rename_file_index(num_digits)?;
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
        if let Some(column) = &self.mask_column {
            column.list_files(files);
        }
    }

    /// Drops the open files without writing what is still buffered, so that
    /// a failed build can remove them on every platform.
    pub(crate) fn abandon(&mut self) {
        self.open = false;
        self.idx.file_mut().abandon();
        self.hdr.file_mut().abandon();
        self.seq.file_mut().abandon();
        if let Some(column) = &mut self.mask_column {
            column.abandon();
        }
    }
}
