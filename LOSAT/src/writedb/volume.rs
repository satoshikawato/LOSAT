//! One volume of a database being written (`CWriteDB_Volume`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp`).

use super::convert::find_sequence_length;
use super::files::{
    make_short_name, orig_num_digits, DbVersion, HeaderFile, IndexFile, OffsetArrays, SequenceFile,
};
use super::id_set::LabelSet;
use super::isam::{GiIndex, IsamType, NumericIsam, PigIsam, StringIsam, StringIsamScratch};
use super::mask_data::column::{Column, ColumnMeta};
use super::WriteDbError;
use crate::seqid::{label, ParsedSeqId};

/// `CWriteDB::EIndexType` as makeblastdb sets it (makeblastdb.cpp:1207-1209): `eFullIndex`
/// for `-parse_seqids`, `eAddHash` for `-hash_index`; `eNoIndex` without both.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Indices {
    /// `eFullIndex` (`-parse_seqids`).
    pub full: bool,
    /// `eAddHash` (`-hash_index`).
    pub hash: bool,
}

impl Indices {
    /// `m_Indices != CWriteDB::eNoIndex`.
    pub fn any(self) -> bool {
        self.full || self.hash
    }
}

/// The memory of a volume's indices that the next volume takes over (as `OffsetArrays`):
/// the duplicate check's set and the record chunks of the string indices, emptied.
#[derive(Debug, Default)]
pub(crate) struct IndexScratch {
    id_set: LabelSet,
    acc: StringIsamScratch,
    hash: StringIsamScratch,
}

/// The ISAM indices, the OID to GI file and the duplicate check of a volume built with
/// indices (`m_Indices != eNoIndex`).
#[derive(Debug)]
struct VolumeIndices {
    /// `m_PigIsam` (protein only).
    pig: Option<PigIsam>,
    /// `m_GiIsam`.
    gi: NumericIsam,
    /// `m_AccIsam` (version 4 only).
    acc: Option<StringIsam>,
    /// `m_HashIsam` (`eAddHash`).
    hash: Option<StringIsam>,
    /// `m_GiIndex`.
    gi_index: GiIndex,
    /// `m_IdSet`.
    id_set: LabelSet,
    /// The label being built.
    label: Vec<u8>,
}

/// `CWriteDB_Volume` without OID masks, and with the mask column as its only column.
#[derive(Debug)]
pub(crate) struct Volume {
    dbname: String,
    protein: bool,
    vol_name: String,
    index: u32,
    idx: IndexFile,
    hdr: HeaderFile,
    seq: SequenceFile,
    /// The indices of `m_Indices != eNoIndex`.
    indices: Option<VolumeIndices>,
    /// `m_OID`: the number of sequences written to this volume.
    oid: u32,
    open: bool,
    /// `m_Columns`: the mask column (`BlastDb/MaskData`) when the database has mask data.
    mask_column: Option<Column>,
}

/// The two blobs of the mask column for one sequence (the offsets in big-endian and in
/// little-endian byte order); empty for a sequence without mask data.
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:86-135
    // ```c++
    //     if (m_Indices != CWriteDB::eNoIndex) {
    //         bool sparse =
    //             (m_Indices & CWriteDB::eSparseIndex) == CWriteDB::eSparseIndex;
    //
    //         if (m_Protein) {
    //             m_PigIsam.Reset(new CWriteDB_Isam(ePig,
    //                                               dbname,
    //                                               protein,
    //                                               index,
    //                                               max_file_size,
    //                                               false));
    //         }
    //
    //         m_GiIsam.Reset(new CWriteDB_Isam(eGi,
    //                                          dbname,
    //                                          protein,
    //                                          index,
    //                                          max_file_size,
    //                                          false));
    //         if(m_DbVersion != eBDB_Version5) {
    //             m_AccIsam.Reset(new CWriteDB_Isam(eAcc,
    //                                           dbname,
    //                                           protein,
    //                                           index,
    //                                           max_file_size,
    //                                           sparse));
    //         }
    //         if (m_Indices & CWriteDB::eAddTrace) {
    //             ...
    //         }
    //
    //         if (m_Indices & CWriteDB::eAddHash) {
    //             m_HashIsam.Reset(new CWriteDB_Isam(eHash,
    //                                                dbname,
    //                                                protein,
    //                                                index,
    //                                                max_file_size,
    //                                                false));
    //         }
    //
    //         m_GiIndex.Reset(new CWriteDB_GiIndex(dbname,
    //                                              protein,
    //                                              index,
    //                                              max_file_size));
    //     }
    // ```
    // makeblastdb never sets `eAddTrace` or `eSparseIndex`.
    /// `arrays`: the offset arrays of the previous volume
    /// (`take_offset_arrays`), or new ones for the first; `sizeof_long` the
    /// platform's `sizeof(long)` (`Platform::sizeof_long`); `indices` makeblastdb's index
    /// type, `scratch` the index memory of the previous volume (`take_index_scratch`).
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
        indices: Indices,
        scratch: IndexScratch,
    ) -> Result<Self, WriteDbError> {
        // NCBI does not test whether the three files opened: a file that cannot be opened is
        // written to nowhere (`super::ofstream`) and the build goes on.
        let idx = IndexFile::create(
            dbname,
            protein,
            title,
            date,
            index,
            max_file_size,
            version,
            sizeof_long,
            arrays,
        );
        let hdr = HeaderFile::create(dbname, protein, index, max_file_size);
        let seq = SequenceFile::create(dbname, protein, index, max_file_size)?;
        let indices = indices.any().then(|| {
            let IndexScratch {
                mut id_set,
                acc,
                hash,
            } = scratch;
            id_set.clear();
            VolumeIndices {
                pig: protein.then(|| PigIsam::new(dbname, index, max_file_size)),
                gi: NumericIsam::new(IsamType::Gi, dbname, protein, index, max_file_size),
                acc: (version != DbVersion::V5).then(|| {
                    StringIsam::new(IsamType::Acc, dbname, protein, index, max_file_size, acc)
                }),
                hash: indices.hash.then(|| {
                    StringIsam::new(IsamType::Hash, dbname, protein, index, max_file_size, hash)
                }),
                gi_index: GiIndex::new(dbname, protein, index, max_file_size),
                id_set,
                label: Vec::new(),
            }
        });
        Ok(Self {
            dbname: dbname.to_string(),
            protein,
            vol_name: make_short_name(dbname, index),
            index,
            idx,
            hdr,
            seq,
            indices,
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
    // `CWriteDB_Volume::CreateColumn`'s `mbo` is true by default. `date` is the creation
    // date of the column's index file. The column files are the same in version 4 and
    // version 5: the LMDB files of version 5 belong to the database, not to a volume.
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

    /// The residues or bases in this volume (`m_Letters` of its index file).
    pub(crate) fn letters(&self) -> u64 {
        self.idx.letters()
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:186-204
    // ```c++
    //     if (m_Indices != CWriteDB::eNoIndex) {
    //
    //         int num = (int)idlist.size();
    //
    //
    //         if (! ( (m_AccIsam.Empty() || m_AccIsam->CanFit(num)) &&
    //                m_GiIsam->CanFit(num) &&
    //                (m_TraceIsam.Empty() || m_TraceIsam->CanFit(num)))) {
    //             overfull = true;
    //         }
    //
    //         if (m_Protein && (! m_PigIsam->CanFit(1))) {
    //             overfull = true;
    //         }
    //
    //         if (m_HashIsam.NotEmpty() && (! m_HashIsam->CanFit(1))) {
    //             overfull = true;
    //         }
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:228-253
    // ```c++
    //     // check the uniqueness of id
    //     if (m_Indices != CWriteDB::eNoIndex) {
    //     	set<string>::size_type orig_size = m_IdSet.size();
    //     	string id;
    //     	pair<set<string>::iterator, bool > rv;
    //         CSeq_id::TLabelFlags label_flags =
    //             CSeq_id::fLabel_Default | CSeq_id::fLabel_UpperCase;
    //         ITERATE(TIdList, iter, idlist) {
    //             id = kEmptyStr;
    //             (*iter)->GetLabel(&id, CSeq_id::eDefault, label_flags);
    //             rv = m_IdSet.insert(id);
    //             if((rv.second == false) && (!(*iter)->IsLocal())) {
    //             	CNcbiOstrstream msg;
    //             	msg << "Error: Duplicate seq_ids are found: " << endl << id << endl;
    //             	NCBI_THROW(CWriteDBException, eArgErr, CNcbiOstrstreamToString(msg));
    //             }
    //         }
    //
    //         if(m_IdSet.size() == orig_size) {
    //         	CNcbiOstrstream msg;
    //         	msg << "Error: Duplicate seq_ids are found: " << endl
    //     	    << id << endl;
    //         	NCBI_THROW(CWriteDBException, eArgErr, CNcbiOstrstreamToString(msg));
    //         }
    //     }
    // ```
    // A record without ids (`idlist` empty) adds nothing to the set and reports the empty
    // label; makeblastdb's FASTA path gives every record at least one id.
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:266-291
    // ```c++
    //     if (m_Indices != CWriteDB::eNoIndex) {
    //         if(m_AccIsam.NotEmpty()) m_AccIsam->AddIds(m_OID, idlist);
    //         m_GiIsam->AddIds(m_OID, idlist);
    //
    //         TGi gi = INVALID_GI;
    //         ITERATE(TIdList, iter, idlist) {
    //             const CSeq_id & seqid = **iter;
    //             if (seqid.IsGi()) {
    //                 gi = seqid.GetGi();
    //                 break;
    //             }
    //         }
    //         m_GiIndex->AddGi(gi);
    //
    //         if (m_Protein && pig) {
    //             m_PigIsam->AddPig(m_OID, pig);
    //         }
    //
    //         if (m_TraceIsam.NotEmpty()) {
    //             m_TraceIsam->AddIds(m_OID, idlist);
    //         }
    //
    //         if (m_HashIsam.NotEmpty()) {
    //             m_HashIsam->AddHash(m_OID, hash);
    //         }
    //     }
    // ```
    // The PIG of FASTA input is 0, so `m_PigIsam` stays empty. `INVALID_GI` is -1.
    /// `WriteSequence`: adds one sequence, or returns `false` when the
    /// volume already holds a sequence and this one does not fit a file. `ids` is the
    /// sequence's id list (`m_Ids`) and `hash` its sequence hash (`m_Hash`), read when the
    /// volume has indices.
    pub(crate) fn write_sequence(
        &mut self,
        seq: &[u8],
        ambig: &[u8],
        binhdr: &[u8],
        ids: &[ParsedSeqId],
        hash: u32,
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
        if let Some(indices) = &self.indices {
            let num = ids.len();
            if !(indices.acc.as_ref().is_none_or(|acc| acc.can_fit(num)) && indices.gi.can_fit(num))
            {
                overfull = true;
            }
            if indices.pig.as_ref().is_some_and(|pig| !pig.can_fit()) {
                overfull = true;
            }
            if indices.hash.as_ref().is_some_and(|hash| !hash.can_fit(1)) {
                overfull = true;
            }
        }
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
        if let Some(indices) = &mut self.indices {
            let orig_size = indices.id_set.len();
            for id in ids {
                indices.label.clear();
                label::both_upper_version(id, &mut indices.label);
                let inserted = indices.id_set.insert(&indices.label);
                if !inserted && !matches!(id, ParsedSeqId::Local(_)) {
                    return Err(duplicate_error(&indices.label));
                }
            }
            if indices.id_set.len() == orig_size {
                if ids.is_empty() {
                    indices.label.clear();
                }
                return Err(duplicate_error(&indices.label));
            }
        }
        let off_hdr = self.hdr.add_sequence(binhdr)?;
        if self.protein {
            let off_seq = self.seq.add_protein(seq)?;
            self.idx.add_protein(length, off_hdr, off_seq);
        } else {
            let (off_seq, off_amb) = self.seq.add_nucleotide(seq, ambig)?;
            self.idx.add_nucleotide(length, off_hdr, off_seq, off_amb);
        }
        if let Some(indices) = &mut self.indices {
            if let Some(acc) = &mut indices.acc {
                acc.add_ids(self.oid, ids)?;
            }
            indices.gi.add_gis(self.oid, ids);
            indices
                .gi_index
                .add_gi(ids.iter().find_map(ParsedSeqId::gi).unwrap_or(-1));
            if let Some(hash_isam) = &mut indices.hash {
                hash_isam.add_hash(self.oid, hash)?;
            }
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:344-361
    // ```c++
    //         if (m_Indices != CWriteDB::eNoIndex) {
    //             if (m_Protein) {
    //                 m_PigIsam->Close();
    //             }
    //             m_GiIsam->Close();
    //             if(m_AccIsam.NotEmpty()) m_AccIsam->Close();
    //             m_GiIndex->Close();
    //
    //             if (m_TraceIsam.NotEmpty()) {
    //                 m_TraceIsam->Close();
    //             }
    //
    //             if (m_HashIsam.NotEmpty()) {
    //                 m_HashIsam->Close();
    //             }
    //             m_IdSet.clear();
    //         }
    // ```
    // The empty PIG index writes nothing.
    pub(crate) fn close(&mut self) -> Result<(), WriteDbError> {
        if self.open {
            self.open = false;
            self.idx.close()?;
            self.hdr.close();
            self.seq.close();
            if let Some(indices) = &mut self.indices {
                indices.gi.close()?;
                if let Some(acc) = &mut indices.acc {
                    acc.close()?;
                }
                indices.gi_index.close()?;
                if let Some(hash) = &mut indices.hash {
                    hash.close()?;
                }
                indices.id_set.clear();
            }
            // The columns are closed after the other files (writedb_volume.cpp:362-367).
            if let Some(column) = &mut self.mask_column {
                column.close()?;
            }
        }
        Ok(())
    }

    /// The emptied index memory, after `close`, for the next volume.
    pub(crate) fn take_index_scratch(&mut self) -> IndexScratch {
        match &mut self.indices {
            Some(indices) => IndexScratch {
                id_set: std::mem::take(&mut indices.id_set),
                acc: indices
                    .acc
                    .as_mut()
                    .map(StringIsam::take_scratch)
                    .unwrap_or_default(),
                hash: indices
                    .hash
                    .as_mut()
                    .map(StringIsam::take_scratch)
                    .unwrap_or_default(),
            },
            None => IndexScratch::default(),
        }
    }

    /// The emptied offset arrays of the index file, after `close`.
    pub(crate) fn take_offset_arrays(&mut self) -> OffsetArrays {
        self.idx.take_offset_arrays()
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:386-402
    // ```c++
    //     if (m_Indices != CWriteDB::eNoIndex) {
    //         if (m_Protein) {
    //             m_PigIsam->RenameSingle();
    //         }
    //         m_GiIsam->RenameSingle();
    //         if(m_AccIsam.NotEmpty()) m_AccIsam->RenameSingle();
    //         m_GiIndex->RenameSingle();
    //
    //         if (m_TraceIsam.NotEmpty()) {
    //             m_TraceIsam->RenameSingle();
    //         }
    //
    //         if (m_HashIsam.NotEmpty()) {
    //             m_HashIsam->RenameSingle();
    //         }
    //     }
    // ```
    // The files of the PIG index are never created; a file of their names that is there
    // is renamed all the same.
    pub(crate) fn rename_single(&mut self) {
        self.vol_name = self.dbname.clone();
        self.idx.file_mut().rename_single();
        self.hdr.file_mut().rename_single();
        self.seq.file_mut().rename_single();
        if let Some(indices) = &mut self.indices {
            if let Some(pig) = &mut indices.pig {
                pig.rename_single();
            }
            indices.gi.rename_single();
            if let Some(acc) = &mut indices.acc {
                acc.rename_single();
            }
            indices.gi_index.file_mut().rename_single();
            if let Some(hash) = &mut indices.hash {
                hash.rename_single();
            }
        }
        if let Some(column) = &mut self.mask_column {
            column.rename_single();
        }
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:430-446
    // ```c++
    //     if (m_Indices != CWriteDB::eNoIndex) {
    //         if (m_Protein) {
    //             m_PigIsam->RenameFileIndex(num_digits);
    //         }
    //         m_GiIsam->RenameFileIndex(num_digits);
    //         if(m_AccIsam.NotEmpty()) m_AccIsam->RenameFileIndex(num_digits);
    //         m_GiIndex->RenameFileIndex(num_digits);
    //
    //         if (m_TraceIsam.NotEmpty()) {
    //             m_TraceIsam->RenameFileIndex(num_digits);
    //         }
    //
    //         if (m_HashIsam.NotEmpty()) {
    //             m_HashIsam->RenameFileIndex(num_digits);
    //         }
    //     }
    // ```
    pub(crate) fn rename_file_index(&mut self, num_digits: u32) {
        self.idx.file_mut().rename_file_index(num_digits);
        self.hdr.file_mut().rename_file_index(num_digits);
        self.seq.file_mut().rename_file_index(num_digits);
        if self.index == 0 || orig_num_digits(self.index) < num_digits {
            let index_filename = self.idx.file().filename();
            if let Some(t) = index_filename.rfind('.') {
                self.vol_name = index_filename[..t].to_string();
            }
        }
        if let Some(indices) = &mut self.indices {
            if let Some(pig) = &mut indices.pig {
                pig.rename_file_index(num_digits);
            }
            indices.gi.rename_file_index(num_digits);
            if let Some(acc) = &mut indices.acc {
                acc.rename_file_index(num_digits);
            }
            indices.gi_index.file_mut().rename_file_index(num_digits);
            if let Some(hash) = &mut indices.hash {
                hash.rename_file_index(num_digits);
            }
        }
        if let Some(column) = &mut self.mask_column {
            column.rename_file_index(num_digits);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:461-465
    // ```c++
    // void CWriteDB_Volume::ListFiles(vector<string> & files) const
    // {
    //     files.push_back(m_Idx->GetFilename());
    //     files.push_back(m_Hdr->GetFilename());
    //     files.push_back(m_Seq->GetFilename());
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:467-490
    // ```c++
    //     if (m_AccIsam.NotEmpty()) {
    //         m_AccIsam->ListFiles(files);
    //     }
    //
    //     if (m_GiIsam.NotEmpty()) {
    //         m_GiIsam->ListFiles(files);
    //     }
    //
    //     if (m_PigIsam.NotEmpty()) {
    //         m_PigIsam->ListFiles(files);
    //     }
    //
    //     if (m_TraceIsam.NotEmpty()) {
    //         m_TraceIsam->ListFiles(files);
    //     }
    //
    //     if (m_HashIsam.NotEmpty()) {
    //         m_HashIsam->ListFiles(files);
    //     }
    //
    //     if (m_GiIndex.NotEmpty()) {
    //         files.push_back(m_GiIndex->GetFilename());
    //     }
    // ```
    pub(crate) fn list_files(&self, files: &mut Vec<String>) {
        files.push(self.idx.file().filename().to_string());
        files.push(self.hdr.file().filename().to_string());
        files.push(self.seq.file().filename().to_string());
        if let Some(indices) = &self.indices {
            if let Some(acc) = &indices.acc {
                acc.list_files(files);
            }
            indices.gi.list_files(files);
            if let Some(hash) = &indices.hash {
                hash.list_files(files);
            }
            files.push(indices.gi_index.filename().to_string());
        }
        if let Some(column) = &self.mask_column {
            column.list_files(files);
        }
    }
}

/// The error of the duplicate check: "Error: Duplicate seq_ids are found: " and the label on
/// a line of its own (its bytes as they are: a label may hold bytes of 0x80 or more).
fn duplicate_error(label: &[u8]) -> WriteDbError {
    let mut message = b"Error: Duplicate seq_ids are found: \n".to_vec();
    message.extend_from_slice(label);
    message.push(b'\n');
    match String::from_utf8(message) {
        Ok(text) => WriteDbError::Creation(text),
        Err(error) => WriteDbError::CreationBytes(error.into_bytes()),
    }
}
