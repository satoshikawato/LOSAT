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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:86-135
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1206-1210
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

    /// `GetVolumeName`: the path of the volume without extension.
    pub(crate) fn volume_name(&self) -> &str {
        &self.vol_name
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:167-185
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:186-204
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:219-226
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:228-253
    // A record without ids (`idlist` empty) adds nothing to the set and reports the empty
    // label; makeblastdb's FASTA path gives every record at least one id.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:254-264
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:321-323
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:266-291
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
        // Every sequence adds its blobs to the column, empty ones too.
        if let Some(column) = &mut self.mask_column {
            column.add_blob(mask.blob, mask.blob2)?;
        }
        self.oid += 1;
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:336-342
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:344-361
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:386-402
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
    // The comparison is in floating point: it holds for index 0
    // (`log10(0)` is minus infinity) and for every index with fewer decimal
    // digits than `num_digits`.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:430-446
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:467-490
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
