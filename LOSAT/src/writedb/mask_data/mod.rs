//! The mask data that `makeblastdb -mask_data` stores in a database: the
//! registry of masking algorithms, the masked ranges of each sequence, and
//! the column files (`.?aa`, `.?ab`, `.?ac`) titled `BlastDb/MaskData`.
//!
//! `losat makeblastdb -mask_data` (and the library's `run_makeblastdb`) reach it through
//! [`process`]. The mask files themselves are read by
//! `writedb::mask_info`; [`process`] joins the two as makeblastdb's
//! `x_ProcessMaskData` and documents the order in which the command line
//! calls everything.
//!
//! * [`process`]: `x_ProcessMaskData`, the mask files into the registry and
//!   the set of ranges.
//! * [`registry`]: `CMaskInfoRegistry`, the id of each algorithm.
//! * [`range_set`]: `CMaskedRangeSet`, the ranges of each sequence.
//! * [`MaskWriter`]: the part of `CWriteDB_Impl` that registers algorithms
//!   and builds the two blobs of a sequence, or with `-gi_mask` hands the ranges to the GI masks
//!   (`writedb::gi_mask`); `WriteDb` holds one.
//!
//! `makeblastdb -gi_mask` decides once, for the whole build, where the masks go
//! ([`MaskWriter::use_gi_masks`], NCBI's `m_UseGiMask`): without it to the mask column of each
//! volume, with it to the GI-indexed files of each algorithm, and then no mask column is made.
//!
//! Which sequence a mask reaches. A mask is stored under the `CSeq_id_Handle` of its Seq-id and a
//! sequence looks up every Seq-id of its `Bioseq::GetId()` in order, as
//! `CMaskedRangeSet::{Insert, GetRanges}` do (masked_range_set.cpp:39-106; build_db.cpp:499-502,
//! 837-851). Two Seq-ids are one handle by the rules of `CSeq_id_Mapper` (seq_id_handle.cpp:113-116,
//! seq_id_tree.cpp:244-275 and the tree of each alternative): `embl|AY100001.1|` masks
//! `gb|AY100001.1|`, accessions and local string ids match without case, `local str "5"` does not
//! match `lcl|5`; `seqid::mask_handle` has the table. With `-parse_seqids` a sequence has the
//! Seq-ids of its defline (those of the other molecule type turned into `local str`,
//! build_db.cpp:837-851), without it `local id N`. The ranges of several Seq-ids are combined as
//! `CSeq_loc::Add` with the flags `fMerge_All` and `fSort` does (Seq_loc.cpp:4283-4693), empty and
//! whole ranges included; `giim` Seq-ids find their masks in the order of the mask files, as NCBI's
//! Giim tree.
//!
//! What is refused (`WriteDbError::Unsupported`: "losat makeblastdb: ...", status 1, naming the
//! Seq-ids and the sequence; `range_set` has the rule and NCBI's lines): only what NCBI orders by
//! memory addresses, ranges of two different handles with the same `m_Packed` written on one
//! strand of a sequence in an order that matters (`CSeq_id_Handle::operator<` falls back to the
//! addresses of the infos, seq_id_handle.hpp:193-207; two `giim` handles of one number are such a
//! pair). That order shows only in the files of a build that succeeds, so the refusal comes at the
//! close of such a build (`BuildDatabase::end_build`), after every error of NCBI that the rest of
//! the run gives, and leaves the files that an error of NCBI's build there leaves. The environment
//! variables `OBJECTS_PACK_TEXTID` and `OBJECTS_PACK_GENERAL` (seq_id_tree.cpp:692-707) change the
//! handles and are refused at the start of the program (`blastinput::ncbi_environment`).
//!
//! A database written without mask data never creates any of these objects
//! (a `MaskWriter` holds empty containers that do not allocate, and its
//! blobs and column description stay `None`).

pub(crate) mod process;
pub(crate) mod range_set;
pub(crate) mod registry;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_cli;
#[cfg(test)]
mod tests_e2e;
#[cfg(test)]
mod tests_ids;
#[cfg(test)]
mod tests_ncbi_unit;

use super::column::{unsupported_late_column, BlobPair, ColumnDates, ColumnMeta};
pub(crate) use range_set::{ranges_are_empty, MaskedRangeSet, MaskedRanges};
pub(crate) use registry::MaskInfoRegistry;

use super::gi_mask::GiMasks;
use super::WriteDbError;

/// `s_EscapeColon`: a ':' of an option or a name becomes the byte 0x01.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1431-1434
fn escape_colon(text: &[u8]) -> Vec<u8> {
    text.iter()
        .map(|&b| if b == b':' { 0x01 } else { b })
        .collect()
}

/// The mask data part of `CWriteDB_Impl`.
#[derive(Debug, Default)]
pub(crate) struct MaskWriter {
    /// `m_MaskAlgoRegistry`.
    registry: MaskInfoRegistry,
    /// `m_ColumnMetas[m_MaskDataColumn]`: the metadata of the mask column,
    /// `Some` once the column exists (`m_MaskDataColumn` is not -1).
    column_meta: Option<ColumnMeta>,
    /// `m_Blobs` of the mask column for the sequence being written, `Some`
    /// once `SetMaskData` has filled both (`m_HaveBlob`).
    blobs: Option<(Vec<u8>, Vec<u8>)>,
    /// `m_UseGiMask` with `m_GiMasks`: `Some` for a build with `-gi_mask`, whose algorithms
    /// write GI-indexed files instead of the mask column.
    gi_masks: Option<GiMasks>,
}

impl MaskWriter {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1905-1911
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1493-1519
    // The mask column is created by the first registration (or, without one,
    // by the first `SetMaskData`). NCBI also creates it in a volume that is
    // already open; LOSAT does not (`volumes_exist`).
    fn column_meta_mut(&mut self, volumes_exist: bool) -> Result<&mut ColumnMeta, WriteDbError> {
        if self.column_meta.is_none() {
            if volumes_exist {
                return Err(unsupported_late_column());
            }
            self.column_meta = Some(ColumnMeta::new());
        }
        match &mut self.column_meta {
            Some(meta) => Ok(meta),
            None => Err(unsupported_late_column()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1436-1462
    // With `-gi_mask` (`m_UseGiMask`) a registration makes no entry of the mask column: the id
    // maps to a new GI mask of the algorithm (`CWriteDB_GiMask`). Its description, in the header
    // of the algorithm's index files, is the value the column would have stored, `name` the name
    // of its files, and `max_file_size` (`m_MaxFileSize`, `-max_file_sz`) the limit of its data
    // volumes. The two index files read the clock (`clock`) in their constructors, the `.gmi`
    // file first.
    /// `RegisterMaskAlgorithm` of a program, its options and a name.
    pub(crate) fn register_program(
        &mut self,
        program: i32,
        options: &[u8],
        name: &[u8],
        volumes_exist: bool,
        max_file_size: u64,
        clock: &ColumnDates,
    ) -> Result<i32, WriteDbError> {
        if self.gi_masks.is_none() && self.column_meta.is_none() && volumes_exist {
            return Err(unsupported_late_column());
        }
        let algorithm_id = self.registry.add_program(program, options, name)?;

        let key = algorithm_id.to_string().into_bytes();
        let mut value = program.to_string().into_bytes();
        value.push(b':');
        value.extend_from_slice(&escape_colon(options));
        if program == registry::FILTER_PROGRAM_OTHER {
            value.push(b':');
            value.extend_from_slice(&escape_colon(name));
            value.push(b':');
        }
        if self.gi_masks.is_some() {
            let Ok(name) = std::str::from_utf8(name) else {
                return Err(WriteDbError::Unsupported(
                    "-gi_mask: the name of a masking algorithm's GI-mask files is not valid UTF-8"
                        .to_string(),
                ));
            };
            let dates = [clock.now()?, clock.now()?];
            if let Some(gi_masks) = self.gi_masks.as_mut() {
                gi_masks.add(algorithm_id, name, &value, max_file_size, dates);
            }
            return Ok(algorithm_id);
        }
        self.column_meta_mut(volumes_exist)?.insert(key, value);
        Ok(algorithm_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1464-1480
    /// `RegisterMaskAlgorithm` of an id, a description and options: the `-mask_id`
    /// form.
    pub(crate) fn register_id(
        &mut self,
        id: &[u8],
        description: &[u8],
        options: &[u8],
        volumes_exist: bool,
    ) -> Result<i32, WriteDbError> {
        if self.column_meta.is_none() && volumes_exist {
            return Err(unsupported_late_column());
        }
        let algorithm_id = self.registry.add_id(id)?;

        let key = algorithm_id.to_string().into_bytes();
        let mut value = b"100:".to_vec();
        value.extend_from_slice(&escape_colon(options));
        value.push(b':');
        value.extend_from_slice(&escape_colon(id));
        value.push(b':');
        value.extend_from_slice(&escape_colon(description));
        self.column_meta_mut(volumes_exist)?.insert(key, value);
        Ok(algorithm_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1339-1429
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1537-1549
    // `seq_length` is `x_ComputeSeqLength()` for the sequence being written,
    // `None` when there is none ("Need sequence data."). All lengths and
    // offsets are 4-byte integers; every size is a multiple of 4, so
    // `WritePadBytes` to a multiple of 4 bytes with `eSimple` adds nothing. With `-gi_mask` the
    // masks of a sequence without a GI (`gis` empty) are dropped before any check (survey D4),
    // and the checked ranges of the others go to the GI masks of their algorithms.
    /// `SetMaskData`: the blobs of the sequence being written, or its ranges under its GIs `gis`
    /// with `-gi_mask`.
    pub(crate) fn set_mask_data(
        &mut self,
        ranges: &[MaskedRanges],
        gis: &[i64],
        seq_length: Option<u32>,
        volumes_exist: bool,
    ) -> Result<(), WriteDbError> {
        if self.gi_masks.is_some() && gis.is_empty() {
            return Ok(());
        }
        let Some(seq_length) = seq_length else {
            return Err(WriteDbError::Creation("Need sequence data.".to_string()));
        };
        let mut range_list_count: u32 = 0;
        for r1 in ranges {
            if r1.offsets.is_empty() {
                continue;
            }
            range_list_count += 1;
            if !self.registry.is_registered(r1.algorithm_id) {
                return Err(WriteDbError::Creation(format!(
                    "Error: Algorithm IDs must be registered before use. Unknown algorithm ID = {}",
                    r1.algorithm_id
                )));
            }
            for &(first, second) in &r1.offsets {
                if first > second || second > seq_length {
                    return Err(WriteDbError::Creation(
                        "Error: Masked data offsets out of bounds.".to_string(),
                    ));
                }
            }
        }

        if let Some(gi_masks) = self.gi_masks.as_mut() {
            for r1 in ranges.iter().filter(|r1| !r1.offsets.is_empty()) {
                gi_masks.add_gi_mask(r1.algorithm_id, gis, &r1.offsets)?;
            }
            return Ok(());
        }
        self.column_meta_mut(volumes_exist)?;
        if self.blobs.is_some() {
            return Err(WriteDbError::Creation(
                "Error: Already have blob for this sequence and column".to_string(),
            ));
        }

        let pairs: usize = ranges.iter().map(|r| r.offsets.len()).sum();
        let capacity = 4 + 8 * ranges.len() + 8 * pairs;
        let mut blob = Vec::with_capacity(capacity);
        let mut blob2 = Vec::with_capacity(capacity);
        blob.extend_from_slice(&range_list_count.to_be_bytes());
        blob2.extend_from_slice(&range_list_count.to_be_bytes());
        for r1 in ranges {
            if r1.offsets.is_empty() {
                continue;
            }
            let algorithm_id = r1.algorithm_id.to_be_bytes();
            let count = (r1.offsets.len() as u32).to_be_bytes();
            blob.extend_from_slice(&algorithm_id);
            blob.extend_from_slice(&count);
            blob2.extend_from_slice(&algorithm_id);
            blob2.extend_from_slice(&count);
            for &(first, second) in &r1.offsets {
                blob.extend_from_slice(&first.to_be_bytes());
                blob.extend_from_slice(&second.to_be_bytes());
                blob2.extend_from_slice(&first.to_le_bytes());
                blob2.extend_from_slice(&second.to_le_bytes());
            }
        }
        self.blobs = Some((blob, blob2));
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:128-135
    /// `x_ResetSequenceData` for the blobs: the next sequence starts without
    /// mask data.
    pub(crate) fn reset_sequence(&mut self) {
        self.blobs = None;
    }

    /// The blobs `x_Publish` hands to `WriteSequence`: empty when the
    /// sequence has no mask data.
    pub(crate) fn blobs(&self) -> BlobPair<'_> {
        match &self.blobs {
            Some((blob, blob2)) => BlobPair { blob, blob2 },
            None => BlobPair::EMPTY,
        }
    }

    /// The metadata of the mask column, `None` when the database has none.
    pub(crate) fn column_meta(&self) -> Option<&ColumnMeta> {
        self.column_meta.as_ref()
    }

    /// The metadata of the mask column to change (`AddColumnMetaData`), `None` when the
    /// database has none.
    pub(crate) fn existing_column_meta_mut(&mut self) -> Option<&mut ColumnMeta> {
        self.column_meta.as_mut()
    }

    /// `m_UseGiMask` set (the `use_gi_mask` argument of `CWriteDB`): the algorithms registered
    /// from now on write GI-indexed files.
    pub(crate) fn use_gi_masks(&mut self) {
        self.gi_masks = Some(GiMasks::default());
    }

    /// `m_UseGiMask`.
    pub(crate) fn uses_gi_masks(&self) -> bool {
        self.gi_masks.is_some()
    }

    /// The close of the GI masks (`CWriteDB_Impl::Close`, after the last volume's).
    pub(crate) fn close_gi_masks(&mut self) {
        if let Some(gi_masks) = self.gi_masks.as_mut() {
            gi_masks.close();
        }
    }

    /// The names of the GI masks after their close (`GetName`, empty for one without GIs), in
    /// the order of registration; none without `-gi_mask`.
    pub(crate) fn gi_mask_names(&self) -> impl Iterator<Item = &str> {
        self.gi_masks.iter().flat_map(GiMasks::names)
    }
}
