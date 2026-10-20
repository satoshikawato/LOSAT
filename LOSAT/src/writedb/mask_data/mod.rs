//! The mask data that `makeblastdb -mask_data` stores in a database: the
//! registry of masking algorithms, the masked ranges of each sequence, and
//! the column files (`.?aa`, `.?ab`, `.?ac`) titled `BlastDb/MaskData`.
//!
//! Line NP, M4 preparation: nothing calls this yet (`LOSAT makedb` still
//! refuses `-mask_data`). The mask files themselves are read by
//! `writedb::mask_info`.
//!
//! * [`registry`]: `CMaskInfoRegistry`, the id of each algorithm.
//! * [`range_set`]: `CMaskedRangeSet`, the ranges of each sequence.
//! * [`column`]: `CWriteDB_Column`, the files of the column of a volume.
//! * [`MaskWriter`]: the part of `CWriteDB_Impl` that registers algorithms
//!   and builds the two blobs of a sequence; `WriteDb` holds one.
//!
//! A database written without mask data never creates any of these objects
//! (a `MaskWriter` holds empty containers that do not allocate, and its
//! blobs and column description stay `None`).

pub(crate) mod column;
pub mod range_set;
pub mod registry;

#[cfg(test)]
mod tests;

use column::ColumnMeta;
pub use range_set::{ranges_are_empty, MaskedRangeSet, MaskedRanges};
pub use registry::MaskInfoRegistry;

use super::volume::MaskBlobs;
use super::WriteDbError;

/// `s_EscapeColon`: a ':' of an option or a name becomes the byte 0x01.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1431-1434
// ```c++
// static const string s_EscapeColon(const string &in) {
//     const char l = 0x1;
//     return NStr::Replace(in, ":", string(l,1));
// }
// ```
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
    /// `Some` once the column exists (`m_MaskDataColumn != -1`).
    column_meta: Option<ColumnMeta>,
    /// `m_Blobs` of the mask column for the sequence being written, `Some`
    /// once `SetMaskData` has filled both (`m_HaveBlob`).
    blobs: Option<(Vec<u8>, Vec<u8>)>,
    /// The creation dates of the column index files by volume, set by a test
    /// (`set_column_date`); a volume without one reads the clock.
    dates: Vec<Option<Vec<u8>>>,
    /// The date of a column when the clock cannot be read.
    date_fallback: Option<Vec<u8>>,
}

fn unsupported_late_column() -> WriteDbError {
    WriteDbError::Other(
        "the mask column of a database is created before its first volume is written; register \
         the masking algorithms before the first sequence is published (NCBI adds blank records \
         to the open volume, which LOSAT does not)"
            .to_string(),
    )
}

impl MaskWriter {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1905-1911
    // ```c++
    // int CWriteDB_Impl::x_GetMaskDataColumnId()
    // {
    //     if (m_MaskDataColumn == -1) {
    //         m_MaskDataColumn = CreateColumn("BlastDb/MaskData", true);
    //     }
    //     return m_MaskDataColumn;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1493-1519
    // ```c++
    // int CWriteDB_Impl::CreateColumn(const string & title, bool mbo)
    // {
    //     ...
    //     m_Blobs       .push_back(new_blob);
    //     m_Blobs       .push_back(new_blob2);
    //     m_HaveBlob    .push_back(0);
    //     m_ColumnTitles.push_back(title);
    //     m_ColumnMetas .push_back(TColumnMeta());
    //
    //     if (m_Volume.NotEmpty()) {
    //         size_t id2 = m_Volume->CreateColumn(title, m_ColumnMetas.back(), mbo);
    //         ...
    //     }
    //
    //     return col_id;
    // }
    // ```
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
    // ```c++
    // int CWriteDB_Impl::
    // RegisterMaskAlgorithm(EBlast_filter_program   program,
    //                       const string          & options,
    //                       const string          & name)
    // {
    //     int algorithm_id = m_MaskAlgoRegistry.Add(program, options, name);
    //
    //     string key = NStr::IntToString(algorithm_id);
    //     string value;
    //     if (program == EBlast_filter_program::eBlast_filter_program_other) {
    //         value = NStr::IntToString((int)program) + ":" +
    //                 s_EscapeColon(options) + ":" +
    //                 s_EscapeColon(name) + ":";
    //     } else {
    //         value = NStr::IntToString((int)program) + ":" + s_EscapeColon(options);
    //     }
    //
    //     if (m_UseGiMask) {
    //         ...
    //     } else {
    //         m_ColumnMetas[x_GetMaskDataColumnId()][key] = value;
    //     }
    //
    //     return algorithm_id;
    // }
    // ```
    // LOSAT has no `-gi_mask` (it needs `-parse_seqids`).
    /// `RegisterMaskAlgorithm(program, options, name)`.
    pub(crate) fn register_program(
        &mut self,
        program: i32,
        options: &[u8],
        name: &[u8],
        volumes_exist: bool,
    ) -> Result<i32, WriteDbError> {
        if self.column_meta.is_none() && volumes_exist {
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
        self.column_meta_mut(volumes_exist)?.insert(key, value);
        Ok(algorithm_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1464-1480
    // ```c++
    // int CWriteDB_Impl::
    // RegisterMaskAlgorithm(const string &id,
    //                       const string &description,
    //                       const string &options)
    // {
    //     int algorithm_id = m_MaskAlgoRegistry.Add(id);
    //
    //     string key = NStr::IntToString(algorithm_id);
    //     string value = "100:" +
    //          s_EscapeColon(options) + ":" +
    //          s_EscapeColon(id) + ":" +
    //          s_EscapeColon(description);
    //
    //     m_ColumnMetas[x_GetMaskDataColumnId()][key] = value;
    //
    //     return algorithm_id;
    // }
    // ```
    /// `RegisterMaskAlgorithm(id, description, options)`: the `-mask_id`
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
    // ```c++
    // void CWriteDB_Impl::SetMaskData(const CMaskedRangesVector & ranges,
    //                                 const vector <TGi>        & gis)
    // {
    //     ...
    //     TSeqPos seq_length = x_ComputeSeqLength();
    //     ...
    //     int range_list_count = 0;
    //     int offset_pairs_count = 0;
    //
    //     ITERATE(CMaskedRangesVector, r1, ranges) {
    //         if (r1->empty()) {
    //             continue;
    //         }
    //
    //         range_list_count ++;
    //         offset_pairs_count += r1->offsets.size();
    //
    //         if ( !m_MaskAlgoRegistry.IsRegistered(r1->algorithm_id) ) {
    //             string msg("Error: Algorithm IDs must be registered before use.");
    //             msg += " Unknown algorithm ID = " +
    //                 NStr::IntToString((int)r1->algorithm_id);
    //             NCBI_THROW(CWriteDBException, eArgErr, msg);
    //         }
    //
    //         ITERATE(TPairVector, r2, r1->offsets) {
    //             if ((r2->first  > r2->second) ||
    //                 (r2->second > seq_length)) {
    //
    //                 NCBI_THROW(CWriteDBException,
    //                            eArgErr,
    //                            "Error: Masked data offsets out of bounds.");
    //             }
    //         }
    //     }
    //     ...
    //     // OID-based masks
    //     const int col_id = x_GetMaskDataColumnId();
    //     CBlastDbBlob & blob = SetBlobData(col_id);
    //     blob.Clear();
    //     blob.WriteInt4(range_list_count);
    //
    //     CBlastDbBlob & blob2 = SetBlobData(col_id);
    //     blob2.Clear();
    //     blob2.WriteInt4(range_list_count);
    //
    //     ITERATE(CMaskedRangesVector, r1, ranges) {
    //         if (r1->offsets.size()) {
    //             blob.WriteInt4(r1->algorithm_id);
    //             blob.WriteInt4(r1->offsets.size());
    //             blob2.WriteInt4(r1->algorithm_id);
    //             blob2.WriteInt4(r1->offsets.size());
    //
    //             ITERATE(TPairVector, r2, r1->offsets) {
    //                 blob.WriteInt4(r2->first);
    //                 blob.WriteInt4(r2->second);
    //                 blob2.WriteInt4_LE(r2->first);
    //                 blob2.WriteInt4_LE(r2->second);
    //             }
    //         }
    //     }
    //
    //     blob.WritePadBytes(4, CBlastDbBlob::eSimple);
    //     blob2.WritePadBytes(4, CBlastDbBlob::eSimple);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1537-1549
    // ```c++
    // CBlastDbBlob & CWriteDB_Impl::SetBlobData(int col_id)
    // {
    //     ...
    //     if (m_HaveBlob[col_id] > 1) {
    //         NCBI_THROW(CWriteDBException, eArgErr,
    //                    "Error: Already have blob for this sequence and column");
    //     }
    //
    //     ++m_HaveBlob[col_id];
    // ```
    // `seq_length` is `x_ComputeSeqLength()` for the sequence being written,
    // `None` when there is none ("Need sequence data."). All lengths and
    // offsets are 4-byte integers; every size is a multiple of 4, so
    // `WritePadBytes(4, eSimple)` adds nothing.
    /// `SetMaskData` (without the GI-mask path): the blobs of the sequence
    /// being written.
    pub(crate) fn set_mask_data(
        &mut self,
        ranges: &[MaskedRanges],
        seq_length: Option<u32>,
        volumes_exist: bool,
    ) -> Result<(), WriteDbError> {
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
    // ```c++
    //     NON_CONST_ITERATE(vector<int>, iter, m_HaveBlob) {
    //         *iter = 0;
    //     }
    //     ...
    //     NON_CONST_ITERATE(vector< CRef<CBlastDbBlob> >, iter, m_Blobs) {
    //         (**iter).Clear();
    //     }
    // ```
    /// `x_ResetSequenceData` for the blobs: the next sequence starts without
    /// mask data.
    pub(crate) fn reset_sequence(&mut self) {
        self.blobs = None;
    }

    /// The blobs `x_Publish` hands to `WriteSequence`: empty when the
    /// sequence has no mask data.
    pub(crate) fn blobs(&self) -> MaskBlobs<'_> {
        match &self.blobs {
            Some((blob, blob2)) => MaskBlobs { blob, blob2 },
            None => MaskBlobs {
                blob: &[],
                blob2: &[],
            },
        }
    }

    /// The metadata of the mask column, `None` when the database has none.
    pub(crate) fn column_meta(&self) -> Option<&ColumnMeta> {
        self.column_meta.as_ref()
    }

    /// Test hook: the creation date of the mask column of volume `index`.
    pub(crate) fn set_column_date(&mut self, index: usize, date: &[u8]) {
        if self.dates.len() <= index {
            self.dates.resize(index + 1, None);
        }
        self.dates[index] = Some(date.to_vec());
    }

    /// The date of the build, for a clock that cannot be read.
    pub(crate) fn set_date_fallback(&mut self, date: &[u8]) {
        self.date_fallback = Some(date.to_vec());
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_column.cpp:199-199
    // ```c++
    //     m_Date = CTime(CTime::eCurrent).AsString();
    // ```
    /// The date stored in the index file of the column of volume `index`:
    /// the one a test set, or `CTime(CTime::eCurrent).AsString()`.
    pub(crate) fn column_date(&self, index: usize) -> Result<Vec<u8>, WriteDbError> {
        if let Some(Some(date)) = self.dates.get(index) {
            return Ok(date.clone());
        }
        if let Some(now) = super::time::LocalTime::now() {
            return Ok(now.log_time().into_bytes());
        }
        self.date_fallback.clone().ok_or_else(|| {
            WriteDbError::Other(
                "the system has no clock for the date of the mask column files".to_string(),
            )
        })
    }
}
