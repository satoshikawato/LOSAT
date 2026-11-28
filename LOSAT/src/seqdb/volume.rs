//! One database volume: its index, sequence and header files, and blocks of
//! sequences read by position.

use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use super::column::{self, Column, MASK_DATA_TITLE};
use super::file::DbFile;
use super::filter::{DeflineList, DeflineRule};
use super::header::{
    decode_defline_taxids, decode_deflines_as, decode_full, decode_full_in, holds_a_seqid, Defline,
};
use super::index::IndexFile;
use super::{SeqBlockResult, SeqDbError, SeqType, XDbIdxRead};
use crate::blastinput::ncbi_environment::ReadVerification;

/// One volume of a database (`NAME.pin` / `.psq` / `.phr`, or the `n` files).
#[derive(Debug)]
pub struct Volume {
    base: PathBuf,
    index: IndexFile,
    seq: OnceLock<DbFile>,
    hdr: OnceLock<DbFile>,
    /// The definition lines a taxid or seqid list keeps (`SeqDb::restrict_taxids`,
    /// `restrict_seqids`), none without one.
    defline_rule: Option<LineRule>,
    /// Whether the alias files' `OID_MASK_TYPE` of the volume erases the definition lines of
    /// model sequences after `defline_rule` (the `fExcludeModel` bit of `m_OidMaskType`,
    /// `oid_mask::excludes_models`; `SeqDb::read_alias_tree`).
    exclude_model: bool,
    /// The id lists of the alias files that reach the volume (`m_VolumeGiLists`,
    /// `alias::AliasOids::lists`), none without one.
    alias_lines: Option<AliasLines>,
    /// The volume's `OID_MASK_TYPE` itself, 0 without one (`m_OidMaskType`): any type makes the
    /// header reads go through NCBI's filter, which needs a set (`Volume::has_oid_mask_type`).
    oid_mask_type: u32,
    /// The volume's `MEMB_BIT`, 0 without one (`m_MemBit`): a non-zero bit filters definition
    /// lines by the first entry of `Blast-def-line.memberships` before the id lists and
    /// `OID_MASK_TYPE` (`Volume::has_memb_bit`).
    memb_bit: u32,
    /// The column files (`.?xa`...), opened on first use (`m_HaveColumns`).
    columns: OnceLock<Vec<Column>>,
    /// The failure of NCBI's OID list that a search meets at its first read of a sequence
    /// (`Volume::fail_sequence_reads`), none when the list can be built.
    read_failure: Option<super::oid_list::OidListError>,
    /// The run's settings of NCBI's checks of what it reads (`SeqDb::refuse_unchecked_reads`),
    /// none for NCBI's default.
    read_verification: Option<Arc<ReadVerification>>,
}

/// A header record in LOSAT's refusal of a setting that turns a check of NCBI's reader off
/// (`ReadVerification::refused`).
const HEADER_RECORD: &str = "a database header record";

/// `error` of a decode of a header record under `setting` (`SeqDb::refuse_unchecked_reads`):
/// NCBI's exception of a check that the setting turns off is LOSAT's refusal, any other error
/// stays. The one place of that choice for every read of a volume's records
/// (`Volume::with_header_record`, `Volume::visit_headers`).
#[cold]
#[inline(never)]
fn refused_read(setting: &ReadVerification, error: SeqDbError) -> SeqDbError {
    match &error {
        SeqDbError::Ncbi(exception) => match setting.refused(exception, HEADER_RECORD) {
            Some(refused) => SeqDbError::Setting(refused),
            None => error,
        },
        _ => error,
    }
}

/// A list's rule for the definition lines of a volume's sequences, with the volume's first OID in
/// the database (`m_VolStart`), which NCBI adds to the tag of each `BL_ORD_ID` before the rule
/// reads the Seq-ids (`Volume::line_mask`).
#[derive(Debug)]
struct LineRule {
    rule: std::sync::Arc<DeflineRule>,
    volume_start: u32,
}

/// The alias files' id lists attached to a volume, with the volume's first OID in the database,
/// which NCBI adds to the tag of each `BL_ORD_ID` before the lists read the Seq-ids.
#[derive(Debug)]
struct AliasLines {
    lists: Vec<std::sync::Arc<super::id_list::AliasIdList>>,
    volume_start: u32,
}

/// File name of `base` with the three-letter extension `{p|n}{ext}`.
pub(super) fn with_extension(base: &Path, seqtype: SeqType, ext: &str) -> PathBuf {
    let mut name = base.as_os_str().to_os_string();
    name.push(".");
    name.push(seqtype.letter());
    name.push(ext);
    PathBuf::from(name)
}

/// The header of a sequence of a volume, with its ordinal ids still relative to the volume
/// (`Volume::header`, `SeqDb::read_header`).
#[derive(Debug)]
pub enum VolumeHeader {
    /// A record whose only Seq-ids are `gnl|BL_ORD_ID|n` of the sequence's number n in its volume
    /// (`Volume::deflines`).
    Ordinal(Vec<Defline>),
    /// A record with parsed Seq-ids (`Volume::defline_set`).
    Parsed(crate::seqid::DeflineSet),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1804
/// The signed 32-bit length that NCBI's `CSeqDBVol::GetSeqLengthProt` and `x_GetSequence` return
/// for a protein of `residues` stored bytes: a defined 32-bit signed truncation of the stored
/// residue count, so a protein of `2^31` or more residues wraps to a negative `i32`. Nucleotide
/// lengths do not pass through here (their four-bases-per-byte overflow in NCBI is undefined
/// behaviour and LOSAT keeps them widened).
#[inline]
pub(crate) const fn ncbi_protein_length(residues: u32) -> i32 {
    residues as i32
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:302-323 (cited at
// `Volume::nucleotide_length_bounds`)
/// `GetSeqLengthExact`'s length of a nucleotide whose packed bases start at `start` and whose
/// ambiguity data start at `amb_start` (`start < amb_start`), the last byte of its bases being
/// `last`: the whole bytes times four plus the count in the last byte's low two bits.
fn packed_length(start: u64, amb_start: u64, last: u8) -> u64 {
    (amb_start - start - 1) * 4 + u64::from(last & 3)
}

impl Volume {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:104-109
    /// Opens the volume `base` (a path without extension) by reading its
    /// index file. The sequence and header files are opened on first use.
    /// EXPERIMENT (LOSAT_X_DBIDXREAD): the index file stays open and its arrays
    /// are read where they are needed (`IndexFile::open`); without the switch
    /// the whole file is read and its arrays are held (`IndexFile::parse`).
    pub fn open(base: &Path, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let mode = super::x_dbidxread();
        if mode != XDbIdxRead::Off {
            crate::utils::xreach::x_reach!("LOSAT_X_DBIDXREAD", true);
        }
        Self::open_as(base, seqtype, mode)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:90
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:3311-3318
    // Opening a volume reads its index file only (seqdbvol.cpp:104-109, above `Volume::open`). Its
    // column files (`.?Xa` index, `.?Xb` / `.?Xc` data: the masks of `makeblastdb -mask_data`,
    // user columns) are opened by `x_OpenAllColumns`, which only `GetColumnBlob`,
    // `GetColumnMetaData`, `ListColumns` and `GetColumnId` call (seqdbvol.cpp:3282-3283,
    // 3300-3301, 3316-3317, 3399-3400). A search calls them only for a masking algorithm that
    // `-db_soft_mask` or `-db_hard_mask` names (blast_args.cpp:2519-2522), or an imported search
    // strategy (`-import_search_strategy`, blast_app_util.cpp:498-511; blast_formatter.cpp:370
    // likewise); without one `GetFilteringAlgorithm` is -1 (uniform_search.cpp:198, 225,
    // 303-308; seqsrc_seqdb.cpp:318-320, 727-729; seqinfosrc_seqdb.cpp:113-123;
    // align_format_util.cpp:396-403). LOSAT refuses the imported strategy and the database masks
    // but BLASTN's (series NP, M4a-2), so a search without them reads no column file, whatever
    // the files hold; BLASTN's masked search, `losat blastdbcmd` and `losat blast_formatter` with
    // an archive's filtering algorithm (`blast_formatter::archive`, the formatter's
    // `GetMaskAlgorithmDetails`) open them on first use ([`Volume::columns`]).
    /// `open` with the index read as `mode` says.
    pub(super) fn open_as(
        base: &Path,
        seqtype: SeqType,
        mode: XDbIdxRead,
    ) -> Result<Self, SeqDbError> {
        let index_path = with_extension(base, seqtype, "in");
        let index_file = DbFile::open(&index_path)?
            .read_in_order()
            .zero_pad_last_page();
        let index = match mode {
            XDbIdxRead::Off => IndexFile::parse(&index_file.read_all()?, seqtype, &index_path)?,
            XDbIdxRead::On => IndexFile::open(index_file, seqtype)?,
        };
        Ok(Self {
            base: base.to_path_buf(),
            index,
            seq: OnceLock::new(),
            hdr: OnceLock::new(),
            defline_rule: None,
            exclude_model: false,
            alias_lines: None,
            oid_mask_type: 0,
            memb_bit: 0,
            columns: OnceLock::new(),
            read_failure: None,
            read_verification: None,
        })
    }

    /// The columns of the volume, opened on first use (`CSeqDBVol::x_OpenAllColumns`); none
    /// for a volume without column files. A failed open is repeated by the next call.
    pub fn columns(&self) -> Result<&[Column], SeqDbError> {
        if let Some(columns) = self.columns.get() {
            return Ok(columns);
        }
        let columns = column::open_all(&self.base, self.index.seqtype, self.index.num_oids)?;
        Ok(self.columns.get_or_init(|| columns))
    }

    /// The column titled `BlastDb/MaskData` (`CSeqDBVol::GetColumnId`), if the volume has one.
    pub fn mask_column(&self) -> Result<Option<&Column>, SeqDbError> {
        Ok(self
            .columns()?
            .iter()
            .find(|column| column.title() == MASK_DATA_TITLE))
    }

    /// Path of the volume without extension.
    pub fn base(&self) -> &Path {
        &self.base
    }

    /// The parsed index file.
    pub fn index(&self) -> &IndexFile {
        &self.index
    }

    /// Number of sequences in the volume.
    pub fn num_oids(&self) -> u32 {
        self.index.num_oids
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:117-134
    fn lazy_file<'a>(
        &self,
        slot: &'a OnceLock<DbFile>,
        ext: &str,
    ) -> Result<&'a DbFile, SeqDbError> {
        if let Some(file) = slot.get() {
            return Ok(file);
        }
        let file = DbFile::open(&with_extension(&self.base, self.index.seqtype, ext))?
            .read_in_order()
            .zero_pad_last_page();
        // Two threads may both open the file; the first one stored is used.
        Ok(slot.get_or_init(|| file))
    }

    fn seq_file(&self) -> Result<&DbFile, SeqDbError> {
        if let Some(failure) = &self.read_failure {
            return Err(failure.to_error());
        }
        self.lazy_file(&self.seq, "sq")
    }

    /// Makes every read of the sequence file fail with `failure`, the failure of NCBI's OID
    /// list of the database, before which a search reads no sequence
    /// (`SeqDb::fail_search_reads_with_oid_list`).
    pub(super) fn fail_sequence_reads(&mut self, failure: super::oid_list::OidListError) {
        self.read_failure = Some(failure);
    }

    /// Opens the sequence file, as the first read of a sequence does (`seq_file`): its error
    /// where it cannot be opened, as NCBI's first read of the volume throws it.
    pub fn open_sequence_file(&self) -> Result<(), SeqDbError> {
        self.seq_file().map(|_| ())
    }

    pub(super) fn hdr_file(&self) -> Result<&DbFile, SeqDbError> {
        self.lazy_file(&self.hdr, "hr")
    }

    fn check_range(&self, first: u32, end: u32) -> Result<(), SeqDbError> {
        if first >= end || end > self.index.num_oids {
            return Err(SeqDbError::Range(format!(
                "sequences {first}..{end} are not a non-empty range of the {} sequences of {}",
                self.index.num_oids,
                self.base.display()
            )));
        }
        Ok(())
    }

    /// Bytes before the first sequence that a block reads with it: the
    /// protein sentinel byte.
    fn lead(&self) -> u64 {
        match self.index.seqtype {
            SeqType::Protein => 1,
            SeqType::Nucleotide => 0,
        }
    }

    /// The number of residues of protein `oid` (volume-relative), read
    /// from the index alone: `read_block(oid, oid + 1, _)`'s `protein(0).len()`.
    pub fn protein_length(&self, oid: u32) -> Result<u32, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein length of a nucleotide volume"
        );
        self.check_range(oid, oid.saturating_add(1))?;
        // One read of the two offsets, as `read_block` reads them, into a pair rather than a
        // vector (no allocation for each sequence); the residues end before the sentinel byte at
        // the start of the next sequence's offset.
        let (start, end) = self.index.sequence_range(oid)?;
        end.checked_sub(start)
            .and_then(|bytes| bytes.checked_sub(1))
            .and_then(|residues| u32::try_from(residues).ok())
            .ok_or_else(|| self.index.sequence_changed(oid, 2))
    }

    fn protein_seq_bounds(&self, oid: u32) -> Result<(u64, u64), SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        self.index.sequence_range(oid)
    }

    /// Opens the volume's sequence file and refuses 0-byte files where NCBI maps a null pointer.
    fn open_non_empty_seq_file(&self) -> Result<&DbFile, SeqDbError> {
        let file = self.seq_file()?;
        if file.is_empty() {
            return Err(SeqDbError::NcbiCrashes {
                path: file.path().to_path_buf(),
                message: "the sequence file has no bytes, which NCBI maps as a null pointer"
                    .to_string(),
            });
        }
        Ok(file)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1784-1804
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:755-758
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:699-707
    /// Relative index of the first protein in `first..end` whose NCBI `x_GetSequence` length is
    /// negative (signed 32-bit `end - 1 - start < 0`), or `None`.
    pub(super) fn x_first_negative_protein_in(
        &self,
        first: u32,
        end: u32,
    ) -> Result<Option<u32>, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein lengths of a nucleotide volume"
        );
        self.check_range(first, end)?;
        let offsets = self.index.sequence_offsets(first, end)?;
        for (i, pair) in offsets.windows(2).enumerate() {
            let len_i32 = (pair[1] as i32)
                .wrapping_sub(1)
                .wrapping_sub(pair[0] as i32);
            if len_i32 < 0 {
                return Ok(Some(i as u32));
            }
        }
        Ok(None)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1817-1841
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:694-710
    /// Relative index of the first nucleotide in `first..end` whose NCBI `x_GetSequence` length is
    /// negative, or `None`. Checks that the byte NCBI reads (`file[start + whole_bytes]`) lies
    /// within the mapped readable length for every scanned OID, reading the byte only when
    /// `whole_bytes >= 2^29` can overflow 32-bit `whole_bytes * 4`.
    pub(super) fn x_first_negative_nucleotide_in(
        &self,
        first: u32,
        end: u32,
        first_oid: u32,
    ) -> Result<Option<u32>, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Nucleotide,
            "nucleotide lengths of a protein volume"
        );
        self.check_range(first, end)?;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        let file = self.open_non_empty_seq_file()?;
        let readable = file.readable_len();
        for (i, &amb) in amb_offsets.iter().enumerate() {
            let start = seq_offsets[i];
            let whole_bytes = (amb as i32).wrapping_sub(start as i32).wrapping_sub(1);
            let last_pos = i64::from(start) + i64::from(whole_bytes);
            if last_pos < 0 {
                return Err(SeqDbError::NcbiCrashes {
                    path: file.path().to_path_buf(),
                    message: format!(
                        "sequence {}: ambiguity offset {amb} causes NCBI to read before the \
                         start of the sequence file mapping",
                        u64::from(first_oid) + i as u64
                    ),
                });
            }
            if (last_pos as u64) >= readable {
                file.x_check_range(last_pos as u64, last_pos as u64 + 1)?;
            }
            if !(0..(1 << 29)).contains(&whole_bytes) {
                if whole_bytes < 0 {
                    return Ok(Some(i as u32));
                }
                let last = file.with_range(last_pos as u64, last_pos as u64 + 1, |b| b[0])?;
                let len_i32 = whole_bytes
                    .wrapping_mul(4)
                    .wrapping_add(i32::from(last & 3));
                if len_i32 < 0 {
                    return Ok(Some(i as u32));
                }
            }
        }
        Ok(None)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1783-1785
    /// Opens the volume's sequence file (matching `x_OpenSeqFile` before `x_GetSequence`) and
    /// returns an empty `SequenceBlock` marked `stopped_at_negative`.
    pub(super) fn empty_stopped_block(
        &self,
        first_ordinal: u32,
        first_oid: u32,
    ) -> Result<SequenceBlock, SeqDbError> {
        let path = self.open_non_empty_seq_file()?.path().to_path_buf();
        Ok(SequenceBlock::empty_stopped(
            self.index.seqtype,
            first_ordinal,
            first_oid,
            path,
        ))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
    /// The signed 32-bit length of protein `oid` (volume-relative) as NCBI's `GetSeqLengthProt`
    /// returns it: `end_offset - 1 - start_offset` in 32-bit signed arithmetic (negative for
    /// `2^31` or more residues or non-ascending offsets).
    pub fn protein_length_i32(&self, oid: u32) -> Result<i32, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein length of a nucleotide volume"
        );
        let (start, end) = self.protein_seq_bounds(oid)?;
        Ok((end as i32).wrapping_sub(1).wrapping_sub(start as i32))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1783-1804
    /// Opens `.psq` and returns NCBI's `x_GetSequence` signed 32-bit length for protein `oid`
    /// (volume-relative). When `length <= 0`, NCBI performs only pointer arithmetic without
    /// dereferencing `.psq`; when `length > 0`, `start == 0` is refused because the consumer reads
    /// the leading sentinel before the mapping, and `[start - 1 .. start + length + 1]` is checked
    /// against the mapped 4 KiB page. Refuses 0-byte `.psq` files where NCBI crashes.
    pub fn protein_sequence_len_i32(&self, oid: u32) -> Result<i32, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein sequence length of a nucleotide volume"
        );
        let bounds = self.protein_seq_bounds(oid)?;
        let (start, end) = bounds;
        let file = self.open_non_empty_seq_file()?;
        let len = (end as i32).wrapping_sub(1).wrapping_sub(start as i32);
        if len <= 0 {
            return Ok(len);
        }
        let Some(lead_start) = start.checked_sub(1) else {
            return Err(SeqDbError::NcbiCrashes {
                path: file.path().to_path_buf(),
                message: format!(
                    "sequence {oid}: start offset 0 causes NCBI to read before the start of the \
                     sequence file mapping"
                ),
            });
        };
        file.x_check_range(lead_start, start + len as u64 + 1)?;
        Ok(len)
    }

    fn checked_protein_range(&self, oid: u32) -> Result<(&DbFile, u64, u64), SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein range check of a nucleotide volume"
        );
        self.check_range(oid, oid.saturating_add(1))?;
        let (start, end) = self.index.sequence_range(oid)?;
        if start < 1 || end <= start {
            return Err(self.index.sequence_changed(oid, 2));
        }
        let file = self.seq_file()?;
        file.x_check_range(start - 1, end)?;
        Ok((file, start, end))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1784-1804
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:755-758
    /// Opens the protein sequence file and validates that protein `oid` (volume-relative) lies
    /// within the file without reading its bytes when `length > 0`; when `length <= 0`, NCBI's
    /// `x_GetSequence` performs only pointer arithmetic without dereferencing `.psq`.
    pub fn check_protein_range(&self, oid: u32) -> Result<(), SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein range check of a nucleotide volume"
        );
        self.check_range(oid, oid.saturating_add(1))?;
        let (start, end) = self.index.sequence_range(oid)?;
        let file = self.open_non_empty_seq_file()?;
        let len = (end as i32).wrapping_sub(1).wrapping_sub(start as i32);
        if len <= 0 {
            return Ok(());
        }
        let Some(lead_start) = start.checked_sub(1) else {
            return Err(SeqDbError::NcbiCrashes {
                path: file.path().to_path_buf(),
                message: format!(
                    "sequence {oid}: start offset 0 causes NCBI to read before the start of the \
                     sequence file mapping"
                ),
            });
        };
        file.x_check_range(lead_start, start + len as u64 + 1)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1630-1635
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1784-1804
    /// Reads `len` NCBIstdaa residue bytes starting at residue `offset` of protein `oid`
    /// (volume-relative) without reading the full sequence: validates the full protein range in
    /// `.psq` first (matching `x_GetSequence`) and then reads only `[offset..offset + len]`.
    pub fn read_protein_slice(
        &self,
        oid: u32,
        offset: u64,
        len: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        let (file, start, end) = self.checked_protein_range(oid)?;
        let raw_len = end - start - 1;
        let slice_end = offset
            .checked_add(len as u64)
            .filter(|&end_off| end_off <= raw_len)
            .ok_or_else(|| {
                SeqDbError::Range(format!(
                    "protein slice {offset}..+{len} exceeds sequence length {raw_len} in {}",
                    self.base.display()
                ))
            })?;
        file.read_range(start + offset, start + slice_end)
    }

    /// Bytes of the sequence file that a block of the sequences `first..end`
    /// occupies (protein: with the sentinel bytes around them; nucleotide:
    /// with their ambiguity data).
    #[cfg(test)]
    pub fn block_bytes(&self, first: u32, end: u32) -> Result<u64, SeqDbError> {
        let (start, stop) = (
            self.index.sequence_start(first)?,
            self.index.sequence_start(end)?,
        );
        // Ascending when the arrays are held (checked at open); read from the file, the two
        // entries can disagree only if the file changed after the open.
        let bytes = stop.checked_sub(start).ok_or_else(|| {
            self.index
                .sequence_changed(first, (end - first) as usize + 1)
        })?;
        Ok(bytes + self.lead())
    }

    /// The end of the block that starts at `first` (`SeqDb::next_block`): the
    /// last end in `first + 1 ..= num_oids` whose `block_bytes` fit
    /// `max_bytes`, at least `first + 1`. The offsets ascend, so the ends that
    /// fit form a prefix.
    pub(super) fn block_end(&self, first: u32, max_bytes: u64) -> Result<u32, SeqDbError> {
        if let Some(offsets) = self.index.held_sequence_offsets() {
            return Ok(held_block_end(offsets, first, max_bytes, self.lead()));
        }
        self.x_read_block_end(first, max_bytes, X_BLOCK_END_WINDOW)
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): `block_end` over the offsets in the
    /// file: the end of the volume when it fits (one entry read), else the
    /// offsets after `first` read forward in windows of `window` entries,
    /// doubled up to `X_BLOCK_END_MAX_WINDOW`, to the first end that does not
    /// fit.
    pub(super) fn x_read_block_end(
        &self,
        first: u32,
        max_bytes: u64,
        window: usize,
    ) -> Result<u32, SeqDbError> {
        let n = self.num_oids();
        let lead = self.lead();
        if max_bytes < lead {
            return Ok(first + 1);
        }
        // `block_bytes(first, end) <= max_bytes` is `sequence_start(end) <= limit`.
        let first_start = self.index.sequence_start(first)?;
        let limit = first_start.saturating_add(max_bytes - lead);
        let last_start = self.index.sequence_start(n)?;
        if last_start >= first_start && last_start <= limit {
            return Ok(n);
        }
        let mut window = window.max(1);
        let mut next = first + 1;
        let mut prev = first_start;
        loop {
            // The end `n` does not fit, so the scan stops at or before it, unless the file
            // changed after the open and the entry of `n` reads differently.
            if next > n {
                return Err(self.index.sequence_changed(first, (n - first) as usize + 1));
            }
            let count = window.min((n + 1 - next) as usize);
            let offsets = self.index.sequence_offsets(next, next + count as u32 - 1)?;
            for (at, &offset) in offsets.iter().enumerate() {
                let offset = u64::from(offset);
                if offset < prev || offset > limit {
                    return Ok((next + at as u32).saturating_sub(1).max(first + 1));
                }
                prev = offset;
            }
            next += count as u32;
            window = (window * 2).min(X_BLOCK_END_MAX_WINDOW);
        }
    }

    /// EXPERIMENT (LOSAT_X_BNDBSMALL): `block_end` for a list of consecutive blocks
    /// (`SeqDb::x_next_blocks`): the sequence offsets from `first` on are read `entries` at a time
    /// into `window` and kept for the blocks after it, where `x_read_block_end` reads two entries
    /// and a window of its own for every block. A block that passes the end of a window that
    /// starts at `first` is cut by `x_read_block_end`. The end is the one `block_end` gives
    /// (the same comparisons on the same entries); a read error leaves `window` empty.
    #[cfg(any(
        test,
        all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )
    ))]
    pub(super) fn x_listed_block_end(
        &self,
        first: u32,
        max_bytes: u64,
        window: &mut XOffsetWindow,
        entries: usize,
    ) -> Result<u32, SeqDbError> {
        if self.index.held_sequence_offsets().is_some() {
            return self.block_end(first, max_bytes);
        }
        let n = self.num_oids();
        let lead = self.lead();
        if max_bytes < lead {
            return Ok(first + 1);
        }
        loop {
            if !window.holds(first) {
                window.offsets.clear();
                let count = entries.max(1).min((n + 1 - first) as usize);
                window.offsets = self
                    .index
                    .sequence_offsets(first, first + count as u32 - 1)?;
                window.first = first;
            }
            let base = (first - window.first) as usize;
            let first_start = u64::from(window.offsets[base]);
            let limit = first_start.saturating_add(max_bytes - lead);
            let after = &window.offsets[base + 1..];
            let mut prev = first_start;
            for (at, &offset) in after.iter().enumerate() {
                let offset = u64::from(offset);
                if offset < prev || offset > limit {
                    return Ok((first + at as u32).max(first + 1));
                }
                prev = offset;
            }
            if window.first as usize + window.offsets.len() == n as usize + 1 {
                return Ok(n);
            }
            if base == 0 {
                return self.x_read_block_end(first, max_bytes, entries);
            }
            window.offsets.clear();
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1838
    // NCBI points into the mapped file one sequence at a time; LOSAT copies
    // the bytes of a contiguous range of sequences with one positioned read.
    // The protein range starts one byte early, so every sequence of the block
    // has its sentinel byte on both sides as in NCBI's buffer.
    /// Reads the sequences `first..end` (volume-relative, non-empty) as one
    /// block. `first_ordinal` is the database ordinal id of `first`, `first_oid` its OID in the
    /// database, which the block's messages name (the sequences of a block have consecutive OIDs).
    pub fn read_block(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_block_in(first, end, first_ordinal, first_oid, BlockBuffer::New)
    }

    /// Reads the protein sequences `first..end` (volume-relative, non-empty) for a database
    /// search: any protein whose NCBI 32-bit length ([`ncbi_protein_length`]) is not positive
    /// (`0` or `>= 2^31` residues, which NCBI's `s_SeqDbGetSequence` skips) has its file range
    /// validated without allocating its bytes and is represented as an empty `[0, 0]` sentinel
    /// pair in the returned block.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1784-1804
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:375-378
    pub fn read_protein_search_block(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_protein_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::New,
            false,
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:339-344,699-707
    /// [`Self::read_protein_search_block`] with `stop_at_negative` for multithreaded DB-split
    /// scans where NCBI's `x_FillSeqBuffer` stops at the first OID whose `GetSequence` length is
    /// negative.
    pub(super) fn read_protein_search_block_mt(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        stop_at_negative: bool,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_protein_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::New,
            stop_at_negative,
        )
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `read_block` with the bytes in a buffer of
    /// `File::x_read_range_mapped`.
    pub fn x_read_block_mapped(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_block_in(first, end, first_ordinal, first_oid, BlockBuffer::Mapped)
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): `read_block` with the bytes in the
    /// worker's buffer `kept` (`File::x_read_range_into`; `SequenceBlock::x_into_kept` gives the
    /// buffer back).
    pub fn x_read_block_into(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        kept: &mut Vec<u8>,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::Given(kept),
        )
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): [`Self::read_protein_search_block`] with
    /// the bytes in the worker's buffer `kept` when no sequence in the block has `>= 2^31`
    /// residues.
    pub fn x_read_protein_search_block_into(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        kept: &mut Vec<u8>,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_protein_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::Given(kept),
            false,
        )
    }

    pub(super) fn x_read_protein_search_block_into_mt(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        kept: &mut Vec<u8>,
        stop_at_negative: bool,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.read_protein_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::Given(kept),
            stop_at_negative,
        )
    }

    fn read_protein_search_block_in(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        buffer: BlockBuffer<'_>,
        stop_at_negative: bool,
    ) -> std::result::Result<SequenceBlock, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein search block of a nucleotide volume"
        );
        self.check_range(first, end)?;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let mut all_positive = seq_offsets[0] >= 1;
        let mut all_non_negative = all_positive;
        for pair in seq_offsets.windows(2) {
            let len = pair[1]
                .checked_sub(pair[0])
                .and_then(|bytes| bytes.checked_sub(1))
                .map(ncbi_protein_length);
            match len {
                Some(n) if n > 0 => {}
                Some(0) => all_positive = false,
                _ => {
                    all_positive = false;
                    all_non_negative = false;
                    break;
                }
            }
        }
        if all_positive
            || (all_non_negative
                && u64::from(seq_offsets[seq_offsets.len() - 1]) <= self.seq_file()?.readable_len())
        {
            return self.read_checked_block_with_offsets(
                first_ordinal,
                first_oid,
                seq_offsets,
                Vec::new(),
                buffer,
            );
        }
        // Cold path: at least one sequence has NCBI 32-bit length < 0 (>= 2^31 residues or
        // non-monotonic offsets), start == 0, or an offset past readable_len. Follow NCBI's
        // `x_GetSequence` / `s_SeqDbGetSequence` / `x_FillSeqBuffer` order: for length <= 0 NCBI
        // only computes pointer arithmetic without dereferencing `.psq`; negative lengths end
        // iteration in multithreaded DB-split mode (`seqdbimpl.cpp:699-707`), length 0 is skipped
        // in both single-thread and multithreaded modes, and positive length with `start == 0` is
        // refused because the consumer reads the leading sentinel before the mapping.
        let file = self.open_non_empty_seq_file()?;
        let count = seq_offsets.len() - 1;
        let mut data = vec![0u8];
        let mut packed_offsets = Vec::with_capacity(seq_offsets.len());
        packed_offsets.push(1u32);
        let mut prev_skipped = false;
        let mut stopped_at_negative = false;
        for i in 0..count {
            let start = seq_offsets[i];
            let next = seq_offsets[i + 1];
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1804
            // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:755-758
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:699-707
            let len_i32 = (next as i32).wrapping_sub(1).wrapping_sub(start as i32);
            if len_i32 < 0 {
                if stop_at_negative {
                    stopped_at_negative = true;
                    break;
                }
                data.push(0);
                packed_offsets.push(data.len() as u32);
                prev_skipped = true;
            } else if len_i32 == 0 {
                data.push(0);
                packed_offsets.push(data.len() as u32);
                prev_skipped = true;
            } else {
                let Some(lead_start) = u64::from(start).checked_sub(1) else {
                    return Err(SeqDbError::NcbiCrashes {
                        path: file.path().to_path_buf(),
                        message: format!(
                            "sequence {}: start offset 0 causes NCBI to read before the start of \
                             the sequence file mapping",
                            u64::from(first_oid) + i as u64
                        ),
                    });
                };
                let raw = file.read_range(lead_start, u64::from(start) + len_i32 as u64 + 1)?;
                if i == 0 || prev_skipped {
                    let last = data.len() - 1;
                    data[last] = raw[0];
                }
                data.extend_from_slice(&raw[1..]);
                packed_offsets.push(data.len() as u32);
                prev_skipped = false;
            }
        }
        Ok(SequenceBlock {
            seqtype: SeqType::Protein,
            first_ordinal,
            first_oid,
            path: file.path().to_path_buf(),
            data,
            data_start: 0,
            seq_offsets: packed_offsets,
            amb_offsets: Vec::new(),
            stopped_at_negative,
            amb_error: None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1805-1841,2350-2363
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:339-344,699-707
    /// Reads the nucleotide sequences `first..end` (volume-relative, non-empty) for a search or
    /// retrieval: normal ascending blocks take the contiguous fast path; irregular ambiguity
    /// offsets (`amb <= start` with `amb >= 1`, or `amb > next`) follow NCBI's `x_GetSequence`
    /// and defer `x_GetAmbChar`'s `File error: could not get ambiguity data.` until ambiguity
    /// data are requested. Offsets that read before the mapping (`amb == 0`) are refused.
    pub(super) fn read_nucleotide_search_block_mt(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        stop_at_negative: bool,
    ) -> SeqBlockResult {
        self.read_nucleotide_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::New,
            stop_at_negative,
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1805-1841
    /// [`Self::read_nucleotide_search_block_mt`] reusing the caller's buffer `kept`.
    pub(super) fn x_read_nucleotide_search_block_into_mt(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        kept: &mut Vec<u8>,
        stop_at_negative: bool,
    ) -> SeqBlockResult {
        self.read_nucleotide_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::Given(kept),
            stop_at_negative,
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1817-1841,2350-2363
    /// [`Self::read_nucleotide_search_block_mt`] into a mapped buffer for traceback / retrieval.
    pub(super) fn x_read_nucleotide_search_block_mapped(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
    ) -> SeqBlockResult {
        self.read_nucleotide_search_block_in(
            first,
            end,
            first_ordinal,
            first_oid,
            BlockBuffer::Mapped,
            false,
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1817-1839
    fn read_nucleotide_search_block_in(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        buffer: BlockBuffer<'_>,
        stop_at_negative: bool,
    ) -> SeqBlockResult {
        assert!(
            self.index.seqtype == SeqType::Nucleotide,
            "nucleotide search block of a protein volume"
        );
        self.check_range(first, end)?;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        // `seq_offsets[i] < amb <= seq_offsets[i + 1]` for every `i` implies `seq_offsets` is
        // strictly ascending, so one pass validates both offset arrays for the contiguous fast path.
        let all_regular = amb_offsets
            .iter()
            .enumerate()
            .all(|(i, &amb)| seq_offsets[i] < amb && amb <= seq_offsets[i + 1]);
        if all_regular {
            return self.read_checked_block_with_offsets(
                first_ordinal,
                first_oid,
                seq_offsets,
                amb_offsets,
                buffer,
            );
        }
        self.read_nucleotide_irregular_block(
            first_ordinal,
            first_oid,
            &seq_offsets,
            &amb_offsets,
            stop_at_negative,
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1817-1841,2358-2363
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:699-707
    fn read_nucleotide_irregular_block(
        &self,
        first_ordinal: u32,
        first_oid: u32,
        seq_offsets: &[u32],
        amb_offsets: &[u32],
        stop_at_negative: bool,
    ) -> Result<SequenceBlock, SeqDbError> {
        let file = self.open_non_empty_seq_file()?;
        let count = seq_offsets.len() - 1;
        let mut data = Vec::new();
        let mut packed_seq = Vec::with_capacity(seq_offsets.len());
        let mut packed_amb = Vec::with_capacity(count);
        let mut amb_invalid = Vec::with_capacity(count);
        let mut any_amb_invalid = false;
        let mut stopped_at_negative = false;
        packed_seq.push(0u32);
        for i in 0..count {
            let (start, amb, next) = (seq_offsets[i], amb_offsets[i], seq_offsets[i + 1]);
            let whole_bytes = (amb as i32).wrapping_sub(start as i32).wrapping_sub(1);
            let last_pos = i64::from(start) + i64::from(whole_bytes);
            if last_pos < 0 {
                return Err(SeqDbError::NcbiCrashes {
                    path: file.path().to_path_buf(),
                    message: format!(
                        "sequence {}: ambiguity offset {amb} causes NCBI to read before the start \
                         of the sequence file mapping",
                        u64::from(first_oid) + i as u64
                    ),
                });
            }
            // NCBI's `x_GetSequence` (`seqdbvol.cpp:1817-1838`) does pointer arithmetic on `start`
            // and dereferences only `file[start + whole_bytes]` to compute the length.
            file.x_check_range(last_pos as u64, last_pos as u64 + 1)?;
            if whole_bytes < 0 {
                if stop_at_negative {
                    stopped_at_negative = true;
                    break;
                }
                data.push(0);
                let end_off = data.len() as u32;
                packed_amb.push(end_off);
                packed_seq.push(end_off);
                amb_invalid.push(false);
            } else {
                let packed_end = u64::from(start) + whole_bytes as u64 + 1;
                let packed = file.read_range(u64::from(start), packed_end)?;
                data.extend_from_slice(&packed);
                packed_amb.push(data.len() as u32);
                if amb <= next {
                    if amb < next {
                        let amb_bytes = file.read_range(u64::from(amb), u64::from(next))?;
                        data.extend_from_slice(&amb_bytes);
                    }
                    amb_invalid.push(false);
                } else {
                    any_amb_invalid = true;
                    amb_invalid.push(true);
                }
                packed_seq.push(data.len() as u32);
            }
        }
        Ok(SequenceBlock {
            seqtype: SeqType::Nucleotide,
            first_ordinal,
            first_oid,
            path: file.path().to_path_buf(),
            data,
            data_start: 0,
            seq_offsets: packed_seq,
            amb_offsets: packed_amb,
            stopped_at_negative,
            amb_error: any_amb_invalid.then(|| amb_invalid.into_boxed_slice()),
        })
    }

    fn read_block_in(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
        buffer: BlockBuffer<'_>,
    ) -> Result<SequenceBlock, SeqDbError> {
        self.check_range(first, end)?;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        self.read_block_with_offsets(
            first,
            first_ordinal,
            first_oid,
            seq_offsets,
            amb_offsets,
            buffer,
        )
    }

    fn read_block_with_offsets(
        &self,
        first: u32,
        first_ordinal: u32,
        first_oid: u32,
        seq_offsets: Vec<u32>,
        amb_offsets: Vec<u32>,
        buffer: BlockBuffer<'_>,
    ) -> Result<SequenceBlock, SeqDbError> {
        let seqtype = self.index.seqtype;
        let valid_seq = match seqtype {
            SeqType::Protein => {
                seq_offsets[0] >= 1 && seq_offsets.windows(2).all(|pair| pair[0] < pair[1])
            }
            SeqType::Nucleotide => seq_offsets.windows(2).all(|pair| pair[0] <= pair[1]),
        };
        if !valid_seq {
            return Err(self.index.sequence_changed(first, seq_offsets.len()));
        }
        // For `amb <= start`, NCBI reads byte `start - 1` (inside the mapping when `amb >= 1`,
        // before the mapping when `amb == 0`); for `amb > next`, NCBI reads into the next
        // sequence's bytes (and multithreaded `blastn` aborts on traceback). `read_block` (used
        // by `tblastn` and `tblastx`) refuses any `!(start < amb <= next)` (see `refusals.md`);
        // `blastn`, `blastdbcmd`, and `blastdbcheck` use `read_nucleotide_search_block` /
        // `nucleotide_length_i32`.
        if !amb_offsets
            .iter()
            .enumerate()
            .all(|(i, &amb)| seq_offsets[i] < amb && amb <= seq_offsets[i + 1])
        {
            return Err(self.index.ambiguity_changed(first, amb_offsets.len()));
        }
        self.read_checked_block_with_offsets(
            first_ordinal,
            first_oid,
            seq_offsets,
            amb_offsets,
            buffer,
        )
    }

    fn read_checked_block_with_offsets(
        &self,
        first_ordinal: u32,
        first_oid: u32,
        seq_offsets: Vec<u32>,
        amb_offsets: Vec<u32>,
        buffer: BlockBuffer<'_>,
    ) -> Result<SequenceBlock, SeqDbError> {
        let seqtype = self.index.seqtype;
        let data_start = u64::from(seq_offsets[0]) - self.lead();
        let data_end = u64::from(seq_offsets[seq_offsets.len() - 1]);
        let file = self.seq_file()?;
        let data = match buffer {
            BlockBuffer::New => file.read_range(data_start, data_end)?,
            BlockBuffer::Mapped => file.x_read_range_mapped(data_start, data_end)?,
            BlockBuffer::Given(kept) => file.x_read_range_into(data_start, data_end, kept)?,
        };
        let block = SequenceBlock {
            seqtype,
            first_ordinal,
            first_oid,
            path: file.path().to_path_buf(),
            data,
            data_start,
            seq_offsets,
            amb_offsets,
            stopped_at_negative: false,
            amb_error: None,
        };
        Ok(block)
    }

    /// What `read_block(first, end, first_ordinal, first_oid)` of one
    /// nucleotide (`end == first + 1`) gives `nucleotide_len(0)` and `ambiguity_runs(0)`, without
    /// its 2-bit data: the index reads and checks of `read_block` in its order and the checks of
    /// the file range it reads (`DbFile::x_check_range`), then the last byte of the 2-bit data and
    /// the ambiguity data only. Returns the length, with `nucleotide_len`'s error, and a block of
    /// the ambiguity data whose `ambiguity_runs(0)` is the one of `read_block`'s block. An I/O
    /// error in the bytes between is not seen here; the search reads them later.
    pub fn x_nucleotide_check(
        &self,
        first: u32,
        end: u32,
        first_ordinal: u32,
        first_oid: u32,
    ) -> Result<(usize, SequenceBlock), SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Nucleotide && end == first.wrapping_add(1),
            "the check of one nucleotide"
        );
        self.check_range(first, end)?;
        let seq_offsets = self.index.sequence_offsets(first, end)?;
        let amb_offsets = self.index.ambiguity_offsets(first, end)?;
        let (start, amb_start, next) = (seq_offsets[0], amb_offsets[0], seq_offsets[1]);
        if start > next {
            return Err(self.index.sequence_changed(first, seq_offsets.len()));
        }
        if start >= amb_start || amb_start > next {
            return Err(self.index.ambiguity_changed(first, amb_offsets.len()));
        }
        let file = self.seq_file()?;
        file.x_check_range(u64::from(start), u64::from(next))?;
        let last = file.read_range(u64::from(amb_start) - 1, u64::from(amb_start))?;
        let length = u64::from(amb_start - start - 1) * 4 + u64::from(last[0] & 3);
        let ambiguity = SequenceBlock {
            seqtype: SeqType::Nucleotide,
            first_ordinal,
            first_oid,
            path: file.path().to_path_buf(),
            data: file.read_range(u64::from(amb_start), u64::from(next))?,
            data_start: u64::from(amb_start),
            // The block holds the ambiguity data only: its letters' range is
            // empty, and nothing but `ambiguity_runs` reads it.
            seq_offsets: vec![amb_start, next],
            amb_offsets: vec![amb_start],
            stopped_at_negative: false,
            amb_error: None,
        };
        if length > i32::MAX as u64 {
            return Err(ambiguity.format(0, "the sequence is longer than a signed 32-bit length"));
        }
        Ok((length as usize, ambiguity))
    }

    /// The numbers of residues of the proteins `first..end`
    /// (volume-relative, non-empty), read from the index alone with one read of their offsets:
    /// `read_block(first, end, _)`'s `protein(i).len()`.
    pub fn x_protein_lengths(&self, first: u32, end: u32) -> Result<Vec<u32>, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein lengths of a nucleotide volume"
        );
        self.check_range(first, end)?;
        let offsets = self.index.sequence_offsets(first, end)?;
        // The residues end before the sentinel byte at the start of the next sequence's offset;
        // read from the file, the entries can break that only if the file changed after the open.
        offsets
            .windows(2)
            .map(|pair| {
                pair[1]
                    .checked_sub(pair[0])
                    .and_then(|bytes| bytes.checked_sub(1))
            })
            .collect::<Option<Vec<u32>>>()
            .ok_or_else(|| self.index.sequence_changed(first, offsets.len()))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
    /// [`Self::x_protein_lengths`] as NCBI's signed 32-bit `GetSeqLengthProt` lengths
    /// (`end - 1 - start`, including negative lengths for `>= 2^31` or non-ascending offsets).
    pub fn x_protein_lengths_i32(&self, first: u32, end: u32) -> Result<Vec<i32>, SeqDbError> {
        assert!(
            self.index.seqtype == SeqType::Protein,
            "protein lengths of a nucleotide volume"
        );
        self.check_range(first, end)?;
        let offsets = self.index.sequence_offsets(first, end)?;
        Ok(offsets
            .windows(2)
            .map(|pair| {
                (pair[1] as i32)
                    .wrapping_sub(1)
                    .wrapping_sub(pair[0] as i32)
            })
            .collect())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2316
    /// The encoded `Blast-def-line-set` of the sequence `oid` (volume-relative;
    /// `with_header_record`'s).
    pub fn header_bytes(&self, oid: u32) -> Result<Vec<u8>, SeqDbError> {
        self.with_header_record(oid, |record, _| Ok(record.to_vec()))
    }

    /// `read` of the encoded `Blast-def-line-set` of the sequence `oid` (volume-relative) and the
    /// path of the header file, which the messages of its decode name: the open file's own, not a
    /// name built again for each record. The record is lent from the thread's buffer of the reads
    /// of a range (`DbFile::with_range`): a read of one header per sequence allocates no record.
    /// Under a setting that turns a check of NCBI's reader off, `read`'s error is the one
    /// `refused_read` makes; NCBI's default reads with `read` itself, so that a record that
    /// decodes costs nothing more than the test of the setting.
    fn with_header_record<T>(
        &self,
        oid: u32,
        read: impl FnOnce(&[u8], &Path) -> Result<T, SeqDbError>,
    ) -> Result<T, SeqDbError> {
        match self.read_verification.as_deref() {
            None => self.read_header_record(oid, read),
            Some(setting) => self.read_header_record(oid, |record, path| {
                read(record, path).map_err(|error| refused_read(setting, error))
            }),
        }
    }

    /// `with_header_record`'s read of the record.
    fn read_header_record<T>(
        &self,
        oid: u32,
        read: impl FnOnce(&[u8], &Path) -> Result<T, SeqDbError>,
    ) -> Result<T, SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let (start, end) = self.index.header_range(oid)?;
        let file = self.hdr_file()?;
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2208-2220
        // When `end < start`, `raw_data.size()` is non-zero (negative length cast to size), so
        // `CObjectIStreamAsnBinary` runs and immediately hits EOF in `Blast-def-line-set`.
        if end < start {
            file.x_check_range(start, start)?;
            return Err(SeqDbError::Ncbi(crate::seqid::ber::negative_length_eof()));
        }
        file.with_range(start, end, |record| read(record, file.path()))?
    }

    /// Cold-path collection of `Error: (803.7)` `Bad char [0x..]` posts emitted by NCBI's binary
    /// ASN.1 reader while reading the header record of `oid` (volume-relative), whether decoding
    /// succeeds or stops at a subsequent BER error.
    #[cold]
    #[inline(never)]
    pub fn header_replaced_char_posts(
        &self,
        oid: u32,
    ) -> Vec<crate::report::serial::encoder::ReplacedChar> {
        let Ok(()) = self.check_range(oid, oid.saturating_add(1)) else {
            return Vec::new();
        };
        let Ok((start, end)) = self.index.header_range(oid) else {
            return Vec::new();
        };
        if end <= start {
            return Vec::new();
        }
        let Ok(file) = self.hdr_file() else {
            return Vec::new();
        };
        file.with_range(start, end, crate::seqid::ber::replaced_char_posts)
            .unwrap_or_default()
    }

    /// Sets the run's settings of NCBI's checks of what it reads (`SeqDb::refuse_unchecked_reads`).
    pub(super) fn refuse_unchecked_reads(&mut self, setting: Arc<ReadVerification>) {
        self.read_verification = Some(setting);
    }

    /// The header of the sequence `oid` (volume-relative) from one read of its record
    /// (`SeqDb::read_header`).
    pub fn header(&self, oid: u32) -> Result<VolumeHeader, SeqDbError> {
        self.with_header_record(oid, |record, path| self.decode_header(oid, record, path))
    }

    /// The `taxid` of each definition line of the sequence `oid` (volume-relative), 0 for none.
    pub fn defline_taxids(&self, oid: u32) -> Result<Vec<i32>, SeqDbError> {
        self.with_header_record(oid, |record, path| {
            self.kept_lines(oid, record, path, decode_defline_taxids(record, path, oid))
        })
    }

    /// The definition lines of the sequence `oid` (volume-relative) with every member and
    /// every Seq-id (`header::decode_full`); a `general` id `BL_ORD_ID` keeps its
    /// volume-relative tag.
    pub fn defline_set(&self, oid: u32) -> Result<crate::seqid::DeflineSet, SeqDbError> {
        self.with_header_record(oid, |record, path| self.decode_kept_set(oid, record, path))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2202-2212
    /// `defline_set` written over `set`, whose buffers are kept (`header::decode_full_in`):
    /// `None` for a stored header of no bytes, which NCBI's read gives as no set
    /// (`x_GetHdrAsn1`), where `defline_set` gives a set of no lines (`SeqDb::header_set_in`);
    /// else whether the definition-line rule dropped a line of the record
    /// (`x_GetFilteredHeader`'s `changed`). With `record`, the stored record the set was decoded
    /// from is written there where no line was dropped (`SeqDb::defline_set_with_header_in`).
    pub fn header_set_in(
        &self,
        oid: u32,
        set: &mut crate::seqid::DeflineSet,
        record: Option<&mut Vec<u8>>,
    ) -> Result<Option<bool>, SeqDbError> {
        self.with_header_record(oid, |stored, path| {
            if stored.is_empty() {
                set.deflines.clear();
                set.utf8 = false;
                return Ok(None);
            }
            let dropped = self.decode_kept_set_in(oid, stored, path, set)?;
            if let Some(record) = record.filter(|_| !dropped) {
                record.clear();
                record.extend_from_slice(stored);
            }
            Ok(Some(dropped))
        })
    }

    /// Restricts the definition lines that the reads of this volume give to those `rule`
    /// keeps (`SeqDb::restrict_taxids`, `restrict_seqids`); `volume_start` is the volume's first
    /// OID in the database.
    pub(super) fn set_defline_rule(
        &mut self,
        rule: std::sync::Arc<DeflineRule>,
        volume_start: u32,
    ) {
        self.defline_rule = Some(LineRule { rule, volume_start });
    }

    /// Sets the volume's `OID_MASK_TYPE` (`SeqDb::read_alias_tree`) and with it whether the reads
    /// of this volume erase the definition lines of model sequences (`oid_mask::excludes_models`).
    pub(super) fn set_oid_mask_type(&mut self, mask_type: u32) {
        self.oid_mask_type = mask_type;
        self.exclude_model = super::oid_mask::excludes_models(mask_type);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1995-1999
    /// Whether the volume has an `OID_MASK_TYPE` other than 0, of any bits: NCBI's header read
    /// then filters the set it read (`SeqDb::filters_header_sets`).
    pub(super) fn has_oid_mask_type(&self) -> bool {
        self.oid_mask_type != 0
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:879-885
    /// Sets the volume's `MEMB_BIT` (`SeqDb::set_tree_oids`).
    pub(super) fn set_memb_bit(&mut self, memb_bit: u32) {
        self.memb_bit = memb_bit;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1995-1999
    /// Whether the volume has a non-zero `MEMB_BIT`: NCBI's header read then filters the set it
    /// read (`SeqDb::filters_header_sets`).
    pub(super) fn has_memb_bit(&self) -> bool {
        self.memb_bit != 0
    }

    /// Attaches the alias files' id lists that reach this volume (`SeqDb::read_alias_tree`), or
    /// clears them when `lists` is empty; `volume_start` is the volume's first OID in the database.
    pub(super) fn set_alias_lists(
        &mut self,
        lists: Vec<std::sync::Arc<super::id_list::AliasIdList>>,
        volume_start: u32,
    ) {
        self.alias_lines = (!lists.is_empty()).then_some(AliasLines {
            lists,
            volume_start,
        });
    }

    /// Whether the alias files attached id lists to this volume (`m_VolumeGiLists` not empty).
    pub(super) fn has_alias_lists(&self) -> bool {
        self.alias_lines.is_some()
    }

    /// Whether the reads of this volume remove definition lines: a list's rule (any list erases
    /// the lines without a Seq-id, a taxid list those of other taxids), the alias files' id
    /// lists, the membership bit (`MEMB_BIT`) or the erasure of model lines. Each may leave a
    /// sequence without a line (`SeqDb::drops_lineless_hits`).
    pub(super) fn filters_lines(&self) -> bool {
        self.defline_rule.is_some()
            || self.alias_lines.is_some()
            || self.exclude_model
            || self.memb_bit != 0
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2008-2022
    // NCBI tests `m_MemBit` on each line before `x_FilterHasId`, taxids and `OID_MASK_TYPE`: a
    // line without memberships or with an empty list fails, and otherwise only the first entry of
    // its memberships list is tested against the bit one below `m_MemBit`. Values above 32 shift
    // past a 32-bit signed integer in C++ (undefined behaviour), refused where a line carries
    // memberships whose test depends on the shift.
    /// Which definition lines of `set` pass the volume's `MEMB_BIT` (`m_MemBit`): `None` when
    /// `m_MemBit` is 0, else one flag per line.
    fn memb_mask(&self, set: &crate::seqid::DeflineSet) -> Result<Option<Vec<bool>>, SeqDbError> {
        if self.memb_bit == 0 {
            return Ok(None);
        }
        let bit_mask = (self.memb_bit <= 32).then(|| 1u32 << (self.memb_bit - 1));
        let mut mask = Vec::with_capacity(set.deflines.len());
        for defline in &set.deflines {
            let Some(&first) = defline.memberships.as_deref().and_then(<[i32]>::first) else {
                mask.push(false);
                continue;
            };
            let Some(bit_mask) = bit_mask else {
                return Err(SeqDbError::Unsupported {
                    path: self.base.clone(),
                    message: format!(
                        "alias file MEMB_BIT {} is past the 32 bits of a membership bitmask \
                         (NCBI shifts a 32-bit integer by {}, which is undefined behaviour)",
                        self.memb_bit,
                        self.memb_bit - 1,
                    ),
                });
            };
            mask.push(((first as u32) & bit_mask) != 0);
        }
        Ok(Some(mask))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1990-1995
    // This is the read of `GetBioseq` (seqdbvol.cpp:1324-1325); the overload with a stream, for
    // `GetSeqIDs`, does the same (2106-2111). The lists read the Seq-ids of a volume whose first
    // OID is not 0 with the tags of `BL_ORD_ID` made OIDs in the database
    // (`super::adjust_ordinal_tags`): a list's `gnl|BL_ORD_ID|n` keeps a line whose stored tag plus
    // the volume's first OID is n. The lookup that selected the OID read the stored number (the
    // accession table of `makeblastdb -input_type blastdb`, or a FASTA id `gnl|BL_ORD_ID|n` of
    // `-parse_seqids`), so on a list of databases or a database of several volumes the OID it
    // found may keep no line, and a line with another id may be kept or removed by the sum
    // (line DBID ordsl, DECISIONS 55). The read adjusts every line before any test, under every
    // list (a taxid list too): a string tag on a line that the list then erases ends the read
    // (review 2 G1 of line RG D1c-1b).
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2023-2034
    // The search's list and the alias files' lists attached to the volume test each line's
    // Seq-ids as NCBI's filtered header read does (`filter::line_mask`): an exception of
    // NCBI's lookups (a `ti` id whose string tag is no number, a `BL_ORD_ID` or `PIG` id with a
    // string tag) ends the read with NCBI's exception (`SeqDbError::Ncbi`). Lines that `memb_kept`
    // already rejected are skipped before any list lookup or taxid test.
    /// Which definition lines of the sequence `oid` (volume-relative), decoded as `set`, the
    /// membership bit (`memb_kept`), the search's list rule and the alias files' id lists keep, as
    /// NCBI's filtered header read tests them: `None` where a seqid rule keeps none; the reader's
    /// error for a tag NCBI's read cannot adjust, NCBI's exception of a lookup.
    fn list_mask(
        &self,
        oid: u32,
        set: &crate::seqid::DeflineSet,
        memb_kept: Option<&[bool]>,
    ) -> Result<Option<Vec<bool>>, SeqDbError> {
        let user = self.defline_rule.as_ref();
        let alias = self.alias_lines.as_ref();
        let volume_start = match (user, alias) {
            (Some(line_rule), _) => line_rule.volume_start,
            (None, Some(lines)) => lines.volume_start,
            (None, None) => {
                return Ok(Some(
                    memb_kept.map_or_else(|| vec![true; set.deflines.len()], <[bool]>::to_vec),
                ));
            }
        };
        let reads_seqids = alias.is_some() || user.is_some_and(|rule| rule.rule.reads_seqids());
        let adjusted;
        let set = if volume_start > 0 && super::has_ordinal_tag(set) {
            let mut copy = set.clone();
            match super::adjust_ordinal_tags(&mut copy, volume_start) {
                Ok(_) => {
                    adjusted = copy;
                    &adjusted
                }
                // A taxid rule compares no Seq-id: NCBI's sum past a signed 32-bit id raises no
                // exception and changes none of its tests, and the reports adjust the kept lines
                // again where they read them.
                Err(super::OrdinalTagError::PastSigned32) if !reads_seqids => set,
                Err(error) => return Err(error.error(&self.base, volume_start + oid)),
            }
        } else {
            set
        };
        let lists = alias.map_or(&[][..], |lines| lines.lists.as_slice());
        super::filter::line_mask(
            user.map(|line_rule| &*line_rule.rule),
            lists,
            set,
            memb_kept,
        )
        .map_err(SeqDbError::Ncbi)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2006-2058
    // `MEMB_BIT` runs first, then the lists' tests on the lines it keeps, and the erasure of model
    // lines last on the lines they keep; a line whose every test passes is kept. A sequence whose
    // every line is erased is read with none, as NCBI's `GetHdr` gives it (probes A2, A7a of line
    // RG D1b): its hits get no Seq-align (`SeqDb::line_read`), and `GetBioseq` gives none. The
    // accession guide's warning about an id (`oid_mask::is_predicted`) is posted at each of NCBI's
    // header reads that is not served by its cache of filtered headers (seqdbvol.cpp:1977-1988, a
    // cache of the main thread): LOSAT does not count those reads and refuses such an id. NCBI
    // decodes the whole header before it removes a line (seqdbvol.cpp:1993, then 2057), and its
    // decoder writes an error for each byte it reads as '#' in a line it then removes too (probe b1
    // of line RG D1b-2): refused as the bytes of a kept line are where they are read.
    /// Which definition lines of the sequence `oid` (volume-relative), decoded as `set`, the reads
    /// keep, as NCBI's filtered header read tests them (`MEMB_BIT`, then the search's list rule
    /// and the alias files' id lists, `list_mask`, then the erasure of model lines): `None` where
    /// a seqid rule keeps none; the reader's error for a tag NCBI's read cannot adjust, NCBI's
    /// exception of a lookup, LOSAT's refusals of a removed line with bytes NCBI reads as '#', of
    /// an id the accession guide warns about, and of `MEMB_BIT` past 32 bits on memberships.
    fn filtered_mask(
        &self,
        oid: u32,
        set: &crate::seqid::DeflineSet,
    ) -> Result<Option<Vec<bool>>, SeqDbError> {
        let memb_kept = self.memb_mask(set)?;
        let mut mask = match self.list_mask(oid, set, memb_kept.as_deref())? {
            Some(mask) => mask,
            None => {
                self.refuse_removed_replaced_bytes(oid, set, &vec![false; set.deflines.len()])?;
                return Ok(None);
            }
        };
        let mut warned = false;
        if self.exclude_model {
            for (keep, defline) in mask.iter_mut().zip(&set.deflines) {
                if *keep {
                    *keep = !super::oid_mask::is_model_line(&defline.ids, &mut |_| warned = true);
                }
            }
        }
        self.refuse_removed_replaced_bytes(oid, set, &mask)?;
        if warned {
            return Err(SeqDbError::Unsupported {
                path: self.base.clone(),
                message: format!(
                    "sequence {oid} of the volume: a Seq-id that NCBI's accession guide warns \
                     about, under the alias files' OID_MASK_TYPE (NCBI writes the warning at each \
                     of its uncached reads of the header)"
                ),
            });
        }
        Ok(Some(mask))
    }

    /// LOSAT's refusal of the sequence `oid` (volume-relative) where a line of `set` that `mask`
    /// removes has bytes NCBI reads as '#' (`filtered_mask`).
    fn refuse_removed_replaced_bytes(
        &self,
        oid: u32,
        set: &crate::seqid::DeflineSet,
        mask: &[bool],
    ) -> Result<(), SeqDbError> {
        let removed_replaced = set.deflines.iter().zip(mask).any(|(defline, &keep)| {
            !keep && (defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        });
        if removed_replaced {
            let what = if self.memb_bit != 0 {
                "a list, MEMB_BIT or the alias files' OID_MASK_TYPE"
            } else {
                "a list or the alias files' OID_MASK_TYPE"
            };
            return Err(SeqDbError::Unsupported {
                path: self.base.clone(),
                message: format!(
                    "sequence {oid} of the volume: a definition line that {what} removes, with \
                     bytes NCBI reads as '#' (it writes an error each time it reads them)"
                ),
            });
        }
        Ok(())
    }

    /// [`Volume::filtered_mask`] for a read of the sequence `oid` (volume-relative): LOSAT's
    /// refusal (`no_line_kept`) where a seqid rule keeps no line.
    fn line_mask(&self, oid: u32, set: &crate::seqid::DeflineSet) -> Result<Vec<bool>, SeqDbError> {
        self.filtered_mask(oid, set)?
            .ok_or_else(|| self.no_line_kept(oid))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1909-1946
    /// Whether the volume's rule fills NCBI's `CSeqDBIdSet` (`x_InitIdSet`): a negative seqid
    /// list (also one version 5 processing emptied), a positive GI list, or a negative GI list
    /// (`SeqDb::excludes_lineless_at_fetch`).
    pub(super) fn has_fetch_exclusion_list(&self) -> bool {
        self.defline_rule.as_ref().is_some_and(|line_rule| {
            matches!(
                &*line_rule.rule,
                super::filter::DeflineRule::NegativeSeqids(_)
                    | super::filter::DeflineRule::Gis(_)
                    | super::filter::DeflineRule::NegativeGis(_)
            )
        })
    }

    /// Whether the volume's reads keep a line with a Seq-id of the sequence `oid`
    /// (volume-relative; `SeqDb::line_read`): false where the list's rule or the alias files'
    /// lists keep none with a Seq-id (an alias `TAXIDLIST` can rescue a line without a Seq-id,
    /// which still gives `GetSeqIDs` nothing) or the erasure of model lines leaves none, and where
    /// the volume filters no line, where none of its lines holds a Seq-id
    /// (`header::holds_a_seqid`); the reader's error for a tag that NCBI's read cannot adjust
    /// (NCBI's `GetSeqIDs` of the Seq-align throws there), LOSAT's refusal of `filtered_mask`.
    /// With it, whether the record holds a UTF8String where a VisibleString is expected
    /// (`DeflineSet::utf8`; false where the volume filters no line, whose walk of the record does
    /// not read the strings).
    pub(super) fn line_read(&self, oid: u32) -> Result<(bool, bool), SeqDbError> {
        if !self.filters_lines() {
            let kept = self.with_header_record(oid, holds_a_seqid)?;
            return Ok((kept, false));
        }
        let set = self.with_header_record(oid, decode_full)?;
        let kept = self.filtered_mask(oid, &set)?.is_some_and(|mask| {
            mask.iter()
                .zip(&set.deflines)
                .any(|(&keep, defline)| keep && !defline.ids.is_empty())
        });
        Ok((kept, set.utf8))
    }

    /// `decoded`, the definition lines of `record` (one sequence's header), without those
    /// that the volume's `DeflineRule` removes; `decoded` itself without a rule.
    fn kept_lines<T: DeflineList>(
        &self,
        oid: u32,
        record: &[u8],
        path: &Path,
        decoded: Result<T, SeqDbError>,
    ) -> Result<T, SeqDbError> {
        if !self.filters_lines() {
            return decoded;
        }
        let mut decoded = decoded?;
        let keep = self.line_mask(oid, &decode_full(record, path)?)?;
        decoded.retain_lines(&keep);
        Ok(decoded)
    }

    /// The header of the sequence `oid` (volume-relative) in the bytes `record` (`header`): the
    /// definition lines of a record without parsed Seq-ids with their `taxid` (`deflines`).
    fn decode_header(
        &self,
        oid: u32,
        record: &[u8],
        path: &Path,
    ) -> Result<VolumeHeader, SeqDbError> {
        match decode_deflines_as(record, path, oid, true) {
            Err(error) if error.is_parsed_seqids() => self
                .decode_kept_set(oid, record, path)
                .map(VolumeHeader::Parsed),
            read => self
                .kept_lines(oid, record, path, read)
                .map(VolumeHeader::Ordinal),
        }
    }

    /// The header of the sequence `oid` (volume-relative) with every definition line it stores:
    /// no list's rule and no erasure of model lines, as NCBI reads the volume searched bare
    /// (`SeqDb::read_stored_header_of_oid`).
    pub fn stored_header(&self, oid: u32) -> Result<VolumeHeader, SeqDbError> {
        self.with_header_record(oid, |record, path| {
            match decode_deflines_as(record, path, oid, true) {
                Err(error) if error.is_parsed_seqids() => {
                    decode_full(record, path).map(VolumeHeader::Parsed)
                }
                read => read.map(VolumeHeader::Ordinal),
            }
        })
    }

    /// Whether the volume's reads remove a definition line of the sequence `oid`
    /// (volume-relative): its list's rule, the alias files' id lists or the erasure of model
    /// lines (`SeqDb::erases_a_line_of_oid`). The bytes of the removed lines and the accession
    /// guide's warnings are not looked at, as a read of the volume searched bare tests none; the
    /// reader's error for a header that cannot be read or a tag NCBI's read cannot adjust stays,
    /// and so does NCBI's exception of a lookup of the id test (`list_mask`).
    pub(super) fn erases_a_line(&self, oid: u32) -> Result<bool, SeqDbError> {
        if !self.filters_lines() {
            return Ok(false);
        }
        let set = self.with_header_record(oid, decode_full)?;
        let memb_kept = self.memb_mask(&set)?;
        let Some(mask) = self.list_mask(oid, &set, memb_kept.as_deref())? else {
            return Ok(true);
        };
        Ok(mask.iter().zip(&set.deflines).any(|(&keep, defline)| {
            !keep
                || (self.exclude_model && super::oid_mask::is_model_line(&defline.ids, &mut |_| {}))
        }))
    }

    /// LOSAT's refusal of the sequence `oid` (volume-relative) whose every definition line
    /// the seqid list removes (`filter::no_line_kept`).
    fn no_line_kept(&self, oid: u32) -> SeqDbError {
        SeqDbError::Unsupported {
            path: self.base.clone(),
            message: super::filter::no_line_kept(oid),
        }
    }

    /// `decode_full(record, path)` without the definition lines that the volume's
    /// `DeflineRule` removes.
    fn decode_kept_set(
        &self,
        oid: u32,
        record: &[u8],
        path: &Path,
    ) -> Result<crate::seqid::DeflineSet, SeqDbError> {
        let mut set = crate::seqid::DeflineSet::default();
        self.decode_kept_set_in(oid, record, path, &mut set)?;
        Ok(set)
    }

    /// `decode_kept_set` written over `set`, whose buffers are kept; whether the rule dropped a
    /// line of the record.
    fn decode_kept_set_in(
        &self,
        oid: u32,
        record: &[u8],
        path: &Path,
        set: &mut crate::seqid::DeflineSet,
    ) -> Result<bool, SeqDbError> {
        decode_full_in(record, path, set)?;
        let stored_lines = set.deflines.len();
        if self.filters_lines() {
            let keep = self.line_mask(oid, set)?;
            set.retain_lines(&keep);
        }
        Ok(set.deflines.len() != stored_lines)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2208-2218
    // NCBI decodes every member of the header for each report (`x_GetFilteredHeader`); the
    // reports read the title and the taxid of a record without parsed Seq-ids (its source,
    // `algorithm::common::database::HeaderReader::stored_source`) from this one decode.
    /// The definition lines of the sequence `oid` (volume-relative), with
    /// their ordinal ids still relative to the volume and their `taxid`.
    pub fn deflines(&self, oid: u32) -> Result<Vec<Defline>, SeqDbError> {
        self.deflines_as(oid, true)
    }

    /// `deflines` with every definition line of the stored header: the volume's `DeflineRule` is
    /// not applied (`SeqDb::has_parsed_seqids`, which checks the stored first header of each
    /// volume).
    pub(super) fn stored_deflines(&self, oid: u32) -> Result<Vec<Defline>, SeqDbError> {
        self.with_header_record(oid, |record, path| {
            decode_deflines_as(record, path, oid, false)
        })
    }

    /// `deflines`, with the `taxid` of each definition line decoded when `want_taxid`, as NCBI's
    /// `x_GetHdrAsn1` reads the record (`header::decode_deflines_as`; `defline_taxids` keeps the
    /// taxids of the same decode). The reads of the stored titles (`deflines`, `visit_deflines_as`)
    /// decode the taxid; only `stored_deflines` skips it.
    pub fn deflines_as(&self, oid: u32, want_taxid: bool) -> Result<Vec<Defline>, SeqDbError> {
        self.with_header_record(oid, |record, path| {
            self.kept_lines(
                oid,
                record,
                path,
                decode_deflines_as(record, path, oid, want_taxid),
            )
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:302-323
    /// The lengths that the index alone gives for the nucleotides `first..end` (volume-relative),
    /// one read of each offset array for the range: `GetSeqLengthExact`'s whole bytes times four
    /// plus 0 and plus 3 (the count that the last byte holds), so the length that
    /// [`Volume::nucleotide_length`] reads lies between them; `None` for a sequence whose length
    /// that read would not give (offsets out of order, past the end of the sequence file, or a
    /// length that may pass a signed 32-bit length), where `nucleotide_length` reports the error.
    /// No read of the sequence file. For the masks' walk of `-db_hard_mask` (`SeqDb::
    /// first_searched_masks`, review 2 S4 of M4a-2b), which reads the last byte only where the
    /// refusal differs within the bounds.
    pub fn nucleotide_length_bounds(
        &self,
        first: u32,
        end: u32,
    ) -> Result<Vec<Option<(usize, usize)>>, SeqDbError> {
        self.check_range(first, end)?;
        let seq = self.index.sequence_offsets(first, end)?;
        let amb = self.index.ambiguity_offsets(first, end)?;
        let file_len = self.seq_file()?.readable_len();
        Ok(seq
            .windows(2)
            .zip(&amb)
            .map(|(pair, &amb_start)| {
                let (start, amb_start, next) =
                    (u64::from(pair[0]), u64::from(amb_start), u64::from(pair[1]));
                if !(start < amb_start && amb_start <= next) || amb_start > file_len {
                    return None;
                }
                let high = packed_length(start, amb_start, 3);
                (high <= i32::MAX as u64)
                    .then(|| (packed_length(start, amb_start, 0) as usize, high as usize))
            })
            .collect())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:302-323
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1817-1841
    /// Signed 32-bit exact length of nucleotide `oid` (volume-relative; `named` its OID in the
    /// database) as NCBI's `GetSeqLengthExact` and `x_GetSequence` compute it:
    /// `whole_bytes * 4 + (file[start + whole_bytes] & 3)` where `whole_bytes = int(amb - start - 1)`.
    /// Refuses offsets whose `start + whole_bytes < 0` (e.g. `amb == 0`), where NCBI reads before
    /// the sequence file mapping and crashes.
    pub fn nucleotide_length_i32(&self, oid: u32, named: u32) -> Result<i32, SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let amb = self.index.ambiguity_offsets(oid, oid + 1)?;
        let seq = self.index.sequence_offsets(oid, oid + 1)?;
        let (start, amb_start) = (seq[0], amb[0]);
        let file = self.seq_file()?;
        if file.is_empty() {
            return Err(SeqDbError::NcbiCrashes {
                path: file.path().to_path_buf(),
                message: "the sequence file has no bytes, which NCBI maps as a null pointer"
                    .to_string(),
            });
        }
        let whole_bytes = (amb_start as i32)
            .wrapping_sub(start as i32)
            .wrapping_sub(1);
        let last_pos = i64::from(start) + i64::from(whole_bytes);
        if last_pos < 0 {
            return Err(SeqDbError::NcbiCrashes {
                path: file.path().to_path_buf(),
                message: format!(
                    "sequence {named}: ambiguity offset {amb_start} causes NCBI to read before \
                     the start of the sequence file mapping"
                ),
            });
        }
        // NCBI's `GetSeqLengthExact` / `x_GetSequence` (`seqdbvol.cpp:312-323,1817-1838`) does
        // pointer arithmetic on `start` and dereferences only `file[start + whole_bytes]`.
        let last = file.with_range(last_pos as u64, last_pos as u64 + 1, |b| b[0])?;
        if whole_bytes >= 0 {
            let length = (whole_bytes as u64) * 4 + u64::from(last & 3);
            if length > i32::MAX as u64 {
                return Err(SeqDbError::Format {
                    path: file.path().to_path_buf(),
                    message: format!(
                        "sequence {named}: the sequence is longer than a signed 32-bit length"
                    ),
                });
            }
            return Ok(length as i32);
        }
        Ok(whole_bytes
            .wrapping_mul(4)
            .wrapping_add(i32::from(last & 3)))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:302-323 (quoted
    // at `nucleotide_length_bounds`)
    /// The length in bases of nucleotide `oid` (volume-relative; `named` its OID in the database,
    /// which the messages name), from the last byte of its packed bases only.
    pub fn nucleotide_length(&self, oid: u32, named: u32) -> Result<usize, SeqDbError> {
        self.check_range(oid, oid.saturating_add(1))?;
        let amb = self.index.ambiguity_offsets(oid, oid + 1)?;
        let seq = self.index.sequence_offsets(oid, oid + 1)?;
        let (start, amb_start, next) = (u64::from(seq[0]), u64::from(amb[0]), u64::from(seq[1]));
        let file = self.seq_file()?;
        if !(start < amb_start && amb_start <= next) {
            return Err(SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!("sequence {named}: the ambiguity offset lies outside its data"),
            });
        }
        let last = file.with_range(amb_start - 1, amb_start, |last| last[0])?;
        let length = packed_length(start, amb_start, last);
        if length > i32::MAX as u64 {
            return Err(SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!(
                    "sequence {named}: the sequence is longer than a signed 32-bit length"
                ),
            });
        }
        Ok(length as usize)
    }

    /// The IUPAC letters `[from, to)` of nucleotide `oid` (volume-relative; `named` its OID in
    /// the database, which the messages name), as its stored data give them: the 2-bit bases of the range and the ambiguity runs that overlap
    /// it. Only the bytes of the range and the sequence's ambiguity data are read.
    pub fn nucleotide_iupac_range(
        &self,
        oid: u32,
        named: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        let ambiguity = self.nucleotide_ambiguity(oid, named)?;
        self.nucleotide_iupac_range_with(oid, named, from, to, &ambiguity)
    }

    /// The IUPAC letters of the whole nucleotide `oid` (volume-relative; `named` its OID in the
    /// database): [`Volume::nucleotide_iupac_range`] of `[0, length)` with the length that
    /// [`Volume::nucleotide_ambiguity`] reads, not read a second time.
    pub fn nucleotide_iupac(&self, oid: u32, named: u32) -> Result<Vec<u8>, SeqDbError> {
        let ambiguity = self.nucleotide_ambiguity(oid, named)?;
        self.nucleotide_iupac_range_with(oid, named, 0, ambiguity.length, &ambiguity)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1619-1647,2350-2386
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:652-662
    /// The length of nucleotide `oid` (volume-relative; `named` its OID in the database) and
    /// its decoded ambiguity runs, read once for any number of
    /// [`Volume::nucleotide_iupac_range_with`] calls on the same sequence.
    pub fn nucleotide_ambiguity(
        &self,
        oid: u32,
        named: u32,
    ) -> Result<NucleotideAmbiguity, SeqDbError> {
        let length_i32 = self.nucleotide_length_i32(oid, named)?;
        let file = self.seq_file()?;
        let seq = self.index.sequence_offsets(oid, oid + 1)?;
        if length_i32 < 1 {
            return Ok(NucleotideAmbiguity::new(0, u64::from(seq[0]), Vec::new()));
        }
        let length = length_i32 as usize;
        let amb = self.index.ambiguity_offsets(oid, oid + 1)?;
        let (packed_start, amb_start, next) = (seq[0], amb[0], seq[1]);
        if amb_start > next {
            return Err(SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: super::AMBIGUITY_DATA_ERROR.to_string(),
            });
        }
        let mut runs = Vec::new();
        if u64::from(amb_start) < u64::from(next) {
            let words: Vec<u32> =
                file.with_range(u64::from(amb_start), u64::from(next), |data| {
                    data.chunks_exact(4)
                        .map(|word| u32::from_be_bytes([word[0], word[1], word[2], word[3]]))
                        .collect()
                })?;
            runs = decode_ambiguity_runs(&words).map_err(|message| SeqDbError::Format {
                path: file.path().to_path_buf(),
                message: format!("sequence {named}: {message}"),
            })?;
        }
        Ok(NucleotideAmbiguity::new(
            length,
            u64::from(packed_start),
            runs,
        ))
    }

    /// [`Volume::nucleotide_iupac_range`] with the length and the ambiguity runs that
    /// [`Volume::nucleotide_ambiguity`] read for the same `oid`: one read of the range's packed
    /// bytes, and a binary search for the first run that can overlap the range.
    pub fn nucleotide_iupac_range_with(
        &self,
        oid: u32,
        named: u32,
        from: usize,
        to: usize,
        ambiguity: &NucleotideAmbiguity,
    ) -> Result<Vec<u8>, SeqDbError> {
        let length = ambiguity.length;
        if from > to || to > length {
            return Err(SeqDbError::Range(format!(
                "letters {from}..{to} are outside sequence {named}, which has {length}"
            )));
        }
        self.check_range(oid, oid.saturating_add(1))?;
        let first_byte = (from / 4) as u64;
        let end_byte = (to as u64).div_ceil(4);
        let packed = self.seq_file()?.read_range(
            ambiguity.packed_start + first_byte,
            ambiguity.packed_start + end_byte,
        )?;
        let mut letters: Vec<u8> = (from..to)
            .map(|at| {
                let byte = packed[at / 4 - first_byte as usize];
                b"ACGT"[usize::from((byte >> (6 - 2 * (at % 4))) & 3)]
            })
            .collect();
        ambiguity.overlay(from, to, &mut letters);
        Ok(letters)
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): calls `visit(oid, self.deflines(oid))`
    /// for the sequences `first..end` (volume-relative, non-empty) in order
    /// until `visit` returns false; returns whether every call returned true.
    /// The header offsets of the range are read at once and the records in
    /// pieces of at most `piece_bytes` (one longer record alone), each with one
    /// read; a record is decoded from its piece as `deflines` decodes it. A
    /// sequence whose record passes the end of the header file, and every
    /// sequence of a range or piece that cannot be read, goes through
    /// `deflines`, which gives its own result.
    #[cfg(test)]
    pub(super) fn visit_deflines(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        visit: &mut dyn FnMut(u32, Result<Vec<Defline>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_deflines_as(first, end, piece_bytes, true, visit)
    }

    /// `visit_deflines`, with the `taxid` of each definition line decoded when `want_taxid`
    /// (the results of `deflines_as`): one pass over the headers gives the titles and the taxids.
    pub(super) fn visit_deflines_as(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        want_taxid: bool,
        visit: &mut dyn FnMut(u32, Result<Vec<Defline>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_headers(
            first,
            end,
            piece_bytes,
            &|oid| self.deflines_as(oid, want_taxid),
            &|oid, record, path| {
                self.kept_lines(
                    oid,
                    record,
                    path,
                    decode_deflines_as(record, path, oid, want_taxid),
                )
            },
            visit,
        )
    }

    /// `visit_deflines` for `defline_taxids` (OF5b: the taxdb pre-pass of BLASTX 15/16 reads the
    /// taxid of every sequence's first definition line): the same pieces and fallbacks, with
    /// the decoder of `defline_taxids`.
    pub(super) fn visit_defline_taxids(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        visit: &mut dyn FnMut(u32, Result<Vec<i32>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_headers(
            first,
            end,
            piece_bytes,
            &|oid| self.defline_taxids(oid),
            &|oid, record, path| {
                self.kept_lines(oid, record, path, decode_defline_taxids(record, path, oid))
            },
            visit,
        )
    }

    /// `visit_defline_taxids` with the decoder of `defline_set` (`header::decode_full`), which
    /// reads the definition lines of any Seq-ids, parsed ones too: the taxid of every definition
    /// line, 0 for none (item T: the check of NCBI's taxonomy cache,
    /// `algorithm::common::database::TaxdbRun::cache_collision`).
    pub(super) fn visit_defline_set_taxids(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        visit: &mut dyn FnMut(u32, Result<Vec<i32>, SeqDbError>) -> bool,
    ) -> bool {
        let taxids = |set: crate::seqid::DeflineSet| -> Vec<i32> {
            set.deflines
                .iter()
                .map(|defline| defline.taxid.unwrap_or(0))
                .collect()
        };
        self.visit_headers(
            first,
            end,
            piece_bytes,
            &|oid| self.defline_set(oid).map(taxids),
            &|oid, record, path| self.decode_kept_set(oid, record, path).map(taxids),
            visit,
        )
    }

    /// The shared walk of `visit_deflines` and `visit_defline_taxids`: `single(oid)` is the
    /// one-sequence read that gives the result for a sequence that is not read in a piece,
    /// `decode(oid, record, path)` decodes a record read in a piece.
    fn visit_headers<T>(
        &self,
        first: u32,
        end: u32,
        piece_bytes: u64,
        single: &dyn Fn(u32) -> Result<T, SeqDbError>,
        decode: &dyn Fn(u32, &[u8], &Path) -> Result<T, SeqDbError>,
        visit: &mut dyn FnMut(u32, Result<T, SeqDbError>) -> bool,
    ) -> bool {
        // As `with_header_record`: under a setting that turns a check of NCBI's reader off, the
        // records read in a piece decode with `refused_read`'s errors (`single` reads through
        // `with_header_record`).
        let refusing;
        let decode = match self.read_verification.as_deref() {
            None => decode,
            Some(setting) => {
                refusing = move |oid: u32, record: &[u8], path: &Path| {
                    decode(oid, record, path).map_err(|error| refused_read(setting, error))
                };
                &refusing as &dyn Fn(u32, &[u8], &Path) -> Result<T, SeqDbError>
            }
        };
        let one_by_one =
            |oids: std::ops::Range<u32>,
             visit: &mut dyn FnMut(u32, Result<T, SeqDbError>) -> bool| {
                oids.into_iter().all(|oid| visit(oid, single(oid)))
            };
        if self.check_range(first, end).is_err() {
            return one_by_one(first..end, visit);
        }
        let (Ok(offsets), Ok(file)) = (self.index.header_offsets(first, end), self.hdr_file())
        else {
            return one_by_one(first..end, visit);
        };
        let hr_path = with_extension(&self.base, self.index.seqtype, "hr");
        let mut i = 0usize;
        let count = (end - first) as usize;
        while i < count {
            let start = u64::from(offsets[i]);
            if offsets[i + 1] < offsets[i] || u64::from(offsets[i + 1]) > file.len() {
                // Non-ascending or zero-page record: `deflines` gives its exact result.
                if !one_by_one(first + i as u32..first + i as u32 + 1, visit) {
                    return false;
                }
                i += 1;
                continue;
            }
            let mut stop = i + 1;
            while stop < count
                && offsets[stop + 1] >= offsets[stop]
                && u64::from(offsets[stop + 1]) <= file.len()
                && u64::from(offsets[stop + 1]) - start <= piece_bytes
            {
                stop += 1;
            }
            let piece_end = u64::from(offsets[stop]);
            let Ok(data) = file.read_range(start, piece_end) else {
                if !one_by_one(first + i as u32..first + stop as u32, visit) {
                    return false;
                }
                i = stop;
                continue;
            };
            for j in i..stop {
                let lo = (u64::from(offsets[j]) - start) as usize;
                let hi = (u64::from(offsets[j + 1]) - start) as usize;
                if !visit(
                    first + j as u32,
                    decode(first + j as u32, &data[lo..hi], &hr_path),
                ) {
                    return false;
                }
            }
            i = stop;
        }
        true
    }
}

/// EXPERIMENT (LOSAT_X_BNDBSMALL): sequence offsets of one volume from `first` on, kept by
/// `SeqDb::x_next_blocks` for `Volume::x_listed_block_end`.
#[derive(Debug, Default)]
#[cfg(any(
    test,
    all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )
))]
pub(super) struct XOffsetWindow {
    first: u32,
    offsets: Vec<u32>,
}

#[cfg(any(
    test,
    all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )
))]
impl XOffsetWindow {
    /// Whether the window holds the offset of `oid`.
    fn holds(&self, oid: u32) -> bool {
        oid >= self.first && ((oid - self.first) as usize) < self.offsets.len()
    }
}

/// Entries of the first window that `Volume::x_read_block_end` reads.
const X_BLOCK_END_WINDOW: usize = 1024;
/// The largest window of `Volume::x_read_block_end` (1 MiB of offsets).
const X_BLOCK_END_MAX_WINDOW: usize = 1 << 18;

/// `Volume::block_end` over the held sequence offsets (`lead`: `Volume::lead`).
pub(super) fn held_block_end(offsets: &[u32], first: u32, max_bytes: u64, lead: u64) -> u32 {
    let n = (offsets.len() - 1) as u32;
    if max_bytes < lead {
        return first + 1;
    }
    let first_start = u64::from(offsets[first as usize]);
    let limit = first_start.saturating_add(max_bytes - lead);
    let last_start = u64::from(offsets[n as usize]);
    if last_start >= first_start && last_start <= limit {
        return n;
    }
    let mut prev = first_start;
    for end in first + 1..=n {
        let offset = u64::from(offsets[end as usize]);
        if offset < prev || offset > limit {
            return (end - 1).max(first + 1);
        }
        prev = offset;
    }
    n
}

/// Where `Volume::read_block_in` puts the bytes of a block.
enum BlockBuffer<'a> {
    /// A new buffer (`read_block`).
    New,
    /// EXPERIMENT (LOSAT_X_BNMAPPED): `File::x_read_range_mapped`'s buffer.
    Mapped,
    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): the worker's buffer
    /// (`File::x_read_range_into`).
    Given(&'a mut Vec<u8>),
}

/// A contiguous range of sequences of one volume, held in memory.
///
/// Sequences are addressed by their position `i` in the block
/// (`0..self.len()`); the database ordinal id of position `i` is
/// `self.first_ordinal() + i`, its OID `first_oid + i`, by which the messages name it.
#[derive(Debug, Clone)]
pub struct SequenceBlock {
    seqtype: SeqType,
    first_ordinal: u32,
    first_oid: u32,
    path: PathBuf,
    data: Vec<u8>,
    data_start: u64,
    /// Start offsets of the block's sequences in the file, plus the end.
    seq_offsets: Vec<u32>,
    /// Ambiguity start offsets of the block's sequences (nucleotide).
    amb_offsets: Vec<u32>,
    /// Whether a multithreaded DB-split scan stopped early at a negative-length sequence
    /// (`CSeqDBImpl::x_FillSeqBuffer`, seqdbimpl.cpp:699-707).
    stopped_at_negative: bool,
    /// Per-sequence flag for nucleotide sequences with `amb > next` whose ambiguity read throws
    /// NCBI's `File error: could not get ambiguity data.` (`seqdbvol.cpp:2358-2363`); `None` on
    /// every normal block.
    amb_error: Option<Box<[bool]>>,
}

/// A run of one ambiguity code in a nucleotide sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AmbiguityRun {
    /// Offset of the first base of the run.
    pub offset: u32,
    /// Number of bases.
    pub length: u32,
    /// The ambiguity code (NCBI4na, 0-15).
    pub ncbi4na: u8,
}

/// The length of one nucleotide sequence and its ambiguity runs, read once by
/// [`Volume::nucleotide_ambiguity`] for any number of range reads
/// ([`Volume::nucleotide_iupac_range_with`]).
#[derive(Debug, Clone)]
pub struct NucleotideAmbiguity {
    length: usize,
    /// Offset in the sequence file of the first packed byte.
    packed_start: u64,
    /// The runs in stored order: a later run overwrites an earlier one.
    runs: Vec<AmbiguityRun>,
    /// Whether the runs are in increasing order of offset and do not overlap, as NCBI writes
    /// them: then the runs that overlap a range are found by binary search. Otherwise every run
    /// is applied in stored order.
    ordered: bool,
}

impl NucleotideAmbiguity {
    /// The ambiguity of a sequence without letters, which nothing reads (`seqmasks_io`'s record of
    /// a sequence that NCBI's `GetBioseq` gives no Bioseq for).
    pub(crate) fn none() -> Self {
        Self::new(0, 0, Vec::new())
    }

    fn new(length: usize, packed_start: u64, runs: Vec<AmbiguityRun>) -> Self {
        let ordered = runs.windows(2).all(|pair| {
            u64::from(pair[0].offset) + u64::from(pair[0].length) <= u64::from(pair[1].offset)
        });
        Self {
            length,
            packed_start,
            runs,
            ordered,
        }
    }

    /// The length of the sequence in bases.
    pub fn length(&self) -> usize {
        self.length
    }

    /// Writes the IUPAC letters of the runs that overlap `[from, to)` into `letters`, which
    /// holds the bases of that range.
    fn overlay(&self, from: usize, to: usize, letters: &mut [u8]) {
        overlay_ambiguity_runs(&self.runs, self.ordered, from, to, letters);
    }
}

/// Writes the IUPAC letters of `runs` over `letters` (the bases `[from, to)`), a later run over
/// an earlier one. With `ordered` runs (increasing offsets, no overlap) only those from the
/// first run that ends after `from` to the first that starts at or after `to` are visited.
fn overlay_ambiguity_runs(
    runs: &[AmbiguityRun],
    ordered: bool,
    from: usize,
    to: usize,
    letters: &mut [u8],
) {
    let first = if ordered {
        runs.partition_point(|run| (run.offset as usize + run.length as usize) <= from)
    } else {
        0
    };
    for run in &runs[first..] {
        if ordered && run.offset as usize >= to {
            break;
        }
        let run_from = (run.offset as usize).max(from);
        let run_to = (run.offset as usize + run.length as usize).min(to);
        for at in run_from..run_to.max(run_from) {
            letters[at - from] = NCBI4NA_TO_IUPAC[usize::from(run.ncbi4na & 0xF)];
        }
    }
}

/// EXPERIMENT (LOSAT_X_BNTBRANGE): the ambiguity codes of `runs` over the bases of `ranges`
/// (sorted, disjoint) in `seq` (the whole sequence's bases): the runs in stored order, with the
/// whole decode's tests (`nucleotide_blastna_in`: a run that ends before the sequence is skipped,
/// the first one that starts at or after its end stops the loop), each run written over the
/// ranges it overlaps. So a range's bytes are the whole decode's for ordered and unordered runs.
/// A run writes `letter[ncbi4na & 0xF]`: `X_NCBI4NA` (its code) for the decode that maps every
/// byte to BLASTNA after it, `NCBINA8_TO_BLASTNA8` (its BLASTNA letter) for `x_blastna_by_byte`.
fn x_overlay_runs_in_ranges(
    runs: &[AmbiguityRun],
    ranges: &[(usize, usize)],
    seq: &mut [u8],
    letter: &[u8; 16],
) {
    let region_end = seq.len() as i64;
    for run in runs {
        let position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= region_end {
            break;
        }
        let run_from = position.max(0) as usize;
        let run_to = (position + row_len).min(region_end) as usize;
        let first = ranges.partition_point(|r| r.1 <= run_from);
        for &(from, to) in &ranges[first..] {
            if from >= run_to {
                break;
            }
            let letter = letter[usize::from(run.ncbi4na & 0xF)];
            for at in run_from.max(from)..run_to.min(to) {
                seq[at] = letter;
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:590-607
const NCBINA8_TO_BLASTNA8: [u8; 16] = [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];

/// No NCBI counterpart: each NCBI4na code as itself (`x_overlay_runs_in_ranges`'s letters for a
/// decode that maps every byte to BLASTNA after the runs).
const X_NCBI4NA: [u8; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

/// EXPERIMENT (LOSAT_X_BNDECODE): the four bases of each stored 2-bit byte, from the high bits
/// down, as the BLASTNA letters NCBI's decode makes of them: `s_SeqDBMapNA2ToNA8` (cited above
/// `nucleotide_blastna`) turns code `p` into `1 << p`, `s_SeqDBMapNcbiNA8ToBlastNA8` that into
/// `NCBINA8_TO_BLASTNA8[1 << p]` (which is `p`: BLASTNA's A, C, G and T are 0 to 3).
const X_PACKED_TO_BLASTNA: [[u8; 4]; 256] = {
    let mut table = [[0u8; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut base = 0;
        while base < 4 {
            table[byte][base] = NCBINA8_TO_BLASTNA8[1 << ((byte >> (6 - 2 * base)) & 3)];
            base += 1;
        }
        byte += 1;
    }
    table
};

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNDECODE, in the strict set: on unless it is `0`)
/// a `-db` sequence's BLASTNA letters (`nucleotide_blastna`, its mapped form and
/// `x_nucleotide_blastna_ranges`) are made four bases per stored byte with their BLASTNA letters
/// (`x_blastna_by_byte`), and each ambiguity run writes its BLASTNA letter; `=0` makes each base's
/// NCBI8na value from its two bits, writes the runs' NCBI4na codes, then maps every byte to
/// BLASTNA. The letters are the same.
fn x_bn_decode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNDECODE"))
}

/// EXPERIMENT (LOSAT_X_BNDECODE): the BLASTNA letters of the bases `from .. to` of a stored
/// sequence (`packed`: four bases per byte, from the high bits down) into `seq[from .. to]`
/// (`seq` holds one byte per base of the sequence): the letters of each stored byte whose four
/// bases are all in the range at once, the bases of a byte the range cuts one by one.
fn x_blastna_bases(packed: &[u8], from: usize, to: usize, seq: &mut [u8]) {
    let head_end = from.next_multiple_of(4).min(to);
    for at in from..head_end {
        seq[at] = X_PACKED_TO_BLASTNA[usize::from(packed[at >> 2])][at & 3];
    }
    let body_end = head_end + (to - head_end) / 4 * 4;
    let bytes = &packed[head_end >> 2..body_end >> 2];
    for (bases, &byte) in seq[head_end..body_end].chunks_exact_mut(4).zip(bytes) {
        bases.copy_from_slice(&X_PACKED_TO_BLASTNA[usize::from(byte)]);
    }
    for at in body_end..to {
        seq[at] = X_PACKED_TO_BLASTNA[usize::from(packed[at >> 2])][at & 3];
    }
}

/// EXPERIMENT (LOSAT_X_BNDECODE): the BLASTNA letters of the bases of `ranges` (sorted, disjoint,
/// non-empty, within `seq`, which holds one byte per base) of a stored sequence with its ambiguity
/// `runs` in stored order: the bases four per stored byte (`x_blastna_bases`), then each run's
/// BLASTNA letter over the ranges it overlaps, with the whole decode's tests
/// (`x_overlay_runs_in_ranges`). The bytes of `=0`'s decode, which writes the bases' NCBI8na
/// values and the runs' NCBI4na codes and then maps each byte with `NCBINA8_TO_BLASTNA8`: that map
/// is per byte, so mapping each letter where it is written gives the same byte.
fn x_blastna_by_byte(
    packed: &[u8],
    runs: &[AmbiguityRun],
    ranges: &[(usize, usize)],
    seq: &mut [u8],
) {
    for &(from, to) in ranges {
        x_blastna_bases(packed, from, to, seq);
    }
    x_overlay_runs_in_ranges(runs, ranges, seq, &NCBINA8_TO_BLASTNA8);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:455-471
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:804-832
/// `s_SeqDBMapNA2ToNA8` and `s_SeqDBRebuildDNA_NA8` over the whole sequence: each base's NCBI8na
/// value `1 << p` from its two bits, then the ambiguity `runs` in stored order (a run that ends
/// before the sequence is skipped, the first that starts at or after its end stops the loop) write
/// their NCBI4na codes (`seq` holds one byte per base).
fn ncbi4na_by_base(packed: &[u8], runs: &[AmbiguityRun], seq: &mut [u8]) {
    for (position, base) in seq.iter_mut().enumerate() {
        let code = (packed[position >> 2] >> (6 - 2 * (position & 3))) & 3;
        *base = 1 << code;
    }
    let region_end = seq.len() as i64;
    for run in runs {
        let mut position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= region_end {
            break;
        }
        for _ in 0..row_len {
            if position >= 0 && position < region_end {
                seq[position as usize] = run.ncbi4na;
            }
            position += 1;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1716-1771
// The ambiguity words are read whole (a last partial word is not read, as
// `SequenceBlock::ambiguity_runs`), and a count that passes their end is an error where NCBI
// reads past its vector. A sequence of no base gives no base before its ambiguity data are read.
/// `SeqDB_UnpackAmbiguities`: the bases of a stored nucleotide (`packed`: four bases per byte,
/// the last byte holding the count of its bases; `ambig`: the stored ambiguity data) as NCBI4na
/// codes, one byte per base. The error is the message of a malformed sequence.
pub(crate) fn unpack_ambiguities(packed: &[u8], ambig: &[u8]) -> Result<Vec<u8>, &'static str> {
    let Some(&last) = packed.last() else {
        return Err("Error: packed sequence data is not valid.");
    };
    let base_length = (packed.len() - 1) * 4 + usize::from(last & 3);
    if base_length == 0 {
        return Ok(Vec::new());
    }
    let words: Vec<u32> = ambig
        .chunks_exact(4)
        .map(|word| u32::from_be_bytes([word[0], word[1], word[2], word[3]]))
        .collect();
    let runs = decode_ambiguity_runs(&words)?;
    let mut na8 = vec![0u8; base_length];
    ncbi4na_by_base(packed, &runs, &mut na8);
    Ok(na8)
}

/// The decode of `nucleotide_blastna` that `LOSAT_X_BNDECODE=0` runs: `ncbi4na_by_base`, then
/// every byte mapped with `NCBINA8_TO_BLASTNA8` (`seq` holds one byte per base).
fn x_blastna_by_base(packed: &[u8], runs: &[AmbiguityRun], seq: &mut [u8]) {
    ncbi4na_by_base(packed, runs, seq);
    for base in seq.iter_mut() {
        *base = NCBINA8_TO_BLASTNA8[usize::from(*base & 0xF)];
    }
}

/// The decode of `x_nucleotide_blastna_ranges` that `LOSAT_X_BNDECODE=0` runs on `ranges`
/// (sorted, disjoint, non-empty, within `seq`): `x_blastna_by_base`'s steps on each range, the runs
/// by `x_overlay_runs_in_ranges`.
fn x_blastna_by_base_ranges(
    packed: &[u8],
    runs: &[AmbiguityRun],
    ranges: &[(usize, usize)],
    seq: &mut [u8],
) {
    for &(from, to) in ranges {
        for (position, base) in seq[from..to].iter_mut().enumerate() {
            let position = from + position;
            let code = (packed[position >> 2] >> (6 - 2 * (position & 3))) & 3;
            *base = 1 << code;
        }
    }
    x_overlay_runs_in_ranges(runs, ranges, seq, &X_NCBI4NA);
    for &(from, to) in ranges {
        for base in &mut seq[from..to] {
            *base = NCBINA8_TO_BLASTNA8[usize::from(*base & 0xF)];
        }
    }
}

/// `NCBI4NA_TO_IUPACNA` (c++/src/algo/blast/core/blast_encoding.c:80-83).
const NCBI4NA_TO_IUPAC: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

/// The four bases of each stored 2-bit byte as two NCBI4na bytes, two bases a byte from the high
/// bits (`s_SeqDBMapNA2ToNA4Setup`: A 1, C 2, G 4, T 8).
const PACKED_TO_NCBI4NA: [[u8; 2]; 256] = {
    const NCBI4NA: [u8; 4] = [1, 2, 4, 8];
    let mut table = [[0u8; 2]; 256];
    let mut byte = 0;
    while byte < 256 {
        let [b0, b1, b2, b3] = [
            NCBI4NA[(byte >> 6) & 3],
            NCBI4NA[(byte >> 4) & 3],
            NCBI4NA[(byte >> 2) & 3],
            NCBI4NA[byte & 3],
        ];
        table[byte] = [(b0 << 4) | b1, (b2 << 4) | b3];
        byte += 1;
    }
    table
};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:728-786
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:803-844
/// Calls `set(position, ncbi4na)` for each base of the ambiguity runs `runs` (in stored order)
/// inside a sequence of `length` bases, as NCBI's rebuilds of the 4-bit (`s_SeqDBRebuildDNA_NA4`,
/// [`SequenceBlock::nucleotide_seq_data`]) and 8-bit letters (`s_SeqDBRebuildDNA_NA8`,
/// [`SequenceBlock::nucleotide_iupac`]) place them: the bases of a run outside it are skipped, and
/// the runs after one that starts at or past its end are not read.
fn for_each_ambiguous_base(runs: &[AmbiguityRun], length: usize, mut set: impl FnMut(usize, u8)) {
    let region_end = length as i64;
    for run in runs {
        let mut position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= region_end {
            break;
        }
        for _ in 0..row_len {
            if position >= 0 && position < region_end {
                set(position as usize, run.ncbi4na & 0xF);
            }
            position += 1;
        }
    }
}

/// The four bases of each stored 2-bit byte, from the high bits down, as IUPAC letters.
const PACKED_TO_IUPAC: [[u8; 4]; 256] = {
    let mut table = [[0u8; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut base = 0;
        while base < 4 {
            table[byte][base] = b"ACGT"[(byte >> (6 - 2 * base)) & 3];
            base += 1;
        }
        byte += 1;
    }
    table
};

/// Sentinel byte NCBI puts before and after a BLASTNA sequence.
pub const BLASTNA_SENTINEL: u8 = 15;

/// EXPERIMENT (LOSAT_X_BNTBRANGE): NCBI's `FENCE_SENTRY` (seqdbgeneral.hpp:1229, blast_util.h:364),
/// the byte around each range of a partly decoded sequence (`x_nucleotide_blastna_ranges`).
pub const FENCE_SENTRY: u8 = 201;

/// The ambiguity runs of the ambiguity words of one sequence, in stored order
/// (`SequenceBlock::ambiguity_runs`); the error is its message.
fn decode_ambiguity_runs(words: &[u32]) -> Result<Vec<AmbiguityRun>, &'static str> {
    let mut runs = Vec::new();
    let Some(&first) = words.first() else {
        return Ok(runs);
    };
    let new_format = first & 0x8000_0000 != 0;
    let amb_num = (first & 0x7FFF_FFFF) as usize;
    let words_per_entry = if new_format { 2 } else { 1 };
    runs.reserve(amb_num.min(words.len()) / words_per_entry);
    let mut index = 1;
    // In the new format the loop above advances `i` by two per entry, so
    // there `amb_num` counts 32-bit words, not entries.
    while index < amb_num + 1 {
        if index + words_per_entry > words.len() {
            return Err("the ambiguity data are shorter than their count");
        }
        let word = words[index];
        let ncbi4na = ((word >> 28) & 0xF) as u8;
        let (length, offset) = if new_format {
            (((word >> 16) & 0xFFF) + 1, words[index + 1])
        } else {
            (((word >> 24) & 0xF) + 1, word & 0x00FF_FFFF)
        };
        runs.push(AmbiguityRun {
            offset,
            length,
            ncbi4na,
        });
        index += words_per_entry;
    }
    Ok(runs)
}

impl SequenceBlock {
    /// EXPERIMENT (LOSAT_X_BNMAPPED): gives the block's bytes to `utils::x_mapped::X_BYTES` (a
    /// buffer of `File::x_read_range_mapped` of 1 MiB or more; smaller ones are dropped). No NCBI
    /// counterpart: only where the bytes go.
    pub fn x_recycle(mut self) {
        crate::utils::x_mapped::X_BYTES.give(&mut self.data);
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): gives the buffer of a block of
    /// `x_read_block_into` back for the worker's next block: into `kept` below 1 MiB, else to
    /// `utils::x_mapped::X_BYTES` (`x_mapped::X_BYTES.x_keep`). No NCBI counterpart: only where
    /// the bytes are kept.
    pub fn x_into_kept(self, kept: &mut Vec<u8>) {
        crate::utils::x_mapped::X_BYTES.x_keep(kept, self.data);
    }

    /// The bytes of a block of one nucleotide cut to its 2-bit
    /// data (`nucleotide_packed(0)`, which starts the block), in the block's buffer with its
    /// capacity. No NCBI counterpart: only where the bytes are kept.
    pub fn x_into_nucleotide_packed(mut self) -> Vec<u8> {
        self.expect(SeqType::Nucleotide);
        assert_eq!(self.len(), 1, "the packed bytes of a block of one sequence");
        assert_eq!(
            u64::from(self.seq_offsets[0]),
            self.data_start,
            "a nucleotide block starts with its first sequence"
        );
        let end = self.nucleotide_packed(0).len();
        self.data.truncate(end);
        self.data
    }

    /// Sequence type of the block.
    #[cfg(test)]
    pub fn seqtype(&self) -> SeqType {
        self.seqtype
    }

    /// Database ordinal id of the first sequence of the block.
    pub fn first_ordinal(&self) -> u32 {
        self.first_ordinal
    }

    /// Number of sequences in the block.
    pub fn len(&self) -> usize {
        self.seq_offsets.len() - 1
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:699-707
    /// Constructs an empty block marked `stopped_at_negative` when a multithreaded scan stopped
    /// before this block (for example at an excluded OID).
    pub(super) fn empty_stopped(
        seqtype: SeqType,
        first_ordinal: u32,
        first_oid: u32,
        path: PathBuf,
    ) -> Self {
        let (data, seq_offsets) = match seqtype {
            SeqType::Protein => (vec![0u8], vec![1u32]),
            SeqType::Nucleotide => (Vec::new(), vec![0u32]),
        };
        Self {
            seqtype,
            first_ordinal,
            first_oid,
            path,
            data,
            data_start: 0,
            seq_offsets,
            amb_offsets: Vec::new(),
            stopped_at_negative: true,
            amb_error: None,
        }
    }

    /// Whether a multithreaded DB-split scan stopped at a negative-length sequence before the
    /// end of the requested range (`CSeqDBImpl::x_FillSeqBuffer`, seqdbimpl.cpp:699-707).
    pub fn stopped_at_negative(&self) -> bool {
        self.stopped_at_negative
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2358-2363
    /// Checks whether nucleotide `i` had `amb > next` (`GetAmbStartEnd` false) and returns NCBI's
    /// `File error: could not get ambiguity data.` (`seqdbvol.cpp:2358-2363`).
    fn check_amb_valid(&self, i: usize) -> Result<(), SeqDbError> {
        if self
            .amb_error
            .as_deref()
            .and_then(|flags| flags.get(i))
            .copied()
            .unwrap_or(false)
        {
            return Err(SeqDbError::Format {
                path: self.path.clone(),
                message: super::AMBIGUITY_DATA_ERROR.to_string(),
            });
        }
        Ok(())
    }

    fn slice(&self, start: u64, end: u64) -> &[u8] {
        let lo = (start - self.data_start) as usize;
        let hi = (end - self.data_start) as usize;
        &self.data[lo..hi]
    }

    fn format(&self, i: usize, message: &str) -> SeqDbError {
        SeqDbError::Format {
            path: self.path.clone(),
            message: format!(
                "sequence {}: {message}",
                u64::from(self.first_oid) + i as u64
            ),
        }
    }

    fn expect(&self, seqtype: SeqType) {
        assert!(
            self.seqtype == seqtype,
            "sequence block of the other sequence type"
        );
    }

    /// The residues of protein `i` in NCBIstdaa, one byte each.
    pub fn protein(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Protein);
        self.slice(
            u64::from(self.seq_offsets[i]),
            u64::from(self.seq_offsets[i + 1]) - 1,
        )
    }

    /// The residues of protein `i`, to be changed in place: the copy that a fetch after the
    /// preliminary stage of a search with `-db_hard_mask` masks
    /// (`algorithm::common::db_mask_search::mask_protein_letters`).
    pub fn protein_mut(&mut self, i: usize) -> &mut [u8] {
        self.expect(SeqType::Protein);
        let lo = (u64::from(self.seq_offsets[i]) - self.data_start) as usize;
        let hi = (u64::from(self.seq_offsets[i + 1]) - 1 - self.data_start) as usize;
        &mut self.data[lo..hi]
    }

    /// Protein `i` with the sentinel byte (NUL in a well-formed volume)
    /// before and after its residues, as NCBI's buffer has them.
    pub fn protein_with_sentinels(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Protein);
        self.slice(
            u64::from(self.seq_offsets[i]) - 1,
            u64::from(self.seq_offsets[i + 1]),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1791-1803
    // NCBI's search reads these bytes in place, and so does LOSAT's (a block
    // is a private copy of them in the layout of the file).
    /// The bytes of a protein block as read: the sentinel byte before the
    /// first sequence, then every sequence followed by its sentinel byte
    /// (`protein_with_sentinels(i)` is a window of them).
    pub fn protein_bytes(&self) -> &[u8] {
        self.expect(SeqType::Protein);
        &self.data
    }

    /// `protein_bytes`, writable: the tests write bytes that a file could hold.
    #[cfg(test)]
    pub fn protein_bytes_mut(&mut self) -> &mut [u8] {
        self.expect(SeqType::Protein);
        &mut self.data
    }

    /// The 2-bit data of nucleotide `i` as stored: four bases per byte from
    /// the high bits down, the last byte holding the remaining bases and,
    /// in its two low bits, their count.
    pub fn nucleotide_packed(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Nucleotide);
        self.slice(
            u64::from(self.seq_offsets[i]),
            u64::from(self.amb_offsets[i]),
        )
    }

    /// Length of nucleotide `i` in bases. A length above NCBI's signed
    /// 32-bit `int` is a format error.
    pub fn nucleotide_len(&self, i: usize) -> Result<usize, SeqDbError> {
        self.nucleotide_packed_len(i).map(|(_, length)| length)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1833-1841
    // NCBI's `x_GetSequence` returns the data (`*buffer`) and the length from one pair of
    // offsets; a database search reads both for every sequence (`blastn` `block_subject`).
    /// `nucleotide_packed(i)` and `nucleotide_len(i)`, from one read of the block's offsets.
    pub fn nucleotide_packed_len(&self, i: usize) -> Result<(&[u8], usize), SeqDbError> {
        let packed = self.nucleotide_packed(i);
        let whole_bytes = packed.len() - 1;
        let remainder = u64::from(packed[whole_bytes] & 3);
        let length = whole_bytes as u64 * 4 + remainder;
        if length > i32::MAX as u64 {
            return Err(self.format(i, "the sequence is longer than a signed 32-bit length"));
        }
        Ok((packed, length as usize))
    }

    /// The ambiguity data of nucleotide `i` as stored (empty when the
    /// sequence has no ambiguous base).
    pub fn ambiguity_bytes(&self, i: usize) -> &[u8] {
        self.expect(SeqType::Nucleotide);
        self.slice(
            u64::from(self.amb_offsets[i]),
            u64::from(self.seq_offsets[i + 1]),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2364-2384
    fn ambiguity_words(&self, i: usize) -> Vec<u32> {
        self.ambiguity_bytes(i)
            .chunks_exact(4)
            .map(|word| u32::from_be_bytes([word[0], word[1], word[2], word[3]]))
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:642-713
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:810-832
    // The entry count comes from the data; a count that passes the end of the
    // ambiguity data, where NCBI reads past its vector, is a format error.
    /// The ambiguity runs of nucleotide `i`, in stored order.
    pub fn ambiguity_runs(&self, i: usize) -> Result<Vec<AmbiguityRun>, SeqDbError> {
        self.check_amb_valid(i)?;
        decode_ambiguity_runs(&self.ambiguity_words(i)).map_err(|message| self.format(i, message))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1619-1707
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:834-842
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:619-625
    // The whole sequence is the range; no mask is applied (LOSAT reads no
    // mask data). `s_SeqDBMapNA2ToNA8` (seqdbvol.cpp:455-471) turns each
    // 2-bit code `p` into the NCBI8na value `1 << p`.
    /// Nucleotide `i` in BLASTNA, one byte per base with the ambiguity codes
    /// restored, and the sentinel byte 15 before and after the bases (the
    /// buffer NCBI's traceback gets): `len + 2` bytes.
    pub fn nucleotide_blastna(&self, i: usize) -> Result<Vec<u8>, SeqDbError> {
        self.nucleotide_blastna_in(i, false)
    }

    // No NCBI counterpart: the formatter's letters of a database sequence (`CSeqVector`) as
    // `nucleotide_blastna` gives them, without the BLASTNA buffer in between.
    /// The IUPAC letters of nucleotide `i`, the letters of [`Self::nucleotide_blastna`] without
    /// its sentinels (the same runs in the same order, the same errors), made from each stored
    /// byte's four bases at once.
    pub fn nucleotide_iupac(&self, i: usize) -> Result<Vec<u8>, SeqDbError> {
        let packed = self.nucleotide_packed(i);
        let (base_length, runs) = self.nucleotide_stored(i)?;
        let mut letters = vec![0u8; base_length];
        for (bases, &byte) in letters.chunks_mut(4).zip(packed) {
            bases.copy_from_slice(&PACKED_TO_IUPAC[usize::from(byte)][..bases.len()]);
        }
        for_each_ambiguous_base(&runs, base_length, |position, ncbi4na| {
            letters[position] = NCBI4NA_TO_IUPAC[usize::from(ncbi4na)];
        });
        Ok(letters)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1380-1424
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2350-2385
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:896-945
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:403-442
    /// The `Seq-data` of nucleotide `i` in the Bioseq that NCBI's `GetBioseq` makes: with fewer
    /// than 4 bytes of ambiguity data (`x_GetAmbChar` reads whole 4-byte words: none), the stored
    /// 2-bit bytes of its bases with the count bits of the last byte cleared and its bits after
    /// the last base as stored (`ncbi2na`, `true`); with a word (even a count of no run), NCBI4na
    /// (`false`), two bases a
    /// byte from the high bits, made from the stored bytes whose bits past the length are
    /// cleared first (`s_SeqDBMapNA2ToNA4`: the half byte after the last base of an odd length
    /// is the A of those bits), with each ambiguity run written over it
    /// (`s_SeqDBRebuildDNA_NA4`, the runs as [`Self::nucleotide_iupac`] places them).
    pub fn nucleotide_seq_data(&self, i: usize) -> Result<(bool, Vec<u8>), SeqDbError> {
        let (packed, length) = self.nucleotide_packed_len(i)?;
        self.check_amb_valid(i)?;
        let whole = length / 4;
        let remainder = length % 4;
        if self.ambiguity_bytes(i).len() < 4 {
            // The last stored byte of a length not a multiple of 4 with its count bits cleared:
            // the bits after the last base stay as stored.
            let mut data = Vec::with_capacity(whole + 1);
            data.extend_from_slice(&packed[..whole]);
            if remainder != 0 {
                data.push(packed[whole] & 0xFC);
            }
            return Ok((true, data));
        }
        // The bases of the last stored byte up to the length, the bits after them cleared.
        let last = (remainder != 0).then(|| packed[whole] & (0xFFu8 << (8 - 2 * remainder)));
        let runs = self.ambiguity_runs(i)?;
        let mut data = Vec::with_capacity(length.div_ceil(2));
        for &byte in &packed[..whole] {
            data.extend_from_slice(&PACKED_TO_NCBI4NA[usize::from(byte)]);
        }
        if let Some(byte) = last {
            let pair = PACKED_TO_NCBI4NA[usize::from(byte)];
            data.push(pair[0]);
            if remainder > 2 {
                data.push(pair[1]);
            }
        }
        for_each_ambiguous_base(&runs, length, |position, ncbi4na| {
            let byte = &mut data[position / 2];
            *byte = if position % 2 == 0 {
                (*byte & 0x0F) | (ncbi4na << 4)
            } else {
                (*byte & 0xF0) | ncbi4na
            };
        });
        Ok((false, data))
    }

    // NCBI: the checks of `x_GetSequence` cited above `nucleotide_blastna`
    // (seqdbvol.cpp:1619-1707: `base_length < 1` is "Error: could not get sequence or range.").
    /// The length of nucleotide `i` and its ambiguity runs in stored order, as
    /// [`Self::nucleotide_iupac`] reads them before it makes the letters (the same errors in the
    /// same order: the length, a sequence without bases, the ambiguity data), for a store of the
    /// stored data (`PackedNucleotide::from_stored`) that gives the same letters.
    pub fn nucleotide_stored(&self, i: usize) -> Result<(usize, Vec<AmbiguityRun>), SeqDbError> {
        let base_length = self.nucleotide_len(i)?;
        if base_length < 1 {
            return Err(self.format(i, "Error: could not get sequence or range."));
        }
        Ok((base_length, self.ambiguity_runs(i)?))
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `nucleotide_blastna` in a buffer taken from or sized
    /// by `utils::x_mapped::X_BYTES.reserve` (a mapped buffer from 1 MiB on, given back with
    /// `x_mapped::X_BYTES.give`). No NCBI counterpart: only where the bytes are kept.
    pub fn x_nucleotide_blastna_mapped(&self, i: usize) -> Result<Vec<u8>, SeqDbError> {
        self.nucleotide_blastna_in(i, true)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1565-1595
    // The ambiguity runs are applied in stored order, each to the ranges it overlaps, with the
    // whole decode's tests (`nucleotide_blastna_in`), so every byte of a range is the whole
    // decode's (NCBI's slices stop at the first run past each slice).
    /// EXPERIMENT (LOSAT_X_BNTBRANGE): `nucleotide_blastna` of the bases of `ranges` only
    /// (`[from, to)`, sorted and disjoint, as NCBI's partial-fetch ranges are; a range past the
    /// end is cut there), with `FENCE_SENTRY` before and after
    /// each range that does not reach an end of the sequence (`GetAmbigPartialSeq`): `len + 2`
    /// bytes, zero outside the ranges. The buffer is zeroed by the allocator; from
    /// `x_mapped::X_MAPPED_FROM` on its capacity is above glibc's mmap ceiling, so glibc maps it
    /// (unless an arena has a free chunk that large) and writes none of it: the pages outside the
    /// ranges are not resident. Other allocators may zero the whole buffer (Wasm's dlmalloc does).
    pub fn x_nucleotide_blastna_ranges(
        &self,
        i: usize,
        ranges: &[(usize, usize)],
    ) -> Result<Vec<u8>, SeqDbError> {
        debug_assert!(ranges.windows(2).all(|pair| pair[0].1 <= pair[1].0));
        let packed = self.nucleotide_packed(i);
        let base_length = self.nucleotide_len(i)?;
        if base_length < 1 {
            return Err(self.format(i, "Error: could not get sequence or range."));
        }
        let mut buffer = vec![0u8; crate::utils::x_mapped::x_mapped_capacity(base_length + 2, 1)];
        buffer.truncate(base_length + 2);
        {
            let seq = &mut buffer[1..=base_length];
            let ranges = ranges
                .iter()
                .map(|&(from, to)| (from.min(base_length), to.min(base_length)));
            for (from, to) in ranges.clone() {
                if from > 0 {
                    seq[from - 1] = FENCE_SENTRY;
                }
                if to < base_length {
                    seq[to] = FENCE_SENTRY;
                }
            }
            let ranges: Vec<(usize, usize)> = ranges.filter(|r| r.0 < r.1).collect();
            let runs = self.ambiguity_runs(i)?;
            if x_bn_decode() {
                x_blastna_by_byte(packed, &runs, &ranges, seq);
            } else {
                x_blastna_by_base_ranges(packed, &runs, &ranges, seq);
            }
        }
        buffer[0] = BLASTNA_SENTINEL;
        buffer[base_length + 1] = BLASTNA_SENTINEL;
        Ok(buffer)
    }

    fn nucleotide_blastna_in(&self, i: usize, mapped: bool) -> Result<Vec<u8>, SeqDbError> {
        let packed = self.nucleotide_packed(i);
        let base_length = self.nucleotide_len(i)?;
        if base_length < 1 {
            return Err(self.format(i, "Error: could not get sequence or range."));
        }
        let mut buffer = if mapped {
            let mut buffer = Vec::new();
            crate::utils::x_mapped::X_BYTES.reserve(&mut buffer, base_length + 2);
            buffer.resize(base_length + 2, 0u8);
            buffer
        } else {
            vec![0u8; base_length + 2]
        };
        {
            let seq = &mut buffer[1..=base_length];
            let runs = self.ambiguity_runs(i)?;
            if x_bn_decode() {
                x_blastna_by_byte(packed, &runs, &[(0, base_length)], seq);
            } else {
                x_blastna_by_base(packed, &runs, seq);
            }
        }
        buffer[0] = BLASTNA_SENTINEL;
        buffer[base_length + 1] = BLASTNA_SENTINEL;
        Ok(buffer)
    }
}

#[cfg(test)]
#[path = "volume_tests.rs"]
mod tests;
