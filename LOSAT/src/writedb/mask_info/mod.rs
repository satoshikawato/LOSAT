//! The mask files of `makeblastdb -mask_data`: `Blast-db-mask-info` objects
//! (blastdb.asn:71-79) in ASN.1 text or binary (BER), as dustmasker,
//! segmasker and windowmasker write them with `-outfmt maskinfo_asn1_text`
//! and `-outfmt maskinfo_asn1_bin`.
//!
//! `losat makeblastdb -mask_data` reads them (`writedb::mask_data::process`), and
//! `writedb::mask_data` turns these objects into the mask column of a database.
//!
//! [`MaskFileReader`] reads such a file as `makeblastdb` does: the format is guessed
//! once per file (`reader.rs`, `guess.rs`), then the objects follow one another in ASN.1 text
//! (`text.rs`) or binary (`ber.rs`); `schema.rs` has the types.
//!
//! What a mask names. Every Seq-id of a file is a [`ParsedSeqId`] of any alternative of the
//! `Seq-id` CHOICE (`local`, `gibbsq`, `gibbmt`, `giim`, `genbank`, `embl`, `pir`, `swissprot`,
//! `patent`, `other`, `general`, `gi`, `ddbj`, `prf`, `pdb`, `tpg`, `tpe`, `tpd`, `gpipe`,
//! `named-annot-track`), built by one builder (`reader.rs`) for text and for binary, as
//! `CSeq_loc::GetId` hands `CMaskedRangeSet::Insert` the Seq-id of a mask's Seq-loc
//! (masked_range_set.cpp:86-106; makeblastdb.cpp:976-995). A Seq-loc of several different Seq-ids
//! has no Seq-id (`CSeq_loc::GetId` fails, "Cannot get masked range Seq-id":
//! makeblastdb.cpp:976-995); the reader reports `id: None` and `process` raises that error as NCBI does. Whether
//! a mask reaches a sequence is decided by `seqid::mask_handle` and `mask_data::range_set`, not
//! here.
//!
//! Every `Seq-loc` alternative is read as `CSeq_loc_CI` and `CSeq_loc::GetId` see it ([`MaskLoc`]):
//! `null`, `empty`, `whole`, `int`, `packed-int`, `pnt`, `packed-pnt`, `mix` and `equiv` (nested),
//! `bond`; a `feat` has no Seq-id ("Cannot get masked range Seq-id", as `null` alone).
//!
//! A `VisibleString` is read as NCBI reads it (`eFNP_ReplaceAndWarn`, `report::serial::visible`):
//! each byte outside 0x20-0x7e is stored as `#`, and the reader hands out NCBI's diagnostic for it
//! with the stack trace of its frames ([`MaskDiagnostic`]: "Bad char [0x..] in string at line N"
//! with the text's chunk, "at byte N" with the binary string); a `UTF8String` tag where a
//! `VisibleString` belongs is read with NCBI's warning (objistrasnb.cpp:276-312), which NCBI posts
//! once per process. `process_mask_data` posts them in the order read.

mod ber;
mod guess;
mod reader;
mod schema;
#[cfg(test)]
mod seq_id_tests;
#[cfg(test)]
mod tests;
mod text;

pub(crate) use reader::{MaskDiagnostic, MaskFileError, MaskFileFormat, MaskFileReader};

use crate::seqid::ParsedSeqId;

/// One `Blast-db-mask-info` object as makeblastdb reads it
/// (makeblastdb.cpp:950-990).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MaskInfo {
    /// `algo-id`: read, not used by makeblastdb.
    pub(crate) algo_id: i32,
    /// `algo-program`, a `Blast-filter-program` value (blastdb.asn:51-59).
    pub(crate) algo_program: i32,
    /// `algo-options`.
    pub(crate) algo_options: Vec<u8>,
    /// `masks`.
    pub(crate) masks: MaskList,
}

/// A `Blast-mask-list` (blastdb.asn:61-67): the one inside a
/// `Blast-db-mask-info`, or a continuation that makeblastdb reads after
/// `more TRUE`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MaskList {
    pub(crate) masks: Vec<MaskLoc>,
    pub(crate) more: bool,
}

/// One Seq-loc of `masks`, reduced to what makeblastdb uses: the Seq-id that
/// `CSeq_loc::GetId()` returns and the ranges that `CSeq_loc_CI` visits,
/// in its order (`eEmpty_Allow`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MaskLoc {
    /// `None` when the Seq-loc has no Seq-id, several different ones or a `feat`;
    /// makeblastdb then stops with "Cannot get masked range Seq-id".
    pub(crate) id: Option<ParsedSeqId>,
    pub(crate) parts: Vec<MaskPart>,
}

#[cfg(test)]
impl MaskLoc {
    /// The intervals among the parts.
    pub(crate) fn intervals(&self) -> Vec<MaskInterval> {
        self.parts
            .iter()
            .filter_map(|p| match p {
                MaskPart::Interval(iv) => Some(*iv),
                _ => None,
            })
            .collect()
    }
}

/// One range of a Seq-loc as `CSeq_loc_CI` visits it. Every part but a `null` one has the Seq-loc's
/// Seq-id (`MaskLoc::id`: they are all `Equals` to it when it is set).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MaskPart {
    /// An interval (`int`, an element of `packed-int`), a point (`pnt`, a point of `packed-pnt`,
    /// a point of a `bond`: `from == to`), or a `whole` Seq-loc (`from` 0, `to` 4294967294: the
    /// same `CRange`, `TRange::GetWhole()`, and no strand).
    Interval(MaskInterval),
    /// An `empty` Seq-loc: skipped where the ranges are written, the empty range inside
    /// `CSeq_loc::Add`.
    Empty,
    /// A `null` Seq-loc: no Seq-id.
    Null,
}

/// A Seq-interval of a mask, `from` and `to` inclusive as stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MaskInterval {
    pub(crate) from: u32,
    pub(crate) to: u32,
    /// The `Na-strand` value (seqloc.asn) when the interval states one.
    pub(crate) strand: Option<u8>,
}
