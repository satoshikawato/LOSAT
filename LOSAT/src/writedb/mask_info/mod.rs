//! The mask files of `makeblastdb -mask_data`: `Blast-db-mask-info` objects
//! (blastdb.asn:71-79) in ASN.1 text or binary (BER), as dustmasker,
//! segmasker and windowmasker write them with `-outfmt maskinfo_asn1_text`
//! and `-outfmt maskinfo_asn1_bin`.
//!
//! Nothing on the command line calls this yet (`losat makeblastdb` refuses `-mask_data`).
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
//! What is refused (`MaskFileError::Unsupported`, "mask data that LOSAT does not read"): a
//! `Seq-loc` that is not an interval, packed interval or point list (`whole`, `pnt`, a mix, a
//! feature: NCBI reads them and LOSAT does not), strings with bytes outside 0x20-0x7e that NCBI
//! repairs and logs (`fix_visible_chars`; a binary Seq-id has no such byte), and a `UTF8String`
//! tag where a `VisibleString` belongs (objistrasnb.cpp:276-312).

mod ber;
mod guess;
mod reader;
mod schema;
#[cfg(test)]
mod seq_id_tests;
#[cfg(test)]
mod tests;
mod text;

pub(crate) use reader::{MaskFileError, MaskFileFormat, MaskFileReader};

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
/// `CSeq_loc::GetId()` returns and the intervals that `CSeq_loc_CI` visits,
/// in its order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MaskLoc {
    /// `None` when the Seq-loc has no Seq-id or several different ones;
    /// makeblastdb then stops with "Cannot get masked range Seq-id".
    pub(crate) id: Option<ParsedSeqId>,
    pub(crate) intervals: Vec<MaskInterval>,
}

/// A Seq-interval of a mask, `from` and `to` inclusive as stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MaskInterval {
    pub(crate) from: u32,
    pub(crate) to: u32,
    /// The `Na-strand` value (seqloc.asn) when the interval states one.
    pub(crate) strand: Option<u8>,
}
