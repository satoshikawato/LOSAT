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

mod ber;
mod guess;
mod reader;
mod schema;
#[cfg(test)]
mod tests;
mod text;

pub(crate) use reader::{MaskFileError, MaskFileFormat, MaskFileReader};

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
    pub(crate) id: Option<MaskSeqId>,
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

/// The Seq-id of a mask. Only `local id N` can name a sequence of a database
/// that LOSAT writes: without `-parse_seqids` the FASTA reader names the N-th
/// sequence of a file `lcl|N`. Every other Seq-id is kept only to compare it
/// with the others of its Seq-loc (`CSerialObject::Equals`).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) enum MaskSeqId {
    /// `local id N` (the `id` choice of Object-id).
    LocalId(i32),
    /// Any other Seq-id, as a canonical encoding of its value: two are equal
    /// exactly when NCBI's `Equals` finds them equal.
    Other(Vec<u8>),
}
