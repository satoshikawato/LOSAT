//! The mask files of `makeblastdb -mask_data`: `Blast-db-mask-info` objects
//! (blastdb.asn:71-79) in ASN.1 text or binary (BER), as dustmasker,
//! segmasker and windowmasker write them with `-outfmt maskinfo_asn1_text`
//! and `-outfmt maskinfo_asn1_bin`.
//!
//! Line NP, M4 preparation: nothing calls this yet (`LOSAT makedb` still
//! refuses `-mask_data`). `writedb::mask_data` turns these objects into the
//! mask column of a database.

/// One `Blast-db-mask-info` object as makeblastdb reads it
/// (makeblastdb.cpp:950-990).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskInfo {
    /// `algo-id`: read, not used by makeblastdb.
    pub algo_id: i32,
    /// `algo-program`, a `Blast-filter-program` value (blastdb.asn:51-59).
    pub algo_program: i32,
    /// `algo-options`.
    pub algo_options: Vec<u8>,
    /// `masks`.
    pub masks: MaskList,
}

/// A `Blast-mask-list` (blastdb.asn:61-67): the one inside a
/// `Blast-db-mask-info`, or a continuation that makeblastdb reads after
/// `more TRUE`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskList {
    pub masks: Vec<MaskLoc>,
    pub more: bool,
}

/// One Seq-loc of `masks`, reduced to what makeblastdb uses: the Seq-id that
/// `CSeq_loc::GetId()` returns and the intervals that `CSeq_loc_CI` visits,
/// in its order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskLoc {
    /// `None` when the Seq-loc has no Seq-id or several different ones;
    /// makeblastdb then stops with "Cannot get masked range Seq-id".
    pub id: Option<MaskSeqId>,
    pub intervals: Vec<MaskInterval>,
}

/// A Seq-interval of a mask, `from` and `to` inclusive as stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaskInterval {
    pub from: u32,
    pub to: u32,
    /// The `Na-strand` value (seqloc.asn) when the interval states one.
    pub strand: Option<u8>,
}

/// The Seq-id of a mask. Only `local id N` can name a sequence of a database
/// that LOSAT writes: without `-parse_seqids` the FASTA reader names the N-th
/// sequence of a file `lcl|N`. Every other Seq-id is kept only to compare it
/// with the others of its Seq-loc (`CSerialObject::Equals`).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MaskSeqId {
    /// `local id N` (the `id` choice of Object-id).
    LocalId(i32),
    /// Any other Seq-id, as a canonical encoding of its value: two are equal
    /// exactly when NCBI's `Equals` finds them equal.
    Other(Vec<u8>),
}
