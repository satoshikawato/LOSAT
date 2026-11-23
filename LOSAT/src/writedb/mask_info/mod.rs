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
//! or binary, read by the reader of the type table (`report::serial::asn_types`, whose table
//! holds these types: `Blast-db-mask-info` and what it reaches).
//!
//! What a mask names. Every Seq-id of a file is a [`ParsedSeqId`] of any alternative of the
//! `Seq-id` CHOICE (`local`, `gibbsq`, `gibbmt`, `giim`, `genbank`, `embl`, `pir`, `swissprot`,
//! `patent`, `other`, `general`, `gi`, `ddbj`, `prf`, `pdb`, `tpg`, `tpe`, `tpd`, `gpipe`,
//! `named-annot-track`), built by the one Seq-id reader (`asn_types::seq_id`) for text and for
//! binary, as
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
//! A string's bytes outside 0x20-0x7e are read as `#` with NCBI's diagnostic, and a `UTF8String`
//! tag where a `VisibleString` belongs with NCBI's warning (once per run): the reader hands them
//! out (`asn_types::Post`, `MaskFileReader::take_posts`) and `process` posts them. A file that NCBI's
//! serial reader cannot read ends with NCBI's exception (`MaskFileError::Ncbi`).
//!
//! NCBI's pointer syntax is read where NCBI reads a `CRef<>` that is empty at that moment: an
//! `OPTIONAL` member held by `CRef<>` (`fuzz-from`, `fuzz-to`, `fuzz`, `Seq-bond.b`,
//! `PDB-seq-id.rel`) and each element of `masks`, `mix`, `equiv` and `packed-int`. `NULL` leaves an
//! optional member unset; `@n` is NCBI's error ("invalid object index: NO_COLLECT defined");
//! `:Type value` naming the declared type is the value, read in a frame of its own.
//!
//! What is refused (`MaskFileError::Unsupported`, "mask data that LOSAT does not read"): objects
//! nested more than 400 levels deep; a string whose length passes the end of the file by an amount
//! that makeblastdb may or may not be able to allocate on the host (1 GiB to 128 TiB); `:Type`
//! naming a type other than the declared one, which makeblastdb reads or not by the types its
//! process has met before; and an element of a list stated as `NULL` that makeblastdb dereferences
//! (`NullElement`), whose `CCoreException` report holds the run's memory addresses and the path of
//! the binary. That last one is refused where makeblastdb dereferences it (`process`), after the
//! masks before it and after any error of the reader in the rest of its object.

mod guess;
mod reader;
#[cfg(test)]
mod seq_id_tests;
#[cfg(test)]
mod tests;

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
    /// The masks that makeblastdb's loop takes before `null_element`, all of them without one.
    pub(crate) masks: Vec<MaskLoc>,
    pub(crate) more: bool,
    /// The element stated as `NULL` that makeblastdb dereferences after `masks`, in the `GetId` of
    /// the next mask: that mask itself, or an element inside it that `CSeq_loc::CheckId` reaches.
    /// LOSAT refuses it there ([`NullElement::refusal`]).
    pub(crate) null_element: Option<NullElement>,
}

/// An element of a list of `CRef`s (`masks`, `mix`, `equiv`, `packed-int`) that the file states
/// as `NULL` (text) or `05 00` (binary). NCBI's reader leaves the `CRef` empty
/// (`CObjectIStream::ReadPointer`) and makeblastdb dereferences it in `x_ProcessMaskData` (the
/// mask through its reference) or in `CSeq_loc::CheckId` (an element of the mix through its
/// reference, for its range or for its own `CheckId`): a `CCoreException`
/// `eNullPtr`, whose report prints a stack trace with the run's library addresses and the path of
/// the installed binary, which no other program can reproduce.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct NullElement {
    /// The list, in words ("a Seq-loc-mix (mix)").
    pub(crate) list: &'static str,
    /// Where the reader read it (`line N`, or `byte N` of its object).
    pub(crate) position: String,
}

impl NullElement {
    /// LOSAT's explicit refusal of the element, where makeblastdb dereferences it.
    pub(crate) fn refusal(&self) -> MaskFileError {
        MaskFileError::Unsupported(format!(
            "an element of {} stated as NULL ({}), which makeblastdb dereferences: its report of \
             the null pointer holds the run's memory addresses",
            self.list, self.position
        ))
    }
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
