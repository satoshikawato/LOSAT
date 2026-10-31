//! The `Blast4-mask` lists of a BLAST archive: the masks of the results and the `LCaseMask`
//! parameters of the request.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1064-1138
//! ```c++
//! static CRef<CBlast4_mask>
//! s_CreateBlastMask(const CPacked_seqint& packed_int, EBlastProgramType program)
//! {
//!     CRef<CBlast4_mask> retval(new CBlast4_mask);
//!     CRef<CSeq_loc> seqloc(new CSeq_loc);
//!     ITERATE(CPacked_seqint::Tdata, masked_region, packed_int.Get()) {
//!         CRef<CSeq_interval> seqint
//!             (new CSeq_interval(const_cast<CSeq_id&>((*masked_region)->GetId()),
//!                           (*masked_region)->GetFrom(),
//!                           (*masked_region)->GetTo()));
//!         ...
//!         seqloc->SetPacked_int().Set().push_back(seqint);
//!     }
//!     retval->SetLocations().push_back(seqloc);
//!     EBlast4_frame_type frame =
//!         (Blast_QueryIsNucleotide(program) || Blast_QueryIsTranslated(program))
//!         ? eBlast4_frame_type_plus1
//!         : eBlast4_frame_type_notset;
//! ```
//! A mask is one `Seq-loc` `packed-int` of the intervals of a query in one frame (a
//! nucleotide query not translated: `plus1`; a protein query: `notset`; a translated query:
//! the frame of its run of intervals).

use crate::report::seqid::SeqId;
use crate::report::serial::encoder::Encoder;

/// `Blast4-frame-type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Frame {
    NotSet,
    Plus1,
    Plus2,
    Plus3,
    Minus1,
    Minus2,
    Minus3,
}

impl Frame {
    fn name_and_value(self) -> (&'static str, i32) {
        match self {
            Frame::NotSet => ("notset", 0),
            Frame::Plus1 => ("plus1", 1),
            Frame::Plus2 => ("plus2", 2),
            Frame::Plus3 => ("plus3", 3),
            Frame::Minus1 => ("minus1", 4),
            Frame::Minus2 => ("minus2", 5),
            Frame::Minus3 => ("minus3", 6),
        }
    }
}

/// A `Blast4-mask`: the intervals (0-based, inclusive) of the query `id` in one frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaskEntry {
    pub id: String,
    /// With `-parse_deflines`: the record's Seq-align id, written instead of `id`.
    pub best: Option<crate::seqid::ParsedSeqId>,
    pub frame: Frame,
    pub intervals: Vec<(i32, i32)>,
}

/// Writes the mask as a `Blast4-mask` value (the caller wrote `element()` or the variant).
pub fn write_mask<E: Encoder>(enc: &mut E, mask: &MaskEntry) {
    enc.begin_sequence("Blast4-mask");
    enc.member("locations", 0);
    enc.begin_list(false);
    enc.element();
    enc.variant("Seq-loc", "packed-int", 4);
    enc.begin_named("Packed-seqint");
    enc.begin_list(false);
    for &(from, to) in &mask.intervals {
        enc.element();
        enc.begin_sequence("Seq-interval");
        enc.member("from", 0);
        enc.int(from);
        enc.end_member();
        enc.member("to", 1);
        enc.int(to);
        enc.end_member();
        enc.member("id", 3);
        match &mask.best {
            Some(id) => SeqId::Db(id).write(enc),
            None => SeqId::Local(&mask.id).write(enc),
        }
        enc.end_member();
        enc.end_sequence();
    }
    enc.end_list();
    enc.end_named();
    enc.end_variant();
    enc.end_list();
    enc.end_member();
    enc.member("frame", 1);
    let (name, value) = mask.frame.name_and_value();
    enc.begin_named("Blast4-frame-type");
    enc.enumerated(name, value);
    enc.end_named();
    enc.end_member();
    enc.end_sequence();
}
