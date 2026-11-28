//! The `Blast4-mask` lists of a BLAST archive: the masks of the results and the `LCaseMask`
//! parameters of the request.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1063-1138
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1965-1988
    /// `FrameNumber2NetworkFrame` for a program that translates the query: the frame of the
    /// query's translation, 1, 2, 3, -1, -2 or -3 (NCBI aborts on any other; no context of LOSAT's
    /// searches and no mask that LOSAT's archive reader takes has one).
    pub fn of_number(frame: i32) -> Self {
        match frame {
            1 => Frame::Plus1,
            2 => Frame::Plus2,
            3 => Frame::Plus3,
            -1 => Frame::Minus1,
            -2 => Frame::Minus2,
            -3 => Frame::Minus3,
            other => unreachable!("FrameNumber2NetworkFrame aborts on the frame {other}"),
        }
    }

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
    let id = match &mask.best {
        Some(id) => SeqId::Db(id),
        None => SeqId::Local(&mask.id),
    };
    // The intervals are positions of the query, never negative (`TSeqPos` in NCBI).
    write_blast4_mask(
        enc,
        mask.frame,
        mask.intervals
            .iter()
            .map(|&(from, to)| (id, from as u32, to as u32)),
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1063-1138
/// A `Blast4-mask` value of `frame` (the caller wrote `element()` or the variant): its one
/// location, a `Seq-loc` `packed-int` of `intervals` (0-based, inclusive), each on its own Seq-id
/// (`ConvertToRemoteMasks` adds each `CSeqLocInfo` with its Seq-id; `s_CreateBlastMask` copies
/// them without a strand).
pub fn write_blast4_mask<'a, E: Encoder>(
    enc: &mut E,
    frame: Frame,
    intervals: impl IntoIterator<Item = (SeqId<'a>, u32, u32)>,
) {
    enc.begin_sequence("Blast4-mask");
    enc.member("locations", 0);
    enc.begin_list(false);
    enc.element();
    enc.variant("Seq-loc", "packed-int", 4);
    enc.begin_named("Packed-seqint");
    enc.begin_list(false);
    for (id, from, to) in intervals {
        enc.element();
        write_interval(enc, &id, from, to);
    }
    enc.end_list();
    enc.end_named();
    enc.end_variant();
    enc.end_list();
    enc.end_member();
    enc.member("frame", 1);
    let (name, value) = frame.name_and_value();
    enc.begin_named("Blast4-frame-type");
    enc.enumerated(name, value);
    enc.end_named();
    enc.end_member();
    enc.end_sequence();
}

/// A `Seq-loc` `packed-int` value: one `Seq-interval` (`from`, `to`, `id`; no strand) for each
/// interval of `intervals` (0-based, inclusive), each with the Seq-id `id` (the `CSeq_loc` of `id`
/// and the ranges, as the Blast4 masks and the mask writers of dustmasker build it).
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_seqloc.cpp:65-79
pub fn write_packed_int<E: Encoder>(
    enc: &mut E,
    id: &SeqId<'_>,
    intervals: impl IntoIterator<Item = (u32, u32)>,
) {
    enc.variant("Seq-loc", "packed-int", 4);
    enc.begin_named("Packed-seqint");
    enc.begin_list(false);
    for (from, to) in intervals {
        enc.element();
        write_interval(enc, id, from, to);
    }
    enc.end_list();
    enc.end_named();
    enc.end_variant();
}

/// The `CSeq_loc` of `id` and the ranges: a Seq-loc `int` for one interval, else `packed-int`
/// (`write_packed_int`); the mask writers of dustmasker build their Seq-locs so.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:173-182
pub fn write_ranges_seq_loc<E: Encoder>(enc: &mut E, id: &SeqId<'_>, intervals: &[(u32, u32)]) {
    if let [(from, to)] = *intervals {
        enc.variant("Seq-loc", "int", 3);
        write_interval(enc, id, from, to);
        enc.end_variant();
    } else {
        write_packed_int(enc, id, intervals.iter().copied());
    }
}

/// One `Seq-interval` (`from`, `to` and `id`; no `strand`). `from` and `to` are `TSeqPos`,
/// unsigned: `CObjectOStream::WriteUint4` writes a position above 2^31 - 1 as the unsigned
/// number (convert2blastmask 2.17.0, `4294967295 - 4294967294` in the three formats; BER
/// `02 05 00 ff ff ff ff`), which `big_int` writes too, and a smaller one as `int` does.
fn write_interval<E: Encoder>(enc: &mut E, id: &SeqId<'_>, from: u32, to: u32) {
    enc.begin_sequence("Seq-interval");
    enc.member("from", 0);
    enc.big_int(i64::from(from));
    enc.end_member();
    enc.member("to", 1);
    enc.big_int(i64::from(to));
    enc.end_member();
    enc.member("id", 3);
    id.write(enc);
    enc.end_member();
    enc.end_sequence();
}
