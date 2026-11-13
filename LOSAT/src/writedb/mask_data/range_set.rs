//! The masked ranges of the sequences, by algorithm and sequence
//! (`CMaskedRangeSet`, `c++/src/app/blastdb/masked_range_set.cpp`), with the
//! combination of several Seq-locs of one sequence (`CSeq_loc::Add` with
//! `fMerge_All | fSort`).
//!
//! A mask finds a sequence the way NCBI's `map<CSeq_id_Handle, CConstRef<CSeq_loc>>` does: through
//! the handle of its Seq-id, and a sequence looks its Seq-ids up in the order of its `GetId()`
//! (`seqid::mask_handle` has the identity of handles: a mask for `embl|AY100001.1` masks
//! `gb|AY100001.1|`, a local string id matches without case, `local str "5"` does not match
//! `lcl|5`, ...). The ranges that several Seq-ids of one sequence find are combined with
//! `CSeq_loc::Add`, which groups them by handle in the handle's order (`m_Packed - 1`, unsigned),
//! sorts and merges each group, and writes the plus strand's groups before the minus strand's.
//!
//! LOSAT refuses (an error that names the Seq-ids, through the mask-data error path) the two
//! cases whose result depends on NCBI's memory addresses: a sequence with a `giim` Seq-id whose
//! number a `giim` Seq-id of the mask data has (NCBI's Giim tree finds the first stored handle of
//! the number that is *not* `Equals`), and ranges of one sequence and algorithm that come from two
//! different handles with the same `m_Packed` (NCBI orders those by the address of the handle's
//! info).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::seqid::mask_handle::{
    mask_handle_key, mask_handle_order, mask_handle_packed, Undecidable,
};
use crate::seqid::ParsedSeqId;
use crate::writedb::mask_info::{MaskInterval, MaskLoc};
use crate::writedb::WriteDbError;

// NCBI reference (598d8ae6): c++/include/objects/blastdb/defline_extra.hpp:78-93
// ```c++
// struct SBlastDbMaskData {
//
//     /// Construct and clear this object.
//     SBlastDbMaskData() : algorithm_id((int)eBlast_filter_program_not_set) {}
//
//     /// Identifies the algorithm used.
//     int algorithm_id;
//
//     /// Start and end offsets of the filtered area.
//     vector< pair<TSeqPos, TSeqPos> > offsets;
//
//     /// Determines whether this object contains no offsets
//     bool empty() const {
//         return offsets.empty();
//     }
// };
// ```
/// `SBlastDbMaskData`: the masked ranges of one sequence for one algorithm,
/// each `(start, end)` with `end` exclusive.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct MaskedRanges {
    pub(crate) algorithm_id: i32,
    pub(crate) offsets: Vec<(u32, u32)>,
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb.hpp:64-79
// ```c++
// class NCBI_XOBJWRITE_EXPORT CMaskedRangesVector: public vector<SBlastDbMaskData>
// {
// public:
//     /// Our parent class
//     typedef vector<SBlastDbMaskData> TParent;
//
//     /// Redefine empty to mean no elements or none of its elements being empty
//     bool empty() const {
//         ITERATE(TParent, itr, *this) {
//             if ( !itr->empty() ) {
//                 return false;
//             }
//         }
//         return true;
//     }
// };
// ```
/// `CMaskedRangesVector::empty`: no algorithm has a range.
pub(crate) fn ranges_are_empty(ranges: &[MaskedRanges]) -> bool {
    ranges.iter().all(|r| r.offsets.is_empty())
}

/// `CMaskedRangeSet`: the Seq-loc of each (algorithm, handle) and the
/// ranges of the algorithms that have any, in the order of their first
/// `insert`.
#[derive(Debug, Default)]
pub(crate) struct MaskedRangeSet {
    /// `m_Ranges`.
    ranges: Vec<MaskedRanges>,
    /// `m_Values`, by algorithm id.
    values: BTreeMap<i32, AlgorithmLocs>,
    /// The numbers of the `giim` Seq-ids of every mask inserted (NCBI's Giim tree is shared by
    /// the algorithms): the handle of a sequence's `giim` id depends on them.
    giim_numbers: BTreeSet<i32>,
}

/// The Seq-locs of one algorithm.
#[derive(Debug, Default)]
struct AlgorithmLocs {
    /// `m_Values[algo_id]` has an entry (any `Insert` makes one, `giim` and `neither accession nor
    /// name` Seq-ids included: `x_Set` creates the map element before it combines).
    used: bool,
    /// The handles that can be found again, by `mask_handle_key`.
    by_handle: HashMap<Vec<u8>, Stored>,
}

/// The Seq-loc of an (algorithm, handle): the intervals NCBI holds, or the
/// error of a combination that LOSAT does not port. NCBI combines Seq-locs
/// for every handle, including ones no sequence has, and finishes the run;
/// the error is therefore raised only when the sequence is looked up.
#[derive(Debug)]
struct Stored {
    /// A number of the handle, unique among the handles of the algorithm: the groups of
    /// ranges of one handle.
    handle: usize,
    /// `m_Packed` of the handle (`mask_handle_packed`).
    packed: i64,
    loc: StoredLoc,
}

#[derive(Debug)]
enum StoredLoc {
    Intervals(Vec<MaskInterval>),
    Unsupported(String),
}

/// One interval of a combination, with the handle of its Seq-id.
#[derive(Clone, Copy)]
struct Piece {
    handle: usize,
    packed: i64,
    interval: MaskInterval,
}

impl Stored {
    /// The error of a Seq-loc that LOSAT did not combine.
    fn refusal(&self) -> Result<(), WriteDbError> {
        match &self.loc {
            StoredLoc::Intervals(_) => Ok(()),
            StoredLoc::Unsupported(message) => Err(WriteDbError::Other(message.clone())),
        }
    }

    fn pieces(&self) -> Result<impl Iterator<Item = Piece> + '_, WriteDbError> {
        match &self.loc {
            StoredLoc::Intervals(intervals) => Ok(intervals.iter().map(|&interval| Piece {
                handle: self.handle,
                packed: self.packed,
                interval,
            })),
            StoredLoc::Unsupported(message) => Err(WriteDbError::Other(message.clone())),
        }
    }
}

/// Why `add_locs` does not give a Seq-loc.
enum AddError {
    /// LOSAT does not port the combination.
    Unsupported(String),
    /// Two different handles with the same `m_Packed`: NCBI orders their ranges by the address of
    /// the handles' infos.
    Tie(usize, usize),
}

fn unsupported(message: &str) -> String {
    format!("{message}; LOSAT does not support this mask")
}

// NCBI reference (598d8ae6): c++/include/util/range.hpp:142-145
// ```c++
//     position_type GetTo(void) const
//         {
//             return GetToOpen()-1;
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4283-4298
// ```c++
// class CRange_Less
// {
// public:
//     CRange_Less(void) {}
//
//     bool operator() (const TRangeWithFuzz& rg1, const TRangeWithFuzz& rg2) const
//     {
//         if ( rg1.IsWhole() ) {
//             return !rg2.IsWhole();
//         }
//         if ( rg1.Empty() ) {
//             return !rg2.Empty()  &&  !rg2.IsWhole();
//         }
//         return !rg2.IsWhole()  &&  !rg2.Empty()  &&  rg1 < rg2;
//     }
// };
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4301-4322
// ```c++
// class CRange_ReverseLess
// {
// public:
//     CRange_ReverseLess(void) {}
//
//     bool operator() (const TRangeWithFuzz& rg1, const TRangeWithFuzz& rg2) const
//     {
//         if ( rg1.IsWhole() ) {
//             return !rg2.IsWhole();
//         }
//         if ( rg1.Empty() ) {
//             return !rg2.Empty()  &&  !rg2.IsWhole();
//         }
//         if ( rg2.IsWhole()  ||  rg2.Empty() ) {
//             return false;
//         }
//         if (rg1.GetTo() != rg2.GetTo()) {
//             return rg1.GetTo() > rg2.GetTo();
//         }
//         return rg1.GetFrom() < rg2.GetFrom();
//     }
// };
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4333-4389
// ```c++
// static
// bool x_MergeRanges(TRangeWithFuzz& rg1, ENa_strand str1,
//                    const TRangeWithFuzz& rg2, ENa_strand str2,
//                    CSeq_loc::TOpFlags flags)
// {
//     if ( !x_MatchStrand(str1, str2, flags) ) {
//         return false;
//     }
//     // Check contained
//     if ( (flags & CSeq_loc::fMerge_Contained) != 0 ) {
//         if (rg1.GetFrom() <= rg2.GetFrom()  &&  rg1.GetTo() >= rg2.GetTo()) {
//             // rg2 already contained in rg1
//             ...
//             return true;
//         }
//         if (rg1.GetFrom() >= rg2.GetFrom()  &&  rg1.GetTo() <= rg2.GetTo()) {
//             // rg1 contained in rg2
//             ...
//             rg1 = rg2;
//             ...
//             return true;
//         }
//     }
//     // Check overlapping
//     if ( (flags & CSeq_loc::fMerge_OverlappingOnly) != 0  &&
//         rg1.IntersectingWith(rg2) ) {
//         rg1 += rg2;
//         return true;
//     }
//     // Check abutting
//     if ((flags & CSeq_loc::fMerge_AbuttingOnly) != 0) {
//         if ( !IsReverse(str1) ) {
//             if ( rg1.GetToOpen() == rg2.GetFrom() ) {
//                 rg1.CopyTo(rg2);
//                 return true;
//             }
//         }
//         else {
//             if (rg1.GetFrom() == rg2.GetToOpen()) {
//                 rg1.CopyFrom(rg2);
//                 return true;
//             }
//         }
//     }
//     return false;
// }
// ```
// The ranges here are `(from, to)` with `to` inclusive. `fMerge_All` sets
// all three merge bits, and the strand of both ranges of a comparison is
// the strand of their group.
fn merge_ranges(rg1: &mut (u32, u32), rg2: (u32, u32), reverse: bool) -> bool {
    // Contained.
    if rg1.0 <= rg2.0 && rg1.1 >= rg2.1 {
        return true;
    }
    if rg1.0 >= rg2.0 && rg1.1 <= rg2.1 {
        *rg1 = rg2;
        return true;
    }
    // Overlapping (`IntersectingWith`: the intersection is not empty).
    if rg1.0 <= rg2.1 && rg2.0 <= rg1.1 {
        rg1.0 = rg1.0.min(rg2.0);
        rg1.1 = rg1.1.max(rg2.1);
        return true;
    }
    // Abutting. `to` is below `u32::MAX` (checked by the caller).
    if !reverse {
        if rg1.1 + 1 == rg2.0 {
            rg1.1 = rg2.1;
            return true;
        }
    } else if rg1.0 == rg2.1 + 1 {
        rg1.0 = rg2.0;
        return true;
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:115-121
// ```asn
// Na-strand ::= ENUMERATED {          -- strand of nucleic acid
//     unknown (0) ,
//     plus (1) ,
//     minus (2) ,
//     both (3) ,                -- in forward orientation
//     both-rev (4) ,            -- in reverse orientation
//     other (255) }
// ```
// NCBI reference (598d8ae6): c++/include/objects/seqloc/Na_strand.hpp:74-79
// ```c++
// inline
// bool IsReverse(ENa_strand s)
// {
//     // treat unknown as forward
//     return (s == eNa_strand_minus  ||  s == eNa_strand_both_rev);
// }
// ```
const NA_STRAND_MINUS: u8 = 2;
const NA_STRAND_BOTH_REV: u8 = 4;

fn is_reverse(strand: Option<u8>) -> bool {
    matches!(strand, Some(NA_STRAND_MINUS | NA_STRAND_BOTH_REV))
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:5060-5084
// ```c++
// CRef<CSeq_loc> CSeq_loc::Add(const CSeq_loc& other,
//                              TOpFlags        flags,
//                              ISynonymMapper* syn_mapper) const
// {
//     ...
//     CRef<CSeq_loc> ret(new CSeq_loc);
//     CSeq_loc tmp;
//     tmp.SetMix().AddSeqLoc(const_cast<CSeq_loc&>(*this));
//     tmp.SetMix().AddSeqLoc(const_cast<CSeq_loc&>(other));
//     if ( (flags & CSeq_loc::fMerge_SingleRange) != 0 ) {
//         x_SingleRange(*ret, tmp, *syn_mapper);
//     }
//     else if ( (flags & CSeq_loc::fSort) == 0 ) {
//         x_MergeNoSort(*ret, tmp, flags, *syn_mapper);
//     }
//     else {
//         x_MergeAndSort(*ret, tmp, flags, *syn_mapper);
//     }
//     return ret;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4640-4693
// ```c++
// static
// void x_MergeAndSort(CSeq_loc& dst,
//                     const CSeq_loc& src,
//                     CSeq_loc::TOpFlags flags,
//                     ISynonymMapper& syn_mapper)
// {
//     bool use_strand = (flags & CSeq_loc::fStrand_Ignore) == 0;
//     // Id -> range map for both strands
//     ...
//     // Split location by by id/strand/range
//     for (CSeq_loc_CI it(src, CSeq_loc_CI::eEmpty_Allow); it; ++it) {
//         CSeq_id_Handle_Wrapper idh(syn_mapper.GetBestSynonym(it.GetSeq_id()), it.GetSeq_id());
//         ENa_strand strand = it.GetStrand();
//         if ( IsReverse(strand) ) {
//             id_map_minus[idh].push_back(TRangeWithFuzz(it));
//             ...
//         }
//         else {
//             id_map_plus[idh].push_back(TRangeWithFuzz(it));
//             ...
//         }
//     }
//
//     x_RangesToSeq_loc(dst, id_map_plus, default_plus, flags);
//     if ( use_strand ) {
//         x_RangesToSeq_loc(dst, id_map_minus, default_minus, flags);
//     }
//     if (dst.Which() == CSeq_loc::e_not_set) {
//         dst.SetNull();
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4546-4594
// ```c++
// static
// void x_RangesToSeq_loc(CSeq_loc& dst,
//                        TIdToRangeMap& id_map,
//                        ENa_strand default_strand,
//                        CSeq_loc::TOpFlags flags)
// {
//     // Iterate ids for each strand
//     NON_CONST_ITERATE(TIdToRangeMap, id_it, id_map) {
//         ...
//         TRanges& ranges = id_it->second;
//         if ( (flags & CSeq_loc::fSort) != 0 ) {
//             if ( !IsReverse(default_strand) ) {
//                 sort(ranges.begin(), ranges.end(), CRange_Less());
//             }
//             else {
//                 sort(ranges.begin(), ranges.end(), CRange_ReverseLess());
//             }
//         }
//         // Merge ranges according to the flags, add to destination
//         TRangeWithFuzz last_rg(TRangeWithFuzz::GetEmpty());
//         bool have_range = false;
//         ITERATE(TRanges, rg, ranges) {
//             if (x_MergeRanges(last_rg, default_strand,
//                             *rg, default_strand,
//                             flags)) {
//                 have_range = true;
//                 continue;
//             }
//             // No merging - push current range, reset last values
//             if (have_range) {
//                 x_PushRange(dst, id_it->first, last_rg, default_strand);
//             }
//             last_rg = *rg;
//             have_range = true;
//         }
//         if (have_range) {
//             x_PushRange(dst, id_it->first, last_rg, default_strand);
//         }
//     }
// }
// ```
// The result is a Seq-loc whose intervals `CSeq_loc_CI` visits in the order
// they were pushed: the ranges of the plus-like strands, ascending, then the
// ranges of the reverse strands, by decreasing end. All intervals of a Seq-loc
// have the same Seq-id (the key of `m_Values`), so the map by Seq-id has one
// element. A range whose end is before its start (an empty `CRange`), one
// that ends at 4294967295 and the one that covers the whole range are
// handled by NCBI through `CRange`'s empty and whole cases; LOSAT does not
// port those and stops.
/// `CSeq_loc::Add(other, fMerge_All | fSort)` on two Seq-locs: the intervals of the result in the
/// order `CSeq_loc_CI` visits them, with their handles. They are grouped by handle (`x_MergeAndSort`:
/// `map<CSeq_id_Handle_Wrapper, TRanges>`, ordered by `CSeq_id_Handle::operator<`), the ranges of
/// each group sorted and merged, the groups of the plus-like strands written first, then those of
/// the reverse strands. The groups of different handles are never merged with each other.
fn add_locs(
    first: impl Iterator<Item = Piece>,
    second: impl Iterator<Item = Piece>,
) -> Result<Vec<Piece>, AddError> {
    // By the handle's place in `operator<`: (handle, packed, plus ranges, minus ranges).
    type Group = (usize, i64, Vec<(u32, u32)>, Vec<(u32, u32)>);
    let mut groups: BTreeMap<u64, Group> = BTreeMap::new();
    for piece in first.chain(second) {
        let iv = piece.interval;
        if iv.to < iv.from {
            return Err(AddError::Unsupported(unsupported(&format!(
                "cannot combine the empty range {}-{} with another Seq-loc of the sequence",
                iv.from, iv.to
            ))));
        }
        if iv.to == u32::MAX || (iv.from == 0 && iv.to == u32::MAX - 1) {
            return Err(AddError::Unsupported(unsupported(&format!(
                "cannot combine the range {}-{} with another Seq-loc of the sequence",
                iv.from, iv.to
            ))));
        }
        let group = groups
            .entry(mask_handle_order(piece.packed))
            .or_insert_with(|| (piece.handle, piece.packed, Vec::new(), Vec::new()));
        if group.0 != piece.handle {
            return Err(AddError::Tie(group.0, piece.handle));
        }
        if is_reverse(iv.strand) {
            group.3.push((iv.from, iv.to));
        } else {
            group.2.push((iv.from, iv.to));
        }
    }
    let mut out = Vec::new();
    // `CRange_Less`: by start, then by end.
    for (handle, packed, plus, _) in groups.values_mut() {
        plus.sort_unstable();
        push_merged(&mut out, *handle, *packed, plus, false);
    }
    // `CRange_ReverseLess`: by decreasing end, then by start.
    for (handle, packed, _, minus) in groups.values_mut() {
        minus.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        push_merged(&mut out, *handle, *packed, minus, true);
    }
    Ok(out)
}

fn push_merged(
    out: &mut Vec<Piece>,
    handle: usize,
    packed: i64,
    ranges: &[(u32, u32)],
    reverse: bool,
) {
    let mut last: Option<(u32, u32)> = None;
    let strand = reverse.then_some(NA_STRAND_MINUS);
    let mut push = |rg: (u32, u32)| {
        out.push(Piece {
            handle,
            packed,
            interval: MaskInterval {
                from: rg.0,
                to: rg.1,
                strand,
            },
        });
    };
    for &rg in ranges {
        match &mut last {
            Some(last_rg) => {
                if merge_ranges(last_rg, rg, reverse) {
                    continue;
                }
                push(*last_rg);
                *last_rg = rg;
            }
            None => last = Some(rg),
        }
    }
    if let Some(last_rg) = last {
        push(last_rg);
    }
}

impl MaskedRangeSet {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/masked_range_set.hpp:48-53
    // ```c++
    //     void Insert(int              algo_id,
    //                 const CSeq_id  & id,
    //                 const CSeq_loc & v)
    //     {
    //         x_CombineLocs(x_Set(algo_id, CSeq_id_Handle::GetHandle(id)), v);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blastdb/masked_range_set.cpp:86-106
    // ```c++
    // void CMaskedRangeSet::x_CombineLocs(CConstRef<CSeq_loc> & L1,
    //                                     const CSeq_loc      & L2)
    // {
    //     if (L1.Empty()) {
    //         L1.Reset(& L2);
    //     } else {
    //         L1 = L1->Add(L2, CSeq_loc::fMerge_All | CSeq_loc::fSort, NULL);
    //     }
    // }
    //
    // CConstRef<CSeq_loc> & CMaskedRangeSet::x_Set(int algo_id, CSeq_id_Handle idh)
    // {
    //     if ((int)m_Values.size() <= algo_id) {
    //         m_Values.resize(algo_id+1);
    //     }
    //     if (m_Values[algo_id].empty()) {
    //         m_Ranges.resize(m_Ranges.size()+1);
    //         m_Ranges[m_Ranges.size()-1].algorithm_id = algo_id;
    //     }
    //     return m_Values[algo_id][idh];
    // }
    // ```
    // The first Seq-loc of an (algorithm, sequence) is kept as it was given:
    // unsorted, overlapping, with its intervals in the order `CSeq_loc_CI`
    // visits them (`loc.intervals`). Each later one is added to the stored
    // one with `CSeq_loc::Add`.
    /// `Insert`: adds the Seq-loc `loc` of the Seq-id `id` for the algorithm
    /// `algo_id`, under the Seq-id's handle. Fails (LOSAT's own error) for a negative
    /// id, and for a combination that LOSAT does not port (`add_locs`).
    pub(crate) fn insert(
        &mut self,
        algo_id: i32,
        id: &ParsedSeqId,
        loc: &MaskLoc,
    ) -> Result<(), WriteDbError> {
        if algo_id < 0 {
            return Err(WriteDbError::Other(unsupported(&format!(
                "the algorithm id {algo_id} is negative"
            ))));
        }
        let algorithm = self.values.entry(algo_id).or_default();
        if !algorithm.used {
            algorithm.used = true;
            self.ranges.push(MaskedRanges {
                algorithm_id: algo_id,
                offsets: Vec::new(),
            });
        }
        let key = match mask_handle_key(id) {
            Ok(Some(key)) => key,
            // A Seq-id with neither accession nor name: a new handle that nothing finds again.
            Ok(None) => return Ok(()),
            // The handle of a `giim` id depends on the order of the ids made before it; a
            // sequence that has such an id of the number is refused when it is looked up.
            Err(Undecidable::Giim { id }) => {
                self.giim_numbers.insert(id);
                return Ok(());
            }
        };
        let next_handle = algorithm.by_handle.len();
        match algorithm.by_handle.get_mut(&key) {
            None => {
                algorithm.by_handle.insert(
                    key,
                    Stored {
                        handle: next_handle,
                        packed: mask_handle_packed(id),
                        loc: StoredLoc::Intervals(loc.intervals.clone()),
                    },
                );
            }
            // Every later Seq-loc of a refused one is combined with it: nothing to do.
            Some(Stored {
                loc: StoredLoc::Unsupported(_),
                ..
            }) => {}
            Some(entry) => {
                let new = Stored {
                    handle: entry.handle,
                    packed: entry.packed,
                    loc: StoredLoc::Intervals(loc.intervals.clone()),
                };
                let combined = add_locs(entry.pieces()?, new.pieces()?);
                entry.loc = match combined {
                    Ok(pieces) => {
                        StoredLoc::Intervals(pieces.into_iter().map(|p| p.interval).collect())
                    }
                    Err(AddError::Unsupported(message)) => StoredLoc::Unsupported(message),
                    // One handle: no two handles to order.
                    Err(AddError::Tie(..)) => {
                        StoredLoc::Unsupported(unsupported("two handles in one Seq-id"))
                    }
                };
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/masked_range_set.cpp:38-69
    // ```c++
    // CMaskedRangesVector &
    // CMaskedRangeSet::GetRanges(const list< CRef<CSeq_id> > & idlist)
    // {
    //     // For each algorithm for which data is provided
    //
    //     NON_CONST_ITERATE(CMaskedRangesVector, algo_iter, m_Ranges) {
    //         algo_iter->offsets.resize(0);
    //         int algo_id = algo_iter->algorithm_id;
    //
    //         // Combine Seq-loc ranges for all provided Seq-ids.
    //
    //         CConstRef<CSeq_loc> oid_ranges;
    //
    //         ITERATE(list< CRef<CSeq_id> >, id_iter, idlist) {
    //             CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(**id_iter);
    //             x_FindAndCombine(oid_ranges, algo_id, idh);
    //         }
    //
    //         if (oid_ranges.Empty()) {
    //             continue;
    //         }
    //
    //         for(CSeq_loc_CI ranges(*oid_ranges); ranges; ++ranges) {
    //             CSeq_loc::TRange rng = ranges.GetRange();
    //             pair<TSeqPos, TSeqPos> pr(rng.GetFrom(), rng.GetToOpen());
    //
    //             algo_iter->offsets.push_back(pr);
    //         }
    //     }
    //
    //     return m_Ranges;
    // }
    // ```
    // `GetToOpen` of a range ending at 4294967295 is 0 (a `TSeqPos` sum),
    // and a range whose end is before its start is passed on as it is
    // (checked by the writer).
    /// `GetRanges`: the ranges of the sequence that has the Seq-ids `ids` (`bs->GetId()`), looked
    /// up in that order: one `MaskedRanges` for each algorithm that has any mask (in the order
    /// of its first `insert`), with no offsets where the sequence has none.
    pub(crate) fn get_ranges(
        &mut self,
        ids: &[ParsedSeqId],
    ) -> Result<&[MaskedRanges], WriteDbError> {
        let Self {
            ranges,
            values,
            giim_numbers,
        } = self;
        // The handle of each Seq-id; a `giim` id of a number that the masks have has none.
        let mut keys = Vec::with_capacity(ids.len());
        for id in ids {
            keys.push(match mask_handle_key(id) {
                Ok(key) => key,
                Err(Undecidable::Giim { id: number }) if giim_numbers.contains(&number) => {
                    return Err(WriteDbError::Other(unsupported(&format!(
                        "the sequence has the Seq-id {} and the mask data have a giim Seq-id with the number \
                         {number}: NCBI's handle for a giim Seq-id depends on the order in which the \
                         Seq-ids were made",
                        id_label(id)
                    ))));
                }
                Err(Undecidable::Giim { .. }) => None,
            });
        }
        for algo in ranges.iter_mut() {
            algo.offsets.clear();
            let Some(locs) = values.get(&algo.algorithm_id) else {
                continue;
            };
            let mut found: Option<Found<'_>> = None;
            let mut found_by: Vec<(usize, &ParsedSeqId)> = Vec::new();
            for (id, key) in ids.iter().zip(&keys) {
                let Some(stored) = key.as_ref().and_then(|key| locs.by_handle.get(key)) else {
                    continue;
                };
                found_by.push((stored.handle, id));
                found = Some(match found {
                    None => {
                        // The first Seq-loc is taken as it is.
                        stored.refusal()?;
                        Found::One(stored)
                    }
                    Some(previous) => {
                        let combined = match &previous {
                            Found::One(one) => add_locs(one.pieces()?, stored.pieces()?),
                            Found::Many(pieces) => {
                                add_locs(pieces.iter().copied(), stored.pieces()?)
                            }
                        };
                        Found::Many(combined.map_err(|error| match error {
                            AddError::Unsupported(message) => WriteDbError::Other(message),
                            AddError::Tie(a, b) => tie_error(&found_by, a, b),
                        })?)
                    }
                });
            }
            match found {
                None => {}
                Some(Found::One(stored)) => {
                    if let StoredLoc::Intervals(intervals) = &stored.loc {
                        algo.offsets
                            .extend(intervals.iter().map(|iv| (iv.from, iv.to.wrapping_add(1))));
                    }
                }
                Some(Found::Many(pieces)) => algo.offsets.extend(
                    pieces
                        .iter()
                        .map(|p| (p.interval.from, p.interval.to.wrapping_add(1))),
                ),
            }
        }
        Ok(ranges)
    }
}

/// The Seq-loc that `GetRanges` has combined so far (`oid_ranges`).
enum Found<'a> {
    /// One stored Seq-loc, as it was given.
    One(&'a Stored),
    /// The result of `CSeq_loc::Add`.
    Many(Vec<Piece>),
}

/// The Seq-id as `AsFastaString()` prints it.
fn id_label(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    crate::seqid::label::fasta(id, &mut out);
    String::from_utf8_lossy(&out).into_owned()
}

/// The refusal of a tie between the handles numbered `a` and `b`, named by the Seq-ids of the
/// sequence that found them.
fn tie_error(found_by: &[(usize, &ParsedSeqId)], a: usize, b: usize) -> WriteDbError {
    let name = |handle: usize| {
        found_by
            .iter()
            .find(|(h, _)| *h == handle)
            .map_or_else(String::new, |(_, id)| id_label(id))
    };
    WriteDbError::Other(unsupported(&format!(
        "the masks of the Seq-ids {} and {} of one sequence combine in an order that NCBI takes \
         from the memory addresses of the two handles (their packed numbers are equal)",
        name(a),
        name(b)
    )))
}
