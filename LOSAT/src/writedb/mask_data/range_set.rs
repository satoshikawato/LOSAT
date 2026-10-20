//! The masked ranges of the sequences, by algorithm and sequence
//! (`CMaskedRangeSet`, `c++/src/app/blastdb/masked_range_set.cpp`), with the
//! combination of several Seq-locs of one sequence (`CSeq_loc::Add` with
//! `fMerge_All | fSort`).

use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};

use crate::writedb::mask_info::{MaskInterval, MaskLoc, MaskSeqId};
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
pub struct MaskedRanges {
    pub algorithm_id: i32,
    pub offsets: Vec<(u32, u32)>,
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
pub fn ranges_are_empty(ranges: &[MaskedRanges]) -> bool {
    ranges.iter().all(|r| r.offsets.is_empty())
}

/// `CMaskedRangeSet`: the Seq-loc of each (algorithm, sequence) and the
/// ranges of the algorithms that have any, in the order of their first
/// `insert`.
#[derive(Debug, Default)]
pub struct MaskedRangeSet {
    /// `m_Ranges`.
    ranges: Vec<MaskedRanges>,
    /// `m_Values`, by algorithm id: the Seq-loc of each sequence as the
    /// intervals `CSeq_loc_CI` visits in it.
    values: BTreeMap<i32, HashMap<MaskSeqId, Vec<MaskInterval>>>,
}

fn unsupported(message: &str) -> WriteDbError {
    WriteDbError::Other(format!("{message}; LOSAT does not support this mask"))
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
/// `CSeq_loc::Add(other, fMerge_All | fSort)` on the intervals of two
/// Seq-locs of one Seq-id: the intervals of the result in the order
/// `CSeq_loc_CI` visits them.
fn add_locs(
    first: &[MaskInterval],
    second: &[MaskInterval],
) -> Result<Vec<MaskInterval>, WriteDbError> {
    let mut plus: Vec<(u32, u32)> = Vec::new();
    let mut minus: Vec<(u32, u32)> = Vec::new();
    for iv in first.iter().chain(second) {
        if iv.to < iv.from {
            return Err(unsupported(&format!(
                "cannot combine the empty range {}-{} with another Seq-loc of the sequence",
                iv.from, iv.to
            )));
        }
        if iv.to == u32::MAX || (iv.from == 0 && iv.to == u32::MAX - 1) {
            return Err(unsupported(&format!(
                "cannot combine the range {}-{} with another Seq-loc of the sequence",
                iv.from, iv.to
            )));
        }
        if is_reverse(iv.strand) {
            minus.push((iv.from, iv.to));
        } else {
            plus.push((iv.from, iv.to));
        }
    }
    let mut out = Vec::with_capacity(plus.len() + minus.len());
    // `CRange_Less`: by start, then by end.
    plus.sort_unstable();
    push_merged(&mut out, &plus, false);
    // `CRange_ReverseLess`: by decreasing end, then by start.
    minus.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    push_merged(&mut out, &minus, true);
    Ok(out)
}

fn push_merged(out: &mut Vec<MaskInterval>, ranges: &[(u32, u32)], reverse: bool) {
    let mut last: Option<(u32, u32)> = None;
    let strand = reverse.then_some(NA_STRAND_MINUS);
    for &rg in ranges {
        match &mut last {
            Some(last_rg) => {
                if merge_ranges(last_rg, rg, reverse) {
                    continue;
                }
                out.push(MaskInterval {
                    from: last_rg.0,
                    to: last_rg.1,
                    strand,
                });
                *last_rg = rg;
            }
            None => last = Some(rg),
        }
    }
    if let Some(last_rg) = last {
        out.push(MaskInterval {
            from: last_rg.0,
            to: last_rg.1,
            strand,
        });
    }
}

impl MaskedRangeSet {
    pub fn new() -> Self {
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
    /// `Insert`: adds the Seq-loc `loc` of the sequence `id` for the algorithm
    /// `algo_id`. Fails (LOSAT's own error) for a negative id, and for a
    /// combination that LOSAT does not port (`add_locs`).
    pub fn insert(
        &mut self,
        algo_id: i32,
        id: &MaskSeqId,
        loc: &MaskLoc,
    ) -> Result<(), WriteDbError> {
        if algo_id < 0 {
            return Err(unsupported(&format!(
                "the algorithm id {algo_id} is negative"
            )));
        }
        let by_id = self.values.entry(algo_id).or_default();
        if by_id.is_empty() {
            self.ranges.push(MaskedRanges {
                algorithm_id: algo_id,
                offsets: Vec::new(),
            });
        }
        match by_id.get_mut(id) {
            None => {
                by_id.insert(id.clone(), loc.intervals.clone());
            }
            Some(stored) => {
                *stored = add_locs(stored, &loc.intervals)?;
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
    /// `GetRanges`: the ranges of the sequence that has the Seq-ids `ids`,
    /// one `MaskedRanges` for each algorithm that has any mask (in the order
    /// of its first `insert`), with no offsets where the sequence has none.
    pub fn get_ranges(&mut self, ids: &[MaskSeqId]) -> Result<&[MaskedRanges], WriteDbError> {
        let Self { ranges, values } = self;
        for algo in ranges.iter_mut() {
            algo.offsets.clear();
            let Some(by_id) = values.get(&algo.algorithm_id) else {
                continue;
            };
            let mut oid_ranges: Option<Cow<'_, [MaskInterval]>> = None;
            for id in ids {
                if let Some(stored) = by_id.get(id) {
                    oid_ranges = Some(match oid_ranges {
                        None => Cow::Borrowed(stored.as_slice()),
                        Some(l1) => Cow::Owned(add_locs(&l1, stored)?),
                    });
                }
            }
            if let Some(loc) = oid_ranges {
                algo.offsets
                    .extend(loc.iter().map(|iv| (iv.from, iv.to.wrapping_add(1))));
            }
        }
        Ok(ranges)
    }
}
