//! The ranges of NCBI's object manager: `COpenRange<TSeqPos>` (`Range`), a Seq-id's ranges of a
//! location (`CHandleRange`, `HandleRange`) and a location's ranges by Seq-id
//! (`CHandleRangeMap`, `HandleRangeMap`), which the annotation index asks for each Seq-id's
//! overlapping range.

use super::view::V;
use crate::seqid::mask_handle::{
    mask_handle_key, mask_handle_order, mask_handle_packed, TreeHandle,
};
use crate::seqid::ParsedSeqId;

/// `GetPositionMax()` of `TSeqPos`.
pub(super) const MAX: u32 = u32::MAX;

// NCBI reference (598d8ae6): c++/include/util/range.hpp:148-151, 255-308, 345-357, 398-401
// ```c++
//     bool Empty(void) const
//         {
//             return GetToOpen() <= GetFrom();
//         }
// ...
//     static TThisType GetWhole(void)
//         {
//             return TThisType(GetWholeFrom(), GetWholeToOpen());
//         }
// ...
//     static TThisType GetEmpty(void)
//         {
//             return TThisType(GetEmptyFrom(), GetEmptyToOpen());
//         }
// ...
//     TThisType& CombineWith(const TThisType& r)
//         {
//             if ( !r.Empty() ) {
//                 if ( !Empty() ) {
//                     m_From = min(m_From, r.GetFrom());
//                     m_ToOpen = max(m_ToOpen, r.GetToOpen());
//                 }
//                 else {
//                     *this = r;
//                 }
//             }
//             return *this;
//         }
// ...
//     CRange(position_type from, position_type to)
//         : TParent(from, to+1)
//         {
//         }
// ```
/// `COpenRange<TSeqPos>`: `from` and the open end `to_open`, empty when `to_open <= from`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Range {
    pub(super) from: u32,
    pub(super) to_open: u32,
}

impl Range {
    /// `GetEmpty()`.
    pub(super) const EMPTY: Range = Range {
        from: MAX,
        to_open: MAX,
    };
    /// `GetWhole()`.
    pub(super) const WHOLE: Range = Range {
        from: 0,
        to_open: MAX,
    };

    /// `CRange<TSeqPos>(from, to)`: the open end `to + 1`, wrapping (`to` 4294967295 makes an
    /// empty range).
    pub(super) fn new(from: u32, to: u32) -> Range {
        Range {
            from,
            to_open: to.wrapping_add(1),
        }
    }

    pub(super) fn is_empty(self) -> bool {
        self.to_open <= self.from
    }

    /// `operator+=` (`CombineWith`).
    pub(super) fn combine(&mut self, other: Range) {
        if !other.is_empty() {
            if !self.is_empty() {
                self.from = self.from.min(other.from);
                self.to_open = self.to_open.max(other.to_open);
            } else {
                *self = other;
            }
        }
    }
}

/// `ENa_strand` (`Na-strand`): unknown 0, plus 1, minus 2, both 3, both-rev 4, other 255; a
/// Seq-table column's value of any `int`.
pub(super) type Strand = i64;
/// `eNa_strand_unknown`.
pub(super) const UNKNOWN: Strand = 0;
const MINUS: Strand = 2;
const BOTH: Strand = 3;
const BOTH_REV: Strand = 4;

// NCBI reference (598d8ae6): c++/include/objmgr/impl/handle_range.hpp:185-199
// ```c++
// bool CHandleRange::x_IncludesPlus(const ENa_strand& strand) const
// {
//     // Anything but "minus" includes "plus"
//     return strand != eNa_strand_minus;
// }
//
//
// inline
// bool CHandleRange::x_IncludesMinus(const ENa_strand& strand) const
// {
//     return strand == eNa_strand_unknown
//         || strand == eNa_strand_minus
//         ||  strand == eNa_strand_both
//         ||  strand == eNa_strand_both_rev;
// }
// ```
fn includes_plus(strand: Strand) -> bool {
    strand != MINUS
}

fn includes_minus(strand: Strand) -> bool {
    matches!(strand, UNKNOWN | MINUS | BOTH | BOTH_REV)
}

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Na_strand.hpp:75-79
// ```c++
// bool IsReverse(ENa_strand s)
// {
//     // treat unknown as forward
//     return (s == eNa_strand_minus  ||  s == eNa_strand_both_rev);
// }
// ```
fn is_reverse(strand: Strand) -> bool {
    matches!(strand, MINUS | BOTH_REV)
}

/// `CHandleRange`: one Seq-id's ranges of a location, with their totals by strand.
#[derive(Clone, Debug)]
pub(super) struct HandleRange {
    ranges: Vec<(Range, Strand)>,
    plus: Range,
    minus: Range,
    pub(super) circular: bool,
    single_strand: bool,
    more_before: bool,
    more_after: bool,
}

impl Default for HandleRange {
    fn default() -> Self {
        HandleRange {
            ranges: Vec::new(),
            plus: Range::EMPTY,
            minus: Range::EMPTY,
            circular: false,
            single_strand: true,
            more_before: false,
            more_after: false,
        }
    }
}

impl HandleRange {
    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:126-206
    // ```c++
    // void CHandleRange::AddRange(TRange range, ENa_strand strand,
    //                             bool more_before, bool more_after, bool circular_rna)
    // {
    //     if ( !m_Ranges.empty()  &&  m_IsSingleStrand ) {
    //         if ( strand != m_Ranges.front().second ) {
    //             // Different strands, the location can not be circular
    //             if ( m_IsCircular ) {
    //                 // Different strands, the location can not be circular
    //                 // Reorganize total ranges by strand
    //                 TRange total_range = m_TotalRanges_plus;
    //                 total_range += m_TotalRanges_minus;
    //                 if ( x_IncludesPlus(m_Ranges.front().second) ) {
    //                     m_TotalRanges_plus = total_range;
    //                 }
    //                 else {
    //                     m_TotalRanges_plus = TRange::GetEmpty();
    //                 }
    //                 if ( x_IncludesMinus(m_Ranges.front().second) ) {
    //                     m_TotalRanges_minus = total_range;
    //                 }
    //                 else {
    //                     m_TotalRanges_minus = TRange::GetEmpty();
    //                 }
    //                 m_IsCircular = false;
    //             }
    //             m_IsSingleStrand = false;
    //         }
    //         else {
    //             // Same strand, but location may become circular
    //             if ( !m_IsCircular && !circular_rna ) {
    //                 // Check if location becomes circular
    //                 NON_CONST_REVERSE_ITERATE ( TRanges, it, m_Ranges ) {
    //                     // compare with last non-empty range
    //                     if ( !it->first.Empty() ) {
    //                         if ( x_IncludesPlus(strand) ) {
    //                             m_IsCircular =
    //                                 range.GetFrom() < it->first.GetFrom();
    //                         }
    //                         else {
    //                             m_IsCircular =
    //                                 range.GetFrom() > it->first.GetFrom();
    //                         }
    //                         break;
    //                     }
    //                 }
    //                 if ( m_IsCircular ) {
    //                     ...
    //                     m_TotalRanges_plus += m_TotalRanges_minus;
    //                     m_TotalRanges_minus = TRange::GetEmpty();
    //                 }
    //                 else {
    //                     ...
    //                     if ( more_after ) {
    //                         m_MoreAfter = true;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     else {
    //         if ( more_before ) {
    //             m_MoreBefore = true;
    //         }
    //         if ( more_after ) {
    //             m_MoreAfter = true;
    //         }
    //     }
    //     m_Ranges.push_back(TRanges::value_type(range, strand));
    //     if ( !m_IsCircular ) {
    //         // Regular location
    //         if ( x_IncludesPlus(strand) ) {
    //             m_TotalRanges_plus += range;
    //         }
    //         if ( x_IncludesMinus(strand) ) {
    //             m_TotalRanges_minus += range;
    //         }
    //     }
    //     else {
    //         // Circular location, second part
    //         m_TotalRanges_minus += range;
    //     }
    // }
    // ```
    // The map's `AddRange` (handle_range_map.cpp:277) passes `more_before` and `more_after`
    // false: it sets the flags itself (`HandleRangeMap::add`).
    /// `AddRange(range, strand, false, false, circular_rna)`.
    pub(super) fn add(&mut self, range: Range, strand: Strand, circular_rna: bool) {
        if let Some(&(_, first)) = self.ranges.first() {
            if self.single_strand {
                if strand != first {
                    if self.circular {
                        let mut total = self.plus;
                        total.combine(self.minus);
                        self.plus = if includes_plus(first) {
                            total
                        } else {
                            Range::EMPTY
                        };
                        self.minus = if includes_minus(first) {
                            total
                        } else {
                            Range::EMPTY
                        };
                        self.circular = false;
                    }
                    self.single_strand = false;
                } else if !self.circular && !circular_rna {
                    if let Some(&(last, _)) = self
                        .ranges
                        .iter()
                        .rev()
                        .find(|(range, _)| !range.is_empty())
                    {
                        self.circular = if includes_plus(strand) {
                            range.from < last.from
                        } else {
                            range.from > last.from
                        };
                    }
                    if self.circular {
                        let minus = self.minus;
                        self.plus.combine(minus);
                        self.minus = Range::EMPTY;
                    }
                }
            }
        }
        self.ranges.push((range, strand));
        if !self.circular {
            if includes_plus(strand) {
                self.plus.combine(range);
            }
            if includes_minus(strand) {
                self.minus.combine(range);
            }
        } else {
            self.minus.combine(range);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:338-376
    // ```c++
    // CHandleRange::GetOverlappingRange(TTotalRangeFlags flags) const
    // {
    //     TRange ret = TRange::GetEmpty();
    //     if (m_IsCircular) {
    //         ETotalRangeFlags circular_strand =
    //             IsReverse(m_Ranges.front().second) ?
    //             eStrandMinus : eStrandPlus;
    //         if (flags & circular_strand) {
    //             ret = TRange::GetWhole();
    //         }
    //         return ret;
    //     }
    //     if (flags & eStrandPlus) {   // == eCircularStart
    //         ret += m_TotalRanges_plus;
    //     }
    //     if (flags & eStrandMinus) {  // == eCircularEnd
    //         ret += m_TotalRanges_minus;
    //     }
    //     if ( m_IsSingleStrand && (m_MoreBefore || m_MoreAfter) ) {
    //         _ASSERT(!m_Ranges.empty());
    //         if ( x_IncludesPlus(m_Ranges.front().second) ) {
    //             if ( (flags & eStrandPlus)  ||
    //                 x_IncludesMinus(m_Ranges.front().second) ) {
    //                 if ( m_MoreBefore ) {
    //                     ret.SetFrom(TRange::GetWholeFrom());
    //                 }
    //                 if ( m_MoreAfter ) {
    //                     ret.SetTo(TRange::GetWholeTo());
    //                 }
    //             }
    //         }
    //         else {
    //             if ( (flags & eStrandMinus) ) {
    //                 if ( m_MoreAfter ) {
    //                     ret.SetFrom(TRange::GetWholeFrom());
    //                 }
    //                 if ( m_MoreBefore ) {
    //                     ret.SetTo(TRange::GetWholeTo());
    //                 }
    //             }
    //         }
    //     }
    //     return ret;
    // }
    // ```
    /// `GetOverlappingRange()` (`flags` `eStrandAny`).
    pub(super) fn overlapping(&self) -> Range {
        if self.circular {
            return Range::WHOLE;
        }
        let mut ret = self.plus;
        ret.combine(self.minus);
        if self.single_strand && (self.more_before || self.more_after) {
            let first = self.ranges[0].1;
            if includes_plus(first) {
                if self.more_before {
                    ret.from = 0;
                }
                if self.more_after {
                    ret.to_open = MAX;
                }
            } else {
                if self.more_after {
                    ret.from = 0;
                }
                if self.more_before {
                    ret.to_open = MAX;
                }
            }
        }
        ret
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:442-445
    // ```c++
    // bool CHandleRange::HasGaps(void) const
    // {
    //     return m_Ranges.size() > 1  ||  m_MoreBefore  ||  m_MoreAfter;
    // }
    // ```
    pub(super) fn has_gaps(&self) -> bool {
        self.ranges.len() > 1 || self.more_before || self.more_after
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:379-417
    // ```c++
    // CHandleRange::GetCircularRangeStart(bool include_origin) const
    // {
    //     _ASSERT(m_IsCircular);
    //     TRange ret = m_TotalRanges_plus;
    //     if ( include_origin ) {
    //         // Adjust start/stop to include cut point
    //         if ( !IsReverse(m_Ranges.front().second) ) {
    //             // Include end
    //             ret.SetTo(TRange::GetWholeTo());
    //         }
    //         else {
    //             // Include start
    //             ret.SetFrom(TRange::GetWholeFrom());
    //         }
    //     }
    //     return ret;
    // }
    // ...
    // CHandleRange::GetCircularRangeEnd(bool include_origin) const
    // {
    //     _ASSERT(m_IsCircular);
    //     TRange ret = m_TotalRanges_minus;
    //     if ( include_origin ) {
    //         // Adjust start/stop to include cut point
    //         if ( !IsReverse(m_Ranges.front().second) ) {
    //             // Include end
    //             ret.SetFrom(TRange::GetWholeFrom());
    //         }
    //         else {
    //             // Include start
    //             ret.SetTo(TRange::GetWholeTo());
    //         }
    //     }
    //     return ret;
    // }
    // ```
    /// `GetCircularRangeStart()` (a circular range).
    pub(super) fn circular_start(&self) -> Range {
        let mut ret = self.plus;
        if !is_reverse(self.ranges[0].1) {
            ret.to_open = MAX;
        } else {
            ret.from = 0;
        }
        ret
    }

    /// `GetCircularRangeEnd()` (a circular range).
    pub(super) fn circular_end(&self) -> Range {
        let mut ret = self.minus;
        if !is_reverse(self.ranges[0].1) {
            ret.from = 0;
        } else {
            ret.to_open = MAX;
        }
        ret
    }
}

/// `CSeq_id_Handle`: the key under which the object manager keeps a Seq-id's ranges
/// (`CSeq_id_Handle::GetHandle`, the tree of the Seq-id's alternative). Two Seq-ids are one handle
/// exactly when `seqid::mask_handle::mask_handle_key` gives them one key (NCBI's identity of every
/// alternative: accessions, names, local and `general` strings, patent numbers without case,
/// versions and releases exactly, GenBank/EMBL/DDBJ accessions alone across the three; `om/` cases
/// `s_i*` of line RG D3-1); a text Seq-id without accession and name, which the tree finds nowhere,
/// is a new handle at each lookup (`Unique`), and so is a `giim` Seq-id here (`Giim`: what NCBI
/// finds for it depends on the `giim` handles alive, `HandleRangeMap::giim_repeats`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum Handle {
    /// The handle of `mask_handle_key`, with its `m_Packed` (`mask_handle_packed`).
    Key { key: Vec<u8>, packed: i64 },
    /// A handle that no other lookup finds.
    Unique(u64),
    /// The lookup of a `giim` Seq-id of the number `id`.
    Giim { id: i32, lookup: u64 },
}

/// The lookups of Seq-ids that no other lookup finds, numbered.
static LOOKUPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1186-1203 (a text Seq-id that
// does not pack: `x_FindStrInfo` looks in `m_ByAcc` when the accession is set, else in `m_ByName`,
// else nowhere; the new info is put in the maps of what is set, so one without accession and name
// is never found again)
// ```c++
//     CSeq_id_Textseq_PlainInfo* info = x_FindStrInfo(id.Which(), tid);
//     CSeq_id_Handle::TVariant variant = 0;
//     if ( !info ) {
//         ...
//         info = new CSeq_id_Textseq_PlainInfo(ref_id, m_Mapper);
//         if ( tid.IsSetAccession() ) {
//             m_ByAcc.insert(TStringMapValue(tid.GetAccession(), info));
//         }
//         if ( tid.IsSetName() ) {
//             m_ByName.insert(TStringMapValue(tid.GetName(), info));
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2941, 2953-2964 (a `giim`
// Seq-id finds the first held handle of its number whose Seq-id is NOT `Equals` to it, else makes
// one: two equal `giim` Seq-ids are two handles, two different ones of a number one handle, unless
// other `giim` handles of the number are alive; `om/` cases `s_i72_giim_same`, `s_i73_giim_db_diff`)
// ```c++
//     TIdMap::const_iterator id_it = m_IdMap.find(gid.GetId());
//     if (id_it == m_IdMap.end())
//         return 0;
//     ITERATE (TGiimList, dbr_it, id_it->second) {
//         CConstRef<CSeq_id> id = (*dbr_it)->GetSeqId();
//         const CGiimport_id& gid2 = id->GetGiim();
//         // Both Db and Release must be equal
//         if ( !gid.Equals(gid2) ) {
//             return *dbr_it;
//         }
//     }
// ...
//     CSeq_id_Info* info = x_FindInfo(gid);
//     if ( !info ) {
//         info = CreateInfo(id);
//         m_IdMap[gid.GetId()].push_back(info);
//     }
//     return CSeq_id_Handle(info);
// ```
impl Handle {
    /// `CSeq_id_Handle::GetHandle(id)`.
    pub(super) fn of(id: &ParsedSeqId) -> Handle {
        let lookup = || LOOKUPS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        match mask_handle_key(id) {
            Ok(Some(key)) => Handle::Key {
                key,
                packed: mask_handle_packed(id),
            },
            Ok(None) => Handle::Unique(lookup()),
            Err(TreeHandle::Giim { id }) => Handle::Giim {
                id,
                lookup: lookup(),
            },
        }
    }

    /// `CSeq_id_Handle::GetGiHandle(gi)`.
    pub(super) fn gi(gi: i64) -> Handle {
        Handle::of(&ParsedSeqId::Gi(gi))
    }

    /// The handle's place in `CSeq_id_Handle::operator<` among handles of another `m_Packed`
    /// (`mask_handle_order`); handles of one `m_Packed` follow the addresses of their infos.
    pub(super) fn order(&self) -> u64 {
        mask_handle_order(match self {
            Handle::Key { packed, .. } => *packed,
            Handle::Unique(_) | Handle::Giim { .. } => 0,
        })
    }
}

/// `CHandleRangeMap::ETransSplicing`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Mode {
    /// `eNoTransSplicing`.
    None,
    /// `eTransSplicing`.
    TransSplicing,
    /// `eCircularRNA`.
    CircularRna,
}

/// `CHandleRangeMap`: a location's ranges by Seq-id, in the order the Seq-ids came; NCBI's map
/// iterates them by `CSeq_id_Handle::operator<` (`HandleRangeMap::ncbi_order`).
#[derive(Default, Debug)]
pub(super) struct HandleRangeMap {
    pub(super) map: Vec<(Handle, HandleRange)>,
}

/// `SAddState`: the Seq-id of the range added before, and the mode.
struct AddState {
    prev: Option<usize>,
    mode: Mode,
}

impl HandleRangeMap {
    /// `AddLocation(loc, mode)` of a fresh map.
    pub(super) fn of_location(loc: V<'_>, mode: Mode) -> HandleRangeMap {
        let mut map = HandleRangeMap::default();
        map.add_location(loc, &mut AddState { prev: None, mode });
        map
    }

    fn index_of(&mut self, handle: Handle) -> usize {
        match self.map.iter().position(|(h, _)| *h == handle) {
            Some(at) => at,
            None => {
                self.map.push((handle, HandleRange::default()));
                self.map.len() - 1
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range_map.cpp:249-281
    // ```c++
    // void CHandleRangeMap::AddRange(const CSeq_id_Handle& h,
    //                                const TRange& range,
    //                                ENa_strand strand,
    //                                SAddState& state)
    // {
    //     CHandleRange& hr = m_LocMap[h];
    //     if ( state.m_TransSplicing == eNoTransSplicing &&
    //          state.m_PrevId && h && state.m_PrevId != h ) {
    //         m_LocMap[state.m_PrevId].m_MoreAfter = true;
    //         hr.m_MoreBefore = true;
    //         if ( m_MasterSeq ) {
    //         ...
    //     }
    //     hr.AddRange(range, strand, false, false, state.m_TransSplicing == eCircularRNA);
    //     state.m_PrevId = h;
    //     state.m_PrevStrand = strand;
    //     state.m_PrevRange = range;
    // }
    // ```
    // The query is a raw Bioseq: its TSE has no master segments (`m_MasterSeq` null).
    fn add(&mut self, handle: Handle, range: Range, strand: Strand, state: &mut AddState) {
        let at = self.index_of(handle);
        if state.mode == Mode::None {
            if let Some(prev) = state.prev.filter(|&prev| prev != at) {
                self.map[prev].1.more_after = true;
                self.map[at].1.more_before = true;
            }
        }
        self.map[at]
            .1
            .add(range, strand, state.mode == Mode::CircularRna);
        state.prev = Some(at);
    }

    /// `AddRange(id, range, strand)` with a state of its own (no Seq-id before it): the rows of an
    /// alignment (annot_object.cpp:591-842).
    pub(super) fn add_range(&mut self, id: &ParsedSeqId, range: Range, strand: Strand) {
        self.add(
            Handle::of(id),
            range,
            strand,
            &mut AddState {
                prev: None,
                mode: Mode::None,
            },
        );
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range_map.cpp:91-205
    // ```c++
    //     switch ( loc.Which() ) {
    //     case CSeq_loc::e_not_set:
    //     case CSeq_loc::e_Null:
    //     {
    //         return;
    //     }
    //     case CSeq_loc::e_Empty:
    //     {
    //         AddRange(loc.GetEmpty(), TRange::GetEmpty(),
    //                  eNa_strand_unknown, state);
    //         return;
    //     }
    //     case CSeq_loc::e_Whole:
    //     {
    //         AddRange(loc.GetWhole(), TRange::GetWhole(),
    //                  eNa_strand_unknown, state);
    //         return;
    //     }
    //     case CSeq_loc::e_Int:
    //     {
    //         const CSeq_interval& i = loc.GetInt();
    //         AddRange(i.GetId(),
    //                  i.GetFrom(),
    //                  i.GetTo(),
    //                  i.IsSetStrand()? i.GetStrand(): eNa_strand_unknown,
    //                  state);
    //         return;
    //     }
    //     case CSeq_loc::e_Pnt:
    //     ...
    //     case CSeq_loc::e_Packed_int:
    //     ...
    //     case CSeq_loc::e_Packed_pnt:
    //     {
    //         // extract each point
    //         const CPacked_seqpnt& pp = loc.GetPacked_pnt();
    //         CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(pp.GetId());
    //         ENa_strand strand =
    //             pp.IsSetStrand()? pp.GetStrand(): eNa_strand_unknown;
    //         ITERATE ( CPacked_seqpnt::TPoints, pi, pp.GetPoints() ) {
    //             AddRange(idh, CRange<TSeqPos>(*pi, *pi), strand, state);
    //         }
    //         return;
    //     }
    //     case CSeq_loc::e_Mix:
    //     {
    //         // extract sub-locations
    //         ITERATE ( CSeq_loc_mix::Tdata, li, loc.GetMix().Get() ) {
    //             AddLocation(**li, state);
    //         }
    //         return;
    //     }
    //     case CSeq_loc::e_Equiv:
    //     {
    //         // extract sub-locations
    //         bool first = true; // allow intron only before the first sub-location
    //         ITERATE ( CSeq_loc_equiv::Tdata, li, loc.GetEquiv().Get() ) {
    //             if ( first ) {
    //                 // the remaining gaps between sub-location aren't introns
    //                 first = false;
    //             }
    //             else {
    //                 // there's no intron between equiv sub-locations
    //                 state.m_PrevId.Reset();
    //             }
    //             AddLocation(**li, state);
    //         }
    //         return;
    //     }
    //     case CSeq_loc::e_Bond:
    //     {
    //         const CSeq_bond& bond = loc.GetBond();
    //         const CSeq_point& pa = bond.GetA();
    //         AddRange(pa.GetId(), pa.GetPoint(), pa.GetPoint(), ..., state);
    //         if ( bond.IsSetB() ) {
    //             const CSeq_point& pb = bond.GetB();
    //             AddRange(pb.GetId(), pb.GetPoint(), pb.GetPoint(), ..., state);
    //         }
    //         return;
    //     }
    //     case CSeq_loc::e_Feat:
    //     {
    //         //### Not implemented (do we need it?)
    //         return;
    //     }
    // ```
    fn add_location(&mut self, loc: V<'_>, state: &mut AddState) {
        let (variant, value) = loc.choice();
        match variant {
            "empty" => self.add(Handle::of(value.seq_id()), Range::EMPTY, UNKNOWN, state),
            "whole" => self.add(Handle::of(value.seq_id()), Range::WHOLE, UNKNOWN, state),
            "int" => self.add_interval(value, state),
            "packed-int" => {
                for interval in value.items() {
                    self.add_interval(interval, state);
                }
            }
            "pnt" => self.add_point(value, state),
            "packed-pnt" => {
                // One handle for all the points (`idh`).
                let handle = Handle::of(value.get("id").expect("mandatory").seq_id());
                let strand = value.get("strand").map_or(UNKNOWN, V::int);
                let points = value.get("points").expect("mandatory");
                for index in 0..points.len() {
                    let point = points.int_at(index).expect("in the list") as u32;
                    self.add(handle.clone(), Range::new(point, point), strand, state);
                }
            }
            "mix" => {
                for part in value.items() {
                    self.add_location(part, state);
                }
            }
            "equiv" => {
                for (index, part) in value.items().enumerate() {
                    if index > 0 {
                        state.prev = None;
                    }
                    self.add_location(part, state);
                }
            }
            "bond" => {
                self.add_point(value.get("a").expect("mandatory"), state);
                if let Some(b) = value.get("b") {
                    self.add_point(b, state);
                }
            }
            // `null`, `feat`.
            _ => {}
        }
    }

    fn add_interval(&mut self, interval: V<'_>, state: &mut AddState) {
        let (id, range, strand) = interval_parts(interval);
        self.add(Handle::of(id), range, strand, state);
    }

    fn add_point(&mut self, point: V<'_>, state: &mut AddState) {
        let (id, range, strand) = point_parts(point);
        self.add(Handle::of(id), range, strand, state);
    }

    /// Whether the maps `maps` hold two lookups of `giim` Seq-ids of one number: NCBI takes them
    /// for one handle or two by the `giim` handles of that number alive at the time (`Handle`),
    /// which come from the query's Seq-ids, the ranges mapped before and the maps that the object
    /// manager keeps from one annotation object to the next (`x_InitFeatKeys` reuses its
    /// `hrmaps` and `key`, seq_annot_info.cpp:726-729).
    pub(super) fn giim_repeats<'m>(maps: impl IntoIterator<Item = &'m HandleRangeMap>) -> bool {
        let mut numbers = Vec::new();
        for map in maps {
            for (handle, _) in &map.map {
                if let Handle::Giim { id, .. } = handle {
                    if numbers.contains(id) {
                        return true;
                    }
                    numbers.push(*id);
                }
            }
        }
        false
    }

    // NCBI reference (598d8ae6): c++/include/objmgr/impl/handle_range_map.hpp:91 (the map's order
    // is `CSeq_id_Handle`'s `operator<`)
    // ```c++
    //     typedef map<CSeq_id_Handle, CHandleRange> TLocMap;
    // ```
    // NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:201-207
    // ```c++
    //     bool operator<  (const CSeq_id_Handle& handle) const
    //         {
    //             // Packed (m_Packed != 0) first:
    //             // zeroes are converted to a highest unsigned value by decrement.
    //             TUintId p1 = INT_ID_TO(TUintId, m_Packed-1);
    //             TUintId p2 = INT_ID_TO(TUintId, handle.m_Packed-1);
    //             return p1 < p2 || (p1 == p2 && m_Info < handle.m_Info);
    //         }
    // ```
    /// The map's entries (indexes into `map`) in NCBI's order, in groups of one `m_Packed`: the
    /// groups in NCBI's order, the handles of a group in the order of their infos' addresses,
    /// which LOSAT does not know (it keeps the order they came in). A gi and the packed accessions
    /// and `general` tags are ordered by their numbers (`om/` cases `s_o*`); the handles of one
    /// number, and every unpacked handle (local, PDB, patent, a text Seq-id with a name or a
    /// release or whose accession does not pack, ...), are one group.
    pub(super) fn ncbi_order(&self) -> Vec<Vec<usize>> {
        let mut entries: Vec<usize> = (0..self.map.len()).collect();
        entries.sort_by_key(|&at| self.map[at].0.order());
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for at in entries {
            match groups.last_mut() {
                Some(group) if self.map[group[0]].0.order() == self.map[at].0.order() => {
                    group.push(at);
                }
                _ => groups.push(vec![at]),
            }
        }
        groups
    }
}

/// A Seq-interval's Seq-id, range `CRange(from, to)` and strand (`unknown` when not set).
pub(super) fn interval_parts(interval: V<'_>) -> (&ParsedSeqId, Range, Strand) {
    let from = interval.get("from").expect("mandatory").int() as u32;
    let to = interval.get("to").expect("mandatory").int() as u32;
    (
        interval.get("id").expect("mandatory").seq_id(),
        Range::new(from, to),
        interval.get("strand").map_or(UNKNOWN, V::int),
    )
}

/// A Seq-point's Seq-id, range `CRange(point, point)` and strand.
pub(super) fn point_parts(point: V<'_>) -> (&ParsedSeqId, Range, Strand) {
    let at = point.get("point").expect("mandatory").int() as u32;
    (
        point.get("id").expect("mandatory").seq_id(),
        Range::new(at, at),
        point.get("strand").map_or(UNKNOWN, V::int),
    )
}

#[cfg(test)]
#[path = "range_tests.rs"]
mod range_tests;
