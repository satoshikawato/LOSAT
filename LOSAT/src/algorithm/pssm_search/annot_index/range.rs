//! The ranges of NCBI's object manager: `COpenRange<TSeqPos>` (`Range`), a Seq-id's ranges of a
//! location (`CHandleRange`, `HandleRange`) and a location's ranges by Seq-id
//! (`CHandleRangeMap`, `HandleRangeMap`), which the annotation index asks for each Seq-id's
//! overlapping range.

use super::view::V;
use crate::seqid::mask_handle::{
    mask_handle_key, mask_handle_order, mask_handle_packed, TreeHandle,
};
use crate::seqid::model::GiimportId;
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
            let first = self.first_strand();
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
    /// The strand of the first range (`m_Ranges.front().second`; the callers' ranges are never
    /// empty).
    fn first_strand(&self) -> Strand {
        self.ranges.first().map_or(UNKNOWN, |&(_, strand)| strand)
    }

    /// `GetCircularRangeStart()` (a circular range).
    pub(super) fn circular_start(&self) -> Range {
        let mut ret = self.plus;
        if !is_reverse(self.first_strand()) {
            ret.to_open = MAX;
        } else {
            ret.from = 0;
        }
        ret
    }

    /// `GetCircularRangeEnd()` (a circular range).
    pub(super) fn circular_end(&self) -> Range {
        let mut ret = self.minus;
        if !is_reverse(self.first_strand()) {
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
/// is a new handle at each lookup (`Unique`); a `giim` Seq-id is the handle that NCBI's inverted
/// `giim` tree gives it (`GiimTree`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum Handle {
    /// The handle of `mask_handle_key`, with its `m_Packed` (`mask_handle_packed`).
    Key { key: Vec<u8>, packed: i64 },
    /// A handle that no other lookup finds.
    Unique(u64),
    /// A `giim` handle of the number `id`.
    Giim { id: i32, handle: GiimHandle },
}

/// The handle of a `giim` lookup (`GiimTree`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum GiimHandle {
    /// The info NCBI's tree gives the lookup (numbered): lookups of one info are one handle.
    Known(u64),
    /// A lookup whose handle depends on `giim` handles whose life LOSAT does not follow
    /// (numbered by the lookup: one handle per lookup, `HandleRangeMap::giim_repeats` refuses a
    /// map where that decides anything).
    Unknown(u64),
}

/// The lookups and infos of Seq-ids that no other lookup finds, numbered.
static LOOKUPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

fn next_lookup() -> u64 {
    LOOKUPS.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1186-1203 (a text Seq-id that
// does not pack: `x_FindStrInfo` looks in `m_ByAcc` when the accession is set, else in `m_ByName`,
// else nowhere; the new info is put in the maps of what is set, so one without accession and name
// is never found again: probe p4 r0261, `om/` case `s_i26_pir_empty`)
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
    /// `CSeq_id_Handle::GetHandle(id)`, a `giim` Seq-id found in `giim`.
    pub(super) fn of(id: &ParsedSeqId, giim: &mut GiimTree) -> Handle {
        match mask_handle_key(id) {
            Ok(Some(key)) => Handle::Key {
                key,
                packed: mask_handle_packed(id),
            },
            Ok(None) => Handle::Unique(next_lookup()),
            Err(TreeHandle::Giim { id: number }) => {
                let ParsedSeqId::Giim(value) = id else {
                    unreachable!("only a giim Seq-id has a tree handle")
                };
                Handle::Giim {
                    id: number,
                    handle: giim.find_or_create(value),
                }
            }
        }
    }

    /// `CSeq_id_Handle::GetHandle(id)` of a lookup whose `giim` handles LOSAT does not follow
    /// (`GiimTree::forget`).
    pub(super) fn of_unfollowed(id: &ParsedSeqId, giim: &mut GiimTree) -> Handle {
        if let ParsedSeqId::Giim(value) = id {
            giim.forget(value.id);
        }
        Handle::of(id, giim)
    }

    /// `CSeq_id_Handle::GetGiHandle(gi)`.
    pub(super) fn gi(gi: i64) -> Handle {
        match mask_handle_key(&ParsedSeqId::Gi(gi)) {
            Ok(Some(key)) => Handle::Key {
                key,
                packed: mask_handle_packed(&ParsedSeqId::Gi(gi)),
            },
            _ => unreachable!("a gi is a key"),
        }
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

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_handle.cpp:89-92 and
// seq_id_tree.cpp:220-236 (an info leaves its tree when its last handle goes:
// `CSeq_id_Info::x_RemoveLastLock` -> `DropInfo` -> `x_Unindex`)
// ```c++
// void CSeq_id_Info::x_RemoveLastLock(void) const
// {
//     GetTree().DropInfo(this);
// }
// ```
// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_info.cpp:300-307 (the query's Seq-ids are
// looked up when its Bioseq joins the scope, before its annotations; their handles live as long
// as the scope)
// ```c++
//     if ( obj.IsSetId() ) {
//         ITERATE ( TObject::TId, it, obj.GetId() ) {
//             m_Id.push_back(CSeq_id_Handle::GetHandle(**it));
//         }
//     }
//     if ( obj.IsSetAnnot() ) {
//         x_SetAnnot();
//     }
// ```
/// NCBI's `giim` tree (`CSeq_id_Giim_Tree`, quoted at `Handle`) during the index. A lookup gets
/// the first live handle of its number whose Seq-id is not `Equals` to it, else a new handle; so
/// the live handles of a number are all made of one value (the first's), a lookup of another value
/// gets the first handle, and a lookup of that value a new one. What is live: the query's Seq-ids
/// (for the whole run, made first: `held`), and within one annotation object every handle its
/// maps took (NCBI builds all of an object's maps before it maps them). Which handles of an
/// earlier object are still alive (its mapped ranges, the `key` and `hrmaps` that
/// `x_InitFeatKeys` keeps from one feature to the next, seq_annot_info.cpp:726-729), and the
/// Seq-ids of the subjects of a `-subject` run (read before the PSSM, blast_args.cpp:2547-2557),
/// LOSAT does not follow: a number looked up by an earlier object, by a Seq-table's columns or a
/// Std-seg's location iterator, and every number in a `-subject` run, is `Unknown`.
#[derive(Debug, Default)]
pub(super) struct GiimTree {
    /// A `-subject` run: every lookup is `Unknown`.
    blind: bool,
    /// The first handle of each number the query's Seq-ids hold, with the Seq-id it was made of.
    held: Vec<(GiimportId, u64)>,
    /// The numbers whose live handles LOSAT does not know.
    unknown: Vec<i32>,
    /// The handles made in the current annotation object, in the order made.
    object: Vec<(GiimportId, u64)>,
}

impl GiimTree {
    /// The tree when the annotations are indexed: after the query's Seq-ids `query_ids`, in a run
    /// with `-subject` sequences or not.
    pub(super) fn new(query_ids: &[ParsedSeqId], subjects: bool) -> GiimTree {
        let mut tree = GiimTree {
            blind: subjects,
            ..GiimTree::default()
        };
        for id in query_ids {
            if let ParsedSeqId::Giim(value) = id {
                if !tree.held.iter().any(|(first, _)| first.id == value.id) {
                    tree.held.push((value.clone(), next_lookup()));
                }
            }
        }
        tree
    }

    /// `CSeq_id_Giim_Tree::FindOrCreate(id)`.
    fn find_or_create(&mut self, value: &GiimportId) -> GiimHandle {
        if self.blind || self.unknown.contains(&value.id) {
            return GiimHandle::Unknown(next_lookup());
        }
        if let Some((first, info)) = self.held.iter().find(|(first, _)| first.id == value.id) {
            return GiimHandle::Known(if first != value { *info } else { next_lookup() });
        }
        match self
            .object
            .iter()
            .find(|(made, _)| made.id == value.id && made != value)
        {
            Some((_, info)) => GiimHandle::Known(*info),
            None => {
                let info = next_lookup();
                self.object.push((value.clone(), info));
                GiimHandle::Known(info)
            }
        }
    }

    /// The handles of the number `id` are no longer followed (a lookup of NCBI's that LOSAT does
    /// not make, or whose handle's life it does not know).
    fn forget(&mut self, id: i32) {
        if !self.unknown.contains(&id) && !self.held.iter().any(|(first, _)| first.id == id) {
            self.unknown.push(id);
        }
    }

    /// The end of an annotation object: the handles it made may live on or not.
    pub(super) fn end_object(&mut self) {
        for (made, _) in std::mem::take(&mut self.object) {
            self.forget(made.id);
        }
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

/// `SAddState`: the Seq-id of the range added before, and the mode; with the `giim` tree of the
/// lookups.
struct AddState<'g> {
    prev: Option<usize>,
    mode: Mode,
    giim: &'g mut GiimTree,
}

impl HandleRangeMap {
    /// `AddLocation(loc, mode)` of a fresh map.
    pub(super) fn of_location(loc: V<'_>, mode: Mode, giim: &mut GiimTree) -> HandleRangeMap {
        let mut map = HandleRangeMap::default();
        map.add_location(
            loc,
            &mut AddState {
                prev: None,
                mode,
                giim,
            },
        );
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
    fn add(&mut self, handle: Handle, range: Range, strand: Strand, state: &mut AddState<'_>) {
        let at = self.index_of(handle);
        if state.mode == Mode::None {
            if let Some(prev) = state.prev.filter(|&prev| prev != at) {
                if let Some((_, hr)) = self.map.get_mut(prev) {
                    hr.more_after = true;
                }
                if let Some((_, hr)) = self.map.get_mut(at) {
                    hr.more_before = true;
                }
            }
        }
        if let Some((_, hr)) = self.map.get_mut(at) {
            hr.add(range, strand, state.mode == Mode::CircularRna);
        }
        state.prev = Some(at);
    }

    /// `AddRange(id, range, strand)` with a state of its own (no Seq-id before it): the rows of an
    /// alignment (annot_object.cpp:591-842).
    pub(super) fn add_range(
        &mut self,
        id: &ParsedSeqId,
        range: Range,
        strand: Strand,
        giim: &mut GiimTree,
    ) {
        let handle = Handle::of(id, giim);
        self.add_handle(handle, range, strand, giim);
    }

    /// `AddRange(handle, range, strand)` with a state of its own.
    pub(super) fn add_handle(
        &mut self,
        handle: Handle,
        range: Range,
        strand: Strand,
        giim: &mut GiimTree,
    ) {
        self.add(
            handle,
            range,
            strand,
            &mut AddState {
                prev: None,
                mode: Mode::None,
                giim,
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
    fn add_location(&mut self, loc: V<'_>, state: &mut AddState<'_>) {
        let (variant, value) = loc.choice();
        match variant {
            "empty" => {
                let handle = Handle::of(value.seq_id(), state.giim);
                self.add(handle, Range::EMPTY, UNKNOWN, state)
            }
            "whole" => {
                let handle = Handle::of(value.seq_id(), state.giim);
                self.add(handle, Range::WHOLE, UNKNOWN, state)
            }
            "int" => self.add_interval(value, state),
            "packed-int" => {
                for interval in value.items() {
                    self.add_interval(interval, state);
                }
            }
            "pnt" => self.add_point(value, state),
            "packed-pnt" => {
                // One handle for all the points (`idh`).
                let handle = Handle::of(value.field("id").seq_id(), state.giim);
                let strand = value.get("strand").map_or(UNKNOWN, V::int);
                for point in value.field("points").ints() {
                    let point = point as u32;
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
                self.add_point(value.field("a"), state);
                if let Some(b) = value.get("b") {
                    self.add_point(b, state);
                }
            }
            // `null`, `feat`.
            _ => {}
        }
    }

    fn add_interval(&mut self, interval: V<'_>, state: &mut AddState<'_>) {
        let (id, range, strand) = interval_parts(interval);
        let handle = Handle::of(id, state.giim);
        self.add(handle, range, strand, state);
    }

    fn add_point(&mut self, point: V<'_>, state: &mut AddState<'_>) {
        let (id, range, strand) = point_parts(point);
        let handle = Handle::of(id, state.giim);
        self.add(handle, range, strand, state);
    }

    /// Whether the maps `maps` hold two `giim` handles of one number, one of them `Unknown`: NCBI
    /// takes their Seq-ids for one handle or two by `giim` handles whose life LOSAT does not
    /// follow (`GiimTree`).
    pub(super) fn giim_repeats<'m>(maps: impl IntoIterator<Item = &'m HandleRangeMap>) -> bool {
        let mut handles: Vec<(i32, GiimHandle)> = Vec::new();
        for map in maps {
            for (handle, _) in &map.map {
                if let Handle::Giim { id, handle } = handle {
                    if !handles.contains(&(*id, *handle)) {
                        handles.push((*id, *handle));
                    }
                }
            }
        }
        handles.iter().any(|&(id, handle)| {
            matches!(handle, GiimHandle::Unknown(_))
                && handles.iter().filter(|&&(other, _)| other == id).count() > 1
        })
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
    /// The map's entries in NCBI's order, in groups of one `m_Packed`: the
    /// groups in NCBI's order, the handles of a group in the order of their infos' addresses,
    /// which LOSAT does not know (it keeps the order they came in). A gi and the packed accessions
    /// and `general` tags are ordered by their numbers (`om/` cases `s_o*`); the handles of one
    /// number, and every unpacked handle (local, PDB, patent, a text Seq-id with a name or a
    /// release or whose accession does not pack, ...), are one group.
    pub(super) fn ncbi_order(&self) -> Vec<Vec<&(Handle, HandleRange)>> {
        let mut entries: Vec<(u64, &(Handle, HandleRange))> = self
            .map
            .iter()
            .map(|entry| (entry.0.order(), entry))
            .collect();
        entries.sort_by_key(|&(order, _)| order);
        entries
            .chunk_by(|a, b| a.0 == b.0)
            .map(|group| group.iter().map(|&(_, entry)| entry).collect())
            .collect()
    }
}

/// A Seq-interval's Seq-id, range `CRange(from, to)` and strand (`unknown` when not set).
pub(super) fn interval_parts(interval: V<'_>) -> (&ParsedSeqId, Range, Strand) {
    let from = interval.field("from").int() as u32;
    let to = interval.field("to").int() as u32;
    (
        interval.field("id").seq_id(),
        Range::new(from, to),
        interval.get("strand").map_or(UNKNOWN, V::int),
    )
}

/// A Seq-point's Seq-id, range `CRange(point, point)` and strand.
pub(super) fn point_parts(point: V<'_>) -> (&ParsedSeqId, Range, Strand) {
    let at = point.field("point").int() as u32;
    (
        point.field("id").seq_id(),
        Range::new(at, at),
        point.get("strand").map_or(UNKNOWN, V::int),
    )
}

#[cfg(test)]
#[path = "range_tests.rs"]
mod range_tests;
