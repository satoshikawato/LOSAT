//! The masked ranges of the sequences, by algorithm and sequence
//! (`CMaskedRangeSet`, `c++/src/app/blastdb/masked_range_set.cpp`), with the
//! combination of several Seq-locs of one sequence (`CSeq_loc::Add` with
//! `fMerge_All | fSort`, Seq_loc.cpp:4283-4693, 5060-5084).
//!
//! A mask finds a sequence the way NCBI's `map<CSeq_id_Handle, CConstRef<CSeq_loc>>` does: through
//! the handle of its Seq-id, and a sequence looks its Seq-ids up in the order of its `GetId()`
//! (`seqid::mask_handle` has the identity of handles: a mask for `embl|AY100001.1` masks
//! `gb|AY100001.1|`, a local string id matches without case, `local str "5"` does not match
//! `lcl|5`, ...). The ranges that several Seq-ids of one sequence find, and the Seq-locs that
//! one handle gets, are combined with `CSeq_loc::Add`, ported here with NCBI's `CRange` arithmetic
//! (`Range`): empty ranges (`to < from`, and `to` = 4294967295 whose open end wraps to 0) and the
//! whole range (`0..4294967295`) take part in the sort and the merge exactly as in NCBI, an empty
//! result becomes an `empty` Seq-loc that the final `CSeq_loc_CI` skips, a whole one a `whole`
//! Seq-loc that the writer's bounds check then refuses as NCBI's does. The sort is libstdc++'s
//! `std::sort` (`stl_sort`), because the order in which it leaves equivalent empty ranges can
//! change the merge.
//!
//! Handles of `giim` Seq-ids are kept in a list per number, as NCBI's `CSeq_id_Giim_Tree`
//! (seq_id_tree.hpp:966-967, seq_id_tree.cpp:2926-2964): a Seq-id finds the first stored handle of
//! its number whose Giimport-id is *not* `Equals` to it, else a new handle is appended, and a
//! handle that nothing holds any more leaves the list (`x_Unindex`, seq_id_tree.cpp:2967-2982;
//! seq_id_handle.cpp:89-92). The handles of the masks (the keys of `m_Values`) stay for the run;
//! the ones `CSeq_loc_CI` (`x_SetId`, Seq_loc.cpp:1422-1427) and `GetBestSynonym` make inside a
//! `CSeq_loc::Add` go when it returns. So which mask a `giim` Seq-id finds follows the order of
//! the mask files, and is computed (`GiimTree`).
//!
//! LOSAT refuses (`WriteDbError::Unsupported`, "losat makeblastdb: ...", exit status 1) only what
//! NCBI decides by memory addresses. `x_MergeAndSort` groups the ranges of a combination in
//! `map<CSeq_id_Handle_Wrapper, TRanges>` (Seq_loc.cpp:4651-4686), ordered by
//! `CSeq_id_Handle::operator<`: `m_Packed - 1` as unsigned, then the address of the handle's info
//! (seq_id_handle.hpp:201-207); the plus-strand groups are written before the minus-strand ones
//! (`x_RangesToSeq_loc` once per map, 4688-4691), and ranges of different handles are never merged.
//! Two different handles with the same `m_Packed` (`gi|5` and `gb|AY000005.1|`, 5; two unpacked
//! handles such as a local id and a `general` string tag, or two `giim` handles, 0) are ordered by
//! address, so LOSAT refuses a combination that writes ranges of two such handles on one strand
//! in an order that matters (the ranges written of one are not those of the other repeated,
//! `note_ties`). One strand each is decided, as plus before minus (review probe p4).
//!
//! `giim` handles need no rule of their own: the first mask Seq-id of a number makes the first
//! handle, which stays (a key of `m_Values`); a Seq-id that is not `Equals` to it finds that
//! handle, one that is finds no handle and makes a new one, made from an equal Seq-id. So every
//! `giim` handle of a number is made from Seq-ids `Equals` to the first one, the lookups of a
//! later `CSeq_loc::Add` give the same groups in any order of the ranges (each range of that
//! Seq-id a group of its own, the others in the first handle's group), and only the address order
//! of those groups is left, which the rule above checks at each combination
//! (`tests/seqdb/mask_ids` n41, n42 refused; n65 decided; n66 refused).
use std::collections::{BTreeMap, HashMap};

use crate::seqid::mask_handle::{
    mask_handle_key, mask_handle_order, mask_handle_packed, TreeHandle,
};
use crate::seqid::model::GiimportId;
use crate::seqid::ParsedSeqId;
use crate::writedb::mask_info::MaskPart;
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

// NCBI reference (598d8ae6): c++/include/util/range.hpp:148-151,178-191,256-308,326-334,358-367
// ```c++
//     bool Empty(void) const
//         {
//             return GetToOpen() <= GetFrom();
//         }
//     TThisType& SetTo(position_type to)
//         {
//             return SetToOpen(to+1);
//         }
//     TThisType& Set(position_type from, position_type to)
//         {
//             return SetFrom(from).SetTo(to);
//         }
//     static position_type GetWholeFrom(void) { return GetPositionMin(); }
//     static position_type GetWholeToOpen(void) { return GetPositionMax(); }
//     bool IsWhole(void) const { return IsWholeFrom() && IsWholeTo(); }
//     static position_type GetEmptyFrom(void) { return GetPositionMax(); }
//     static position_type GetEmptyToOpen(void) { return GetPositionMax(); }
//     TThisType IntersectionWith(const TThisType& r) const
//         {
//             return TThisType(max(GetFrom(), r.GetFrom()),
//                              min(GetToOpen(), r.GetToOpen()));
//         }
//     bool IntersectingWith(const TThisType& r) const
//         {
//             return IntersectionWith(r).NotEmpty();
//         }
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
// ```
/// NCBI's `CRange<TSeqPos>`: `from` and the open end `to_open`, in `TSeqPos` (unsigned 32-bit,
/// wrapping) arithmetic. `GetTo()` is `to_open - 1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Range {
    from: u32,
    to_open: u32,
}

impl Range {
    /// `TRange::GetEmpty()`.
    const EMPTY: Range = Range {
        from: u32::MAX,
        to_open: u32::MAX,
    };
    /// `TRange::GetWhole()`.
    const WHOLE: Range = Range {
        from: 0,
        to_open: u32::MAX,
    };

    /// `Set(from, to)`, as `CSeq_loc_CI` makes the range of an interval or a point.
    fn new(from: u32, to: u32) -> Self {
        Range {
            from,
            to_open: to.wrapping_add(1),
        }
    }

    fn to(self) -> u32 {
        self.to_open.wrapping_sub(1)
    }

    fn is_empty(self) -> bool {
        self.to_open <= self.from
    }

    fn is_whole(self) -> bool {
        self == Self::WHOLE
    }

    fn intersecting(self, other: Range) -> bool {
        self.from.max(other.from) < self.to_open.min(other.to_open)
    }

    fn combine(&mut self, other: Range) {
        if !other.is_empty() {
            if self.is_empty() {
                *self = other;
            } else {
                self.from = self.from.min(other.from);
                self.to_open = self.to_open.max(other.to_open);
            }
        }
    }
}

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
/// `CRange_Less` (or `CRange_ReverseLess` with `reverse`): whole ranges first, then empty ones,
/// then the others.
fn range_less(rg1: &Range, rg2: &Range, reverse: bool) -> bool {
    if rg1.is_whole() {
        return !rg2.is_whole();
    }
    if rg1.is_empty() {
        return !rg2.is_empty() && !rg2.is_whole();
    }
    if rg2.is_whole() || rg2.is_empty() {
        return false;
    }
    if reverse {
        if rg1.to() != rg2.to() {
            return rg1.to() > rg2.to();
        }
        rg1.from < rg2.from
    } else {
        // `COpenRange::operator<`.
        rg1.from < rg2.from || (rg1.from == rg2.from && rg1.to_open < rg2.to_open)
    }
}

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
// `fMerge_All` sets all three merge bits; the strand of both ranges of a comparison is the default
// strand of their group, whose reverse-ness is `reverse`. Fuzz is not kept (only the positions
// reach the database).
fn merge_ranges(rg1: &mut Range, rg2: Range, reverse: bool) -> bool {
    // Contained.
    if rg1.from <= rg2.from && rg1.to() >= rg2.to() {
        return true;
    }
    if rg1.from >= rg2.from && rg1.to() <= rg2.to() {
        *rg1 = rg2;
        return true;
    }
    // Overlapping: `rg1 += rg2` (`CombineWith`).
    if rg1.intersecting(rg2) {
        rg1.combine(rg2);
        return true;
    }
    // Abutting: `CopyTo` sets `to_open` to `rg2.GetTo() + 1`, `CopyFrom` sets `from`.
    if !reverse {
        if rg1.to_open == rg2.from {
            rg1.to_open = rg2.to_open;
            return true;
        }
    } else if rg1.from == rg2.to_open {
        rg1.from = rg2.from;
        return true;
    }
    false
}

// NCBI reference (598d8ae6): /usr/include/c++/8/bits/stl_algo.h (libstdc++ of GCC 8.5, which
// built makeblastdb 2.17.0: `strings makeblastdb` = "GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-26)")
// ```c++
// enum { _S_threshold = 16 };
// __sort(__first, __last, __comp) {
//   if (__first != __last) {
//     std::__introsort_loop(__first, __last, std::__lg(__last - __first) * 2, __comp);
//     std::__final_insertion_sort(__first, __last, __comp);
//   }
// }
// __introsort_loop(__first, __last, __depth_limit, __comp) {
//   while (__last - __first > int(_S_threshold)) {
//     if (__depth_limit == 0) {
//       std::__partial_sort(__first, __last, __last, __comp);   // __heap_select + __sort_heap
//       return;
//     }
//     --__depth_limit;
//     _RandomAccessIterator __cut = std::__unguarded_partition_pivot(__first, __last, __comp);
//     std::__introsort_loop(__cut, __last, __depth_limit, __comp);
//     __last = __cut;
//   }
// }
// __unguarded_partition_pivot(__first, __last, __comp) {
//   _RandomAccessIterator __mid = __first + (__last - __first) / 2;
//   std::__move_median_to_first(__first, __first + 1, __mid, __last - 1, __comp);
//   return std::__unguarded_partition(__first + 1, __last, __first, __comp);
// }
// __move_median_to_first(__result, __a, __b, __c, __comp) {
//   if (__comp(__a, __b)) {
//     if (__comp(__b, __c)) std::iter_swap(__result, __b);
//     else if (__comp(__a, __c)) std::iter_swap(__result, __c);
//     else std::iter_swap(__result, __a);
//   }
//   else if (__comp(__a, __c)) std::iter_swap(__result, __a);
//   else if (__comp(__b, __c)) std::iter_swap(__result, __c);
//   else std::iter_swap(__result, __b);
// }
// __unguarded_partition(__first, __last, __pivot, __comp) {
//   while (true) {
//     while (__comp(__first, __pivot)) ++__first;
//     --__last;
//     while (__comp(__pivot, __last)) --__last;
//     if (!(__first < __last)) return __first;
//     std::iter_swap(__first, __last);
//     ++__first;
//   }
// }
// __final_insertion_sort(__first, __last, __comp) {
//   if (__last - __first > int(_S_threshold)) {
//     std::__insertion_sort(__first, __first + int(_S_threshold), __comp);
//     std::__unguarded_insertion_sort(__first + int(_S_threshold), __last, __comp);
//   }
//   else std::__insertion_sort(__first, __last, __comp);
// }
// __insertion_sort(__first, __last, __comp) {
//   if (__first == __last) return;
//   for (__i = __first + 1; __i != __last; ++__i) {
//     if (__comp(__i, __first)) {
//       __val = *__i; move_backward(__first, __i, __i + 1); *__first = __val;
//     }
//     else std::__unguarded_linear_insert(__i, __comp);
//   }
// }
// __unguarded_linear_insert(__last, __comp) {
//   __val = *__last; __next = __last; --__next;
//   while (__comp(__val, __next)) { *__last = *__next; __last = __next; --__next; }
//   *__last = __val;
// }
// ```
// NCBI reference (598d8ae6): /usr/include/c++/8/bits/stl_heap.h (__make_heap, __adjust_heap,
// __push_heap, __pop_heap, __sort_heap), and stl_algo.h __heap_select
// ```c++
// __heap_select(__first, __middle, __last, __comp) {
//   std::__make_heap(__first, __middle, __comp);
//   for (__i = __middle; __i < __last; ++__i)
//     if (__comp(__i, __first)) std::__pop_heap(__first, __middle, __i, __comp);
// }
// __make_heap(__first, __last, __comp) {
//   if (__last - __first < 2) return;
//   __len = __last - __first; __parent = (__len - 2) / 2;
//   while (true) {
//     __value = *(__first + __parent);
//     std::__adjust_heap(__first, __parent, __len, __value, __comp);
//     if (__parent == 0) return;
//     __parent--;
//   }
// }
// __adjust_heap(__first, __holeIndex, __len, __value, __comp) {
//   __topIndex = __holeIndex; __secondChild = __holeIndex;
//   while (__secondChild < (__len - 1) / 2) {
//     __secondChild = 2 * (__secondChild + 1);
//     if (__comp(__first + __secondChild, __first + (__secondChild - 1))) __secondChild--;
//     *(__first + __holeIndex) = *(__first + __secondChild);
//     __holeIndex = __secondChild;
//   }
//   if ((__len & 1) == 0 && __secondChild == (__len - 2) / 2) {
//     __secondChild = 2 * (__secondChild + 1);
//     *(__first + __holeIndex) = *(__first + (__secondChild - 1));
//     __holeIndex = __secondChild - 1;
//   }
//   std::__push_heap(__first, __holeIndex, __topIndex, __value, __comp);
// }
// __push_heap(__first, __holeIndex, __topIndex, __value, __comp) {
//   __parent = (__holeIndex - 1) / 2;
//   while (__holeIndex > __topIndex && __comp(__first + __parent, __value)) {
//     *(__first + __holeIndex) = *(__first + __parent);
//     __holeIndex = __parent; __parent = (__holeIndex - 1) / 2;
//   }
//   *(__first + __holeIndex) = __value;
// }
// __pop_heap(__first, __last, __result, __comp) {
//   __value = *__result; *__result = *__first;
//   std::__adjust_heap(__first, 0, __last - __first, __value, __comp);
// }
// __sort_heap(__first, __last, __comp) {
//   while (__last - __first > 1) { --__last; std::__pop_heap(__first, __last, __last, __comp); }
// }
// ```
/// `std::sort` of libstdc++ (introsort): the same permutation as NCBI's `sort(..)`, also of
/// elements that compare equivalent and differ (two empty ranges).
fn stl_sort<T: Copy>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    if v.is_empty() {
        return;
    }
    let depth = 2 * (usize::BITS - 1 - v.len().leading_zeros()) as usize;
    introsort_loop(v, 0, v.len(), depth, less);
    // `__final_insertion_sort`.
    if v.len() > STL_THRESHOLD {
        insertion_sort(&mut v[..STL_THRESHOLD], less);
        for i in STL_THRESHOLD..v.len() {
            unguarded_linear_insert(v, i, less);
        }
    } else {
        insertion_sort(v, less);
    }
}

const STL_THRESHOLD: usize = 16;

fn introsort_loop<T: Copy>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    while last - first > STL_THRESHOLD {
        if depth == 0 {
            heap_sort(&mut v[first..last], less);
            return;
        }
        depth -= 1;
        let mid = first + (last - first) / 2;
        move_median_to_first(v, first, first + 1, mid, last - 1, less);
        let cut = unguarded_partition(v, first + 1, last, first, less);
        introsort_loop(v, cut, last, depth, less);
        last = cut;
    }
}

fn move_median_to_first<T: Copy>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    let pick = if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            b
        } else if less(&v[a], &v[c]) {
            c
        } else {
            a
        }
    } else if less(&v[a], &v[c]) {
        a
    } else if less(&v[b], &v[c]) {
        c
    } else {
        b
    };
    v.swap(result, pick);
}

fn unguarded_partition<T: Copy>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    less: &impl Fn(&T, &T) -> bool,
) -> usize {
    loop {
        while less(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while less(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn insertion_sort<T: Copy>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    for i in 1..v.len() {
        if less(&v[i], &v[0]) {
            let val = v[i];
            v.copy_within(0..i, 1);
            v[0] = val;
        } else {
            unguarded_linear_insert(v, i, less);
        }
    }
}

fn unguarded_linear_insert<T: Copy>(v: &mut [T], mut last: usize, less: &impl Fn(&T, &T) -> bool) {
    let val = v[last];
    while less(&val, &v[last - 1]) {
        v[last] = v[last - 1];
        last -= 1;
    }
    v[last] = val;
}

/// `__partial_sort(first, last, last)`: `__heap_select` (a heap of all of it) and `__sort_heap`.
fn heap_sort<T: Copy>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    let len = v.len();
    if len >= 2 {
        let mut parent = (len - 2) / 2;
        loop {
            let value = v[parent];
            adjust_heap(v, parent, len, value, less);
            if parent == 0 {
                break;
            }
            parent -= 1;
        }
    }
    let mut last = len;
    while last > 1 {
        last -= 1;
        let value = v[last];
        v[last] = v[0];
        adjust_heap(&mut v[..last], 0, last, value, less);
    }
}

fn adjust_heap<T: Copy>(
    v: &mut [T],
    mut hole: usize,
    len: usize,
    value: T,
    less: &impl Fn(&T, &T) -> bool,
) {
    let top = hole;
    let mut second = hole;
    while len > 0 && second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if less(&v[second], &v[second - 1]) {
            second -= 1;
        }
        v[hole] = v[second];
        hole = second;
    }
    if len & 1 == 0 && len >= 2 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[hole] = v[second - 1];
        hole = second - 1;
    }
    // `__push_heap`.
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&v[parent], &value) {
            break;
        }
        v[hole] = v[parent];
        hole = parent;
    }
    v[hole] = value;
}

/// What `CSeq_id_Handle::GetHandle` makes of a Seq-id.
#[derive(Clone, Debug)]
enum IdHandle {
    /// A handle that a pure function knows (`mask_handle_key`), and its `m_Packed`.
    Key { key: Vec<u8>, packed: i64 },
    /// A Textseq-id with neither accession nor name: a new handle on every call, which nothing
    /// finds again (seq_id_tree.cpp:1187-1198).
    Fresh,
    /// A `giim` Seq-id: found in the [`GiimTree`].
    Giim(GiimportId),
}

/// A Seq-id of a stored Seq-loc, with what its handle is.
#[derive(Clone, Debug)]
struct LocId {
    seq_id: ParsedSeqId,
    handle: IdHandle,
}

impl LocId {
    fn new(seq_id: ParsedSeqId) -> Self {
        let handle = match mask_handle_key(&seq_id) {
            Ok(Some(key)) => IdHandle::Key {
                key,
                packed: mask_handle_packed(&seq_id),
            },
            Ok(None) => IdHandle::Fresh,
            Err(TreeHandle::Giim { .. }) => match &seq_id {
                ParsedSeqId::Giim(giim) => IdHandle::Giim(giim.clone()),
                _ => unreachable!("only a giim Seq-id has a tree handle"),
            },
        };
        LocId { seq_id, handle }
    }
}

/// One range of a Seq-loc as `CSeq_loc_CI` visits it.
#[derive(Clone, Copy, Debug)]
struct Part {
    /// The Seq-id (an index into the Seq-loc's `ids`); `None` for a `null` Seq-loc.
    id: Option<u32>,
    /// `None` for an `empty` or `null` Seq-loc: skipped by `CSeq_loc_CI` (`eEmpty_Skip`, the
    /// final iteration), `TRange::GetEmpty()` with `eEmpty_Allow` (inside `CSeq_loc::Add`).
    range: Option<Range>,
    /// `IsReverse(strand)`: the strand is used for nothing else.
    reverse: bool,
}

/// A Seq-loc as `CMaskedRangeSet` holds it: the one given, or the result of `CSeq_loc::Add`.
#[derive(Clone, Debug, Default)]
struct Loc {
    ids: Vec<LocId>,
    parts: Vec<Part>,
    /// Set when the ranges NCBI writes of this Seq-loc depend on memory addresses.
    written: Option<String>,
}

impl Loc {
    /// The Seq-loc given for `id`: the parts `CSeq_loc_CI` visits (`MaskPart`).
    fn given(id: ParsedSeqId, parts: Vec<MaskPart>) -> Self {
        let parts = parts
            .into_iter()
            .map(|part| match part {
                MaskPart::Interval(iv) => Part {
                    id: Some(0),
                    range: Some(Range::new(iv.from, iv.to)),
                    reverse: is_reverse(iv.strand),
                },
                MaskPart::Empty => Part {
                    id: Some(0),
                    range: None,
                    reverse: false,
                },
                MaskPart::Null => Part {
                    id: None,
                    range: None,
                    reverse: false,
                },
            })
            .collect();
        Loc {
            ids: vec![LocId::new(id)],
            parts,
            written: None,
        }
    }

    /// The offsets that `GetRanges` writes (`CSeq_loc_CI`, `eEmpty_Skip`; `GetToOpen`).
    fn offsets(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        self.parts
            .iter()
            .filter_map(|p| p.range.map(|r| (r.from, r.to_open)))
    }
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

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2964
// ```c++
// CSeq_id_Info* CSeq_id_Giim_Tree::x_FindInfo(const CGiimport_id& gid) const
// {
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
//     // Not found
//     return 0;
// }
// ...
// CSeq_id_Handle CSeq_id_Giim_Tree::FindOrCreate(const CSeq_id& id)
// {
//     ...
//     CSeq_id_Info* info = x_FindInfo(gid);
//     if ( !info ) {
//         info = CreateInfo(id);
//         m_IdMap[gid.GetId()].push_back(info);
//     }
//     return CSeq_id_Handle(info);
// }
// ```
/// NCBI's `CSeq_id_Giim_Tree`: the live handles of each `giim` number, in the order they were
/// made, each with the Giimport-id it was made from (`GetSeqId()` of the info).
#[derive(Debug, Default)]
struct GiimTree {
    lists: HashMap<i32, Vec<(usize, GiimportId)>>,
    made: usize,
}

impl GiimTree {
    /// `x_FindInfo`: the first live handle whose Giimport-id is not `Equals` to `id`.
    fn find(&self, id: &GiimportId) -> Option<usize> {
        self.lists
            .get(&id.id)?
            .iter()
            .find(|(_, made_from)| made_from != id)
            .map(|(handle, _)| *handle)
    }

    /// `FindOrCreate`: the handle, and whether it is new.
    fn find_or_create(&mut self, id: &GiimportId) -> (usize, bool) {
        if let Some(handle) = self.find(id) {
            return (handle, false);
        }
        self.made += 1;
        self.lists
            .entry(id.id)
            .or_default()
            .push((self.made, id.clone()));
        (self.made, true)
    }

    /// The Giimport-id the handle was made from.
    fn made_from(&self, number: i32, handle: usize) -> &GiimportId {
        let list = &self.lists[&number];
        &list
            .iter()
            .find(|(h, _)| *h == handle)
            .expect("a live handle")
            .1
    }

    /// `x_Unindex` of a handle that nothing holds any more.
    fn drop_handle(&mut self, number: i32, handle: usize) {
        if let Some(list) = self.lists.get_mut(&number) {
            list.retain(|(h, _)| *h != handle);
            if list.is_empty() {
                self.lists.remove(&number);
            }
        }
    }
}

/// The handle of a group of `x_MergeAndSort`.
#[derive(Clone, Debug, PartialEq, Eq)]
enum GroupHandle {
    /// No Seq-id: the null handle (`m_Packed` 0, no info: before the other unpacked handles).
    Null,
    Key(Vec<u8>),
    Giim(usize),
    Fresh(usize),
}

/// One group of `x_MergeAndSort`: a handle, its `m_Packed`, the Seq-id its Seq-locs get
/// (`idh.GetSeqId()`), and its ranges of each strand.
struct Group {
    handle: GroupHandle,
    packed: i64,
    id: Option<LocId>,
    ranges: [Vec<Range>; 2],
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
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1439-1471,1500-1531
// ```c++
// void CSeq_loc_CI_Impl::x_ProcessLocation(const CSeq_loc& loc)
// {
//     switch ( loc.Which() ) {
//     case CSeq_loc::e_not_set:
//     case CSeq_loc::e_Null:
//     case CSeq_loc::e_Empty:
//         {
//             if (m_EmptyFlag == CSeq_loc_CI::eEmpty_Allow) {
//                 SSeq_loc_CI_RangeInfo info;
//                 if (loc.Which() == CSeq_loc::e_Empty) {
//                     x_SetId(info, loc.GetEmpty());
//                 }
//                 else {
//                     info.m_Id.Reset(new CSeq_id);
//                 }
//                 info.m_Range = TRange::GetEmpty();
//                 info.m_Loc = &loc;
//                 m_Ranges.push_back(info);
//             }
//             return;
//         }
//     case CSeq_loc::e_Whole:
//         {
//             SSeq_loc_CI_RangeInfo info;
//             x_SetId(info, loc.GetWhole());
//             info.m_Range = TRange::GetWhole();
//             ...
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4392-4437
// ```c++
// static
// void x_PushRange(CSeq_loc& dst,
//                  const CSeq_id_Handle_Wrapper& idh,
//                  const TRangeWithFuzz& rg,
//                  ENa_strand strand)
// {
//     ...
//     if ( !idh ) {
//         // NULL
//         ... dst.SetMix().Set().push_back(Ref(new CSeq_loc(CSeq_loc::e_Null)));
//         return;
//     }
//     CRef<CSeq_id> id(new CSeq_id());
//     id->Assign(*idh.GetSeqId());
//     if ( rg.IsWhole() ) {
//         ... whole->SetWhole(*id);
//     }
//     else if ( rg.Empty() ) {
//         ... empty->SetEmpty(*id);
//     }
//     else ... (a point or an interval with `strand` unless it is unknown)
// ```
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:335-338
// ```c++
// CSeq_id_Handle CSeq_id_not_set_Tree::FindOrCreate(const CSeq_id& /*id*/)
// {
//     return null;
// }
// ```
// `CSeq_loc_CI` (`eEmpty_Allow`) first makes the handle of every range's Seq-id (`x_SetId`), then
// the loop takes `GetBestSynonym` = `GetHandle` of each again; both sets of handles go when
// `x_MergeAndSort` returns. A `null` Seq-loc has an unset Seq-id, whose handle is the null one:
// its group writes one `null` Seq-loc. A whole range becomes a `whole` Seq-loc and an empty one an
// `empty` Seq-loc, both without strand (the next iteration sees them on the plus strand).
/// `L1->Add(L2, fMerge_All | fSort, NULL)`: the Seq-loc `x_MergeAndSort` writes, with what in
/// it depends on memory addresses (`Loc::written`).
fn add_locs(tree: &mut GiimTree, l1: &Loc, l2: &Loc) -> Loc {
    let parts = || {
        l1.parts
            .iter()
            .map(move |p| (l1, p))
            .chain(l2.parts.iter().map(move |p| (l2, p)))
    };
    let mut made: Vec<(i32, usize)> = Vec::new();
    // The handles of `CSeq_loc_CI` (`x_SetId`).
    for (loc, part) in parts() {
        if let Some(index) = part.id {
            if let IdHandle::Giim(giim) = &loc.ids[index as usize].handle {
                let (handle, new) = tree.find_or_create(giim);
                if new {
                    made.push((giim.id, handle));
                }
            }
        }
    }
    // The groups, in the order their first range comes (`GetBestSynonym`).
    let mut groups: Vec<Group> = Vec::new();
    let mut fresh = 0usize;
    for (loc, part) in parts() {
        let range = part.range.unwrap_or(Range::EMPTY);
        let side = usize::from(part.reverse);
        let (handle, packed, id) = match part.id {
            None => (GroupHandle::Null, 0, None),
            Some(index) => {
                let loc_id = &loc.ids[index as usize];
                match &loc_id.handle {
                    IdHandle::Key { key, packed } => {
                        let same = |g: &Group| matches!(&g.handle, GroupHandle::Key(k) if k == key);
                        if let Some(group) = groups.iter_mut().find(|g| same(g)) {
                            group.ranges[side].push(range);
                            continue;
                        }
                        (GroupHandle::Key(key.clone()), *packed, Some(loc_id.clone()))
                    }
                    IdHandle::Fresh => {
                        fresh += 1;
                        (GroupHandle::Fresh(fresh), 0, Some(loc_id.clone()))
                    }
                    IdHandle::Giim(giim) => {
                        let (handle, new) = tree.find_or_create(giim);
                        if new {
                            made.push((giim.id, handle));
                        }
                        let made_from = tree.made_from(giim.id, handle).clone();
                        let seq_id = ParsedSeqId::Giim(made_from.clone());
                        (
                            GroupHandle::Giim(handle),
                            0,
                            Some(LocId {
                                seq_id,
                                handle: IdHandle::Giim(made_from),
                            }),
                        )
                    }
                }
            }
        };
        match groups.iter_mut().find(|g| g.handle == handle) {
            Some(group) => group.ranges[side].push(range),
            None => {
                let mut ranges = [Vec::new(), Vec::new()];
                ranges[side].push(range);
                groups.push(Group {
                    handle,
                    packed,
                    id,
                    ranges,
                });
            }
        }
    }
    for (number, handle) in made {
        tree.drop_handle(number, handle);
    }
    // `operator<`: `m_Packed - 1` unsigned, then the info's address (the null handle's is 0).
    // Groups that the address orders keep the order of their first range here.
    groups.sort_by_key(|g| (mask_handle_order(g.packed), g.handle != GroupHandle::Null));
    let mut out = Loc::default();
    // The index of each group's Seq-id in `out.ids` (the null group has none).
    let id_index: Vec<Option<u32>> = groups
        .iter()
        .map(|group| {
            group.id.as_ref().map(|id| {
                out.ids.push(id.clone());
                (out.ids.len() - 1) as u32
            })
        })
        .collect();
    // Per strand: the parts each group writes.
    for side in 0..2 {
        let reverse = side == 1;
        let mut written: Vec<(usize, Vec<Part>)> = Vec::new();
        for (index, group) in groups.iter_mut().enumerate() {
            let ranges = &mut group.ranges[side];
            if ranges.is_empty() {
                continue;
            }
            let mut pushed = Vec::new();
            if group.handle == GroupHandle::Null {
                // "All NULLs merged".
                pushed.push(Part {
                    id: None,
                    range: None,
                    reverse: false,
                });
            } else {
                stl_sort(ranges, &|a: &Range, b: &Range| range_less(a, b, reverse));
                let id = id_index[index];
                let mut push = |rg: Range| {
                    pushed.push(if rg.is_whole() {
                        Part {
                            id,
                            range: Some(Range::WHOLE),
                            reverse: false,
                        }
                    } else if rg.is_empty() {
                        Part {
                            id,
                            range: None,
                            reverse: false,
                        }
                    } else {
                        Part {
                            id,
                            range: Some(rg),
                            reverse,
                        }
                    });
                };
                let mut last = Range::EMPTY;
                let mut have_range = false;
                for &rg in ranges.iter() {
                    if merge_ranges(&mut last, rg, reverse) {
                        have_range = true;
                        continue;
                    }
                    if have_range {
                        push(last);
                    }
                    last = rg;
                    have_range = true;
                }
                if have_range {
                    push(last);
                }
            }
            written.push((index, pushed));
        }
        note_ties(&groups, &written, &mut out);
        for (_, pushed) in written {
            out.parts.extend(pushed);
        }
    }
    out
}

/// Sets `out.written` when the groups that the address orders write ranges on one strand
/// (`written`: the parts each group writes, in the order of `groups`) whose order matters: the
/// offsets of one are not those of another repeated (a list `x` and a list `y` give the same
/// offsets in both orders exactly when `x + y == y + x`).
fn note_ties(groups: &[Group], written: &[(usize, Vec<Part>)], out: &mut Loc) {
    let mut i = 0;
    while i < written.len() {
        let order = mask_handle_order(groups[written[i].0].packed);
        let mut j = i + 1;
        while j < written.len() && mask_handle_order(groups[written[j].0].packed) == order {
            j += 1;
        }
        let tied: Vec<&(usize, Vec<Part>)> = written[i..j]
            .iter()
            .filter(|(g, _)| groups[*g].handle != GroupHandle::Null)
            .collect();
        if tied.len() >= 2 {
            let names = || {
                tied.iter()
                    .filter_map(|(g, _)| groups[*g].id.as_ref())
                    .map(|id| id_label(&id.seq_id))
                    .collect::<Vec<_>>()
                    .join(" and ")
            };
            let shown: Vec<Vec<(u32, u32)>> = tied
                .iter()
                .map(|(_, parts)| {
                    parts
                        .iter()
                        .filter_map(|p| p.range.map(|r| (r.from, r.to_open)))
                        .collect::<Vec<_>>()
                })
                .filter(|offsets| !offsets.is_empty())
                .collect();
            let commute = shown.iter().enumerate().all(|(a, x)| {
                shown[a + 1..]
                    .iter()
                    .all(|y| [x.as_slice(), y].concat() == [y.as_slice(), x].concat())
            });
            if !commute && out.written.is_none() {
                out.written = Some(format!(
                    "the masks of the Seq-ids {} of one sequence are written in an order that NCBI \
                     takes from the memory addresses of their handles (equal packed numbers)",
                    names()
                ));
            }
        }
        i = j;
    }
}

/// The key under which `m_Values[algo]` finds a handle.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum HandleKey {
    Key(Vec<u8>),
    Giim(usize),
}

/// `CMaskedRangeSet`: the Seq-loc of each (algorithm, handle) and the
/// ranges of the algorithms that have any, in the order of their first
/// `insert`.
#[derive(Debug, Default)]
pub(crate) struct MaskedRangeSet {
    /// `m_Ranges`.
    ranges: Vec<MaskedRanges>,
    /// `m_Values`, by algorithm id: the Seq-loc of each handle that can be found again.
    values: BTreeMap<i32, HashMap<HandleKey, Loc>>,
    /// The `giim` handles of the masks (shared by the algorithms, as NCBI's tree is).
    giim: GiimTree,
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
    // The first Seq-loc of an (algorithm, handle) is kept as it was given: unsorted,
    // overlapping, with its ranges in the order `CSeq_loc_CI` visits them. Each later one is added
    // to the stored one with `CSeq_loc::Add`. `x_Set` makes the map element (and the algorithm's
    // entry of `m_Ranges`) before anything is combined.
    /// `Insert`: adds the Seq-loc (`parts`) of the Seq-id `id` for the algorithm `algo_id`, under
    /// the Seq-id's handle. Fails (LOSAT's refusal) only for a negative algorithm id, which NCBI's
    /// vector cannot take.
    pub(crate) fn insert(
        &mut self,
        algo_id: i32,
        id: ParsedSeqId,
        parts: Vec<MaskPart>,
    ) -> Result<(), WriteDbError> {
        if algo_id < 0 {
            return Err(WriteDbError::Unsupported(format!(
                "the masking algorithm id {algo_id} is negative"
            )));
        }
        if !self.ranges.iter().any(|r| r.algorithm_id == algo_id) {
            self.ranges.push(MaskedRanges {
                algorithm_id: algo_id,
                offsets: Vec::new(),
            });
        }
        let algorithm = self.values.entry(algo_id).or_default();
        let loc = Loc::given(id, parts);
        let key = match &loc.ids[0].handle {
            IdHandle::Key { key, .. } => HandleKey::Key(key.clone()),
            // A new handle that nothing finds again.
            IdHandle::Fresh => return Ok(()),
            IdHandle::Giim(giim) => HandleKey::Giim(self.giim.find_or_create(giim).0),
        };
        match algorithm.remove(&key) {
            None => {
                algorithm.insert(key, loc);
            }
            Some(stored) => {
                let combined = add_locs(&mut self.giim, &stored, &loc);
                algorithm.insert(key, combined);
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
    /// `GetRanges`: the ranges of the sequence that has the Seq-ids `ids` (`bs->GetId()`), looked
    /// up in that order: one `MaskedRanges` for each algorithm that has any mask (in the order
    /// of its first `insert`), with no offsets where the sequence has none. A range whose end is
    /// before its start or beyond the sequence is passed on as it is (the writer checks it).
    pub(crate) fn get_ranges(
        &mut self,
        ids: &[ParsedSeqId],
    ) -> Result<&[MaskedRanges], WriteDbError> {
        let Self {
            ranges,
            values,
            giim,
        } = self;
        // The handle of each Seq-id (`None`: one that no mask has).
        let keys: Vec<Option<HandleKey>> = ids
            .iter()
            .map(|id| match LocId::new(id.clone()).handle {
                IdHandle::Key { key, .. } => Some(HandleKey::Key(key)),
                IdHandle::Fresh => None,
                IdHandle::Giim(id) => giim.find(&id).map(HandleKey::Giim),
            })
            .collect();
        for algo in ranges.iter_mut() {
            algo.offsets.clear();
            let Some(locs) = values.get(&algo.algorithm_id) else {
                continue;
            };
            let mut found: Option<std::borrow::Cow<'_, Loc>> = None;
            for key in keys.iter().flatten() {
                let Some(stored) = locs.get(key) else {
                    continue;
                };
                found = Some(match found {
                    None => std::borrow::Cow::Borrowed(stored),
                    Some(previous) => std::borrow::Cow::Owned(add_locs(giim, &previous, stored)),
                });
            }
            if let Some(loc) = found {
                if let Some(why) = &loc.written {
                    return Err(WriteDbError::Unsupported(format!(
                        "{why}; LOSAT does not write this mask"
                    )));
                }
                algo.offsets.extend(loc.offsets());
            }
        }
        Ok(ranges)
    }
}

/// The Seq-id as `AsFastaString()` prints it.
fn id_label(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    crate::seqid::label::fasta(id, &mut out);
    String::from_utf8_lossy(&out).into_owned()
}

#[cfg(test)]
#[path = "range_set_tests.rs"]
mod tests;
