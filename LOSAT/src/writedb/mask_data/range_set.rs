//! The masked ranges of the sequences, by algorithm and sequence
//! (`CMaskedRangeSet`, `c++/src/app/blastdb/masked_range_set.cpp`), with the
//! combination of several Seq-locs of one sequence (`CSeq_loc::Add` with
//! the flags `fMerge_All` and `fSort`, Seq_loc.cpp:4283-4693, 5060-5084).
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
//! Seq-loc that the writer's bounds check then refuses as NCBI's does. The order in which the sort
//! leaves different empty ranges, which the comparators do not order, can change the merge, so the
//! sort (`sort_as_observed`) makes the comparator calls and leaves the order that the C++ `sort` of
//! NCBI's build was observed to (`tests/unit/sort_observed/README.md`).
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
//! `CSeq_id_Handle::operator<`: `m_Packed` less 1 as unsigned, then the address of the handle's
//! info (seq_id_handle.hpp:201-207); the plus-strand groups are written before the minus-strand
//! ones (`x_RangesToSeq_loc` once per map, 4688-4691), and ranges of different handles are never
//! merged. Two different handles with the same `m_Packed` (`gi|5` and `gb|AY000005.1|`, 5; two
//! unpacked handles such as a local id and a `general` string tag, or two `giim` handles, 0) are
//! ordered by address, so LOSAT refuses a combination that writes ranges of two such handles on one
//! strand in an order that matters (the ranges written of one are not those of the other repeated,
//! `note_ties`). One strand each is decided, as plus before minus (review probe p4).
//!
//! The refusal waits for the end of the build. `get_ranges` writes LOSAT's order (tied groups in
//! the order of their first range) and names the Seq-ids of the first algorithm whose combination
//! is tied (`SequenceMasks::address_ordered`); `BuildDatabase` keeps the first such note and raises
//! it (`address_order_refusal`) only when the build would otherwise succeed, at its close. NCBI's
//! order is visible only in the column files that a successful build keeps, so every error of
//! NCBI that comes later in the run decides the outcome first: the bounds check of the same
//! sequence over all its algorithms (`SetMaskData`, writedb_impl.cpp:1339-1386), and any error of
//! a later sequence or file, after which `EndBuild` with `erase` true removes the files
//! (build_db.cpp:1415-1417). A refusal at the close leaves the files that such an error leaves
//! (`BuildDatabase::end_build`).
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
    mask_handle_key, mask_handle_order, mask_handle_packed, GiimTree, TreeHandle,
};
use crate::seqid::model::GiimportId;
use crate::seqid::ParsedSeqId;
use crate::writedb::mask_info::MaskPart;
use crate::writedb::WriteDbError;

// NCBI reference (598d8ae6): c++/include/objects/blastdb/defline_extra.hpp:78-93
/// `SBlastDbMaskData`: the masked ranges of one sequence for one algorithm,
/// each `(start, end)` with `end` exclusive.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct MaskedRanges {
    pub(crate) algorithm_id: i32,
    pub(crate) offsets: Vec<(u32, u32)>,
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb.hpp:64-79
/// `CMaskedRangesVector::empty`: no algorithm has a range.
pub(crate) fn ranges_are_empty(ranges: &[MaskedRanges]) -> bool {
    ranges.iter().all(|r| r.offsets.is_empty())
}

// NCBI reference (598d8ae6): c++/include/util/range.hpp:148-151,178-191,256-308,312-353,358-367
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
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4301-4322
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
    // Abutting: `CopyTo` sets `to_open` to one more than `rg2.GetTo()`, `CopyFrom` sets `from`.
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

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4565-4572
// NCBI sorts the ranges of a group with the C++ library's `sort`, whose order of the elements that
// the comparator leaves unordered is not specified. The functions below were written from what the
// `sort` of g++ 13.3 was seen to do with a comparator that recorded each call (the two elements it
// was handed, in argument order, and its answer) and with the order it left: every answer sequence
// of a strict weak order for 0 to 9 elements, random strict weak orders with many ties, NCBI's range
// comparators, and comparators that decide lazily to force long chains of splits
// (`tests/unit/sort_observed/`). They make the same calls, in the same order, and leave the same
// order; no library source was read.

/// Sorts the ranges of one group and strand with `range_less` (`CRange_ReverseLess` when
/// `reverse`). Each strand has its own comparator with nothing captured, so the sort neither
/// tests the strand nor reads it from memory at each comparison.
fn sort_ranges(ranges: &mut [Range], reverse: bool) {
    if reverse {
        sort_as_observed(ranges, &mut |a: &Range, b: &Range| range_less(a, b, true));
    } else {
        sort_as_observed(ranges, &mut |a: &Range, b: &Range| range_less(a, b, false));
    }
}

/// Stretches of at most this many elements are not split; the closing insertion pass orders them.
const SPLIT_ABOVE: usize = 16;

/// Sorts `v` by `less` as observed: stretches longer than [`SPLIT_ABOVE`] are split around a pivot
/// (`split_stretch`) at most `2 * floor(log2(len))` times along any chain of splits, a stretch that
/// would need more is ordered through a heap (`heap_order`), and an insertion pass over the whole
/// slice finishes (`insertion_pass`). Fewer than two elements: no call.
fn sort_as_observed<T: Copy>(v: &mut [T], less: &mut impl FnMut(&T, &T) -> bool) {
    let len = v.len();
    if len < 2 {
        return;
    }
    split_stretch(v, 2 * len.ilog2() as usize, less);
    insertion_pass(v, less);
}

/// While `v` is longer than [`SPLIT_ABOVE`]: with no split left, orders it through a heap; else
/// moves a pivot to its front, splits it, handles the part right of the cut first (with the splits
/// left at this level) and then goes on with the part left of it.
fn split_stretch<T: Copy>(v: &mut [T], mut splits: usize, less: &mut impl FnMut(&T, &T) -> bool) {
    let mut end = v.len();
    while end > SPLIT_ABOVE {
        if splits == 0 {
            heap_order(&mut v[..end], less);
            return;
        }
        splits -= 1;
        pivot_to_front(&mut v[..end], less);
        let cut = split_around_front(&mut v[..end], less);
        split_stretch(&mut v[cut..end], splits, less);
        end = cut;
    }
}

/// Swaps into `v[0]` the middle one of `v[1]`, `v[len / 2]` and `v[len - 1]`. The calls: first
/// `(v[1], v[len / 2])`; if true, `(v[len / 2], v[len - 1])` (true: `v[len / 2]`), then
/// `(v[1], v[len - 1])` (true: `v[len - 1]`, false: `v[1]`); if false, `(v[1], v[len - 1])` (true:
/// `v[1]`), then `(v[len / 2], v[len - 1])` (true: `v[len - 1]`, false: `v[len / 2]`).
fn pivot_to_front<T: Copy>(v: &mut [T], less: &mut impl FnMut(&T, &T) -> bool) {
    let (a, b, c) = (1, v.len() / 2, v.len() - 1);
    let chosen = if less(&v[a], &v[b]) {
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
    v.swap(0, chosen);
}

/// Splits `v[1..]` around the pivot `v[0]`, which stays in place, and returns the cut. A left
/// mark from 1 moves right while `(element, pivot)` is true; a right mark from the end moves left,
/// one step first and then while `(pivot, element)` is true. When the left mark is not before the
/// right one, the left mark is the cut; else the two elements are swapped, the left mark steps
/// right and the scan goes on. A strict weak order always stops both marks inside `v`.
fn split_around_front<T: Copy>(v: &mut [T], less: &mut impl FnMut(&T, &T) -> bool) -> usize {
    let mut left = 1;
    let mut right = v.len();
    loop {
        while less(&v[left], &v[0]) {
            left += 1;
        }
        right -= 1;
        while less(&v[0], &v[right]) {
            right -= 1;
        }
        if left >= right {
            return left;
        }
        v.swap(left, right);
        left += 1;
    }
}

/// Orders `v` through a heap with its greatest element at the front, the way `sort` did on a
/// stretch it could not split any more and `partial_sort` does over a whole range: the heap is
/// built by sinking the elements from position `(len - 2) / 2` down to 0, then the front is swapped
/// to the end of the heap, which shrinks by one, until one element is left.
fn heap_order<T: Copy>(v: &mut [T], less: &mut impl FnMut(&T, &T) -> bool) {
    let len = v.len();
    if len < 2 {
        return;
    }
    let mut top = (len - 2) / 2;
    loop {
        let value = v[top];
        sink(v, top, value, less);
        if top == 0 {
            break;
        }
        top -= 1;
    }
    for end in (1..len).rev() {
        let value = v[end];
        v[end] = v[0];
        sink(&mut v[..end], 0, value, less);
    }
}

/// Puts `value` into the heap `heap` whose position `top` is free: the free position moves down,
/// each step to the greater child (the right one unless `(right, left)` is true), until it has no
/// two children; at an even length a free position whose only child is the last element takes that
/// child. Then `value` climbs from there towards `top` while `(parent, value)` is true.
fn sink<T: Copy>(heap: &mut [T], top: usize, value: T, less: &mut impl FnMut(&T, &T) -> bool) {
    let len = heap.len();
    let mut free = top;
    let mut child = top;
    while child < (len - 1) / 2 {
        child = 2 * (child + 1);
        if less(&heap[child], &heap[child - 1]) {
            child -= 1;
        }
        heap[free] = heap[child];
        free = child;
    }
    if len.is_multiple_of(2) && child == (len - 2) / 2 {
        child = 2 * (child + 1);
        heap[free] = heap[child - 1];
        free = child - 1;
    }
    while free > top {
        let parent = (free - 1) / 2;
        if !less(&heap[parent], &value) {
            break;
        }
        heap[free] = heap[parent];
        free = parent;
    }
    heap[free] = value;
}

/// Inserts each element, from the second on, into the ordered elements before it. Among the
/// first [`SPLIT_ABOVE`] positions an element is first compared with the front, `(element,
/// v[0])`, and goes to the front if that is true; otherwise (and always from position
/// [`SPLIT_ABOVE`] on) it moves left past its neighbours while `(element, neighbour)` is true.
/// After the splits a least element is among the first [`SPLIT_ABOVE`] positions, so under a
/// strict weak order the walk left stops inside `v`.
fn insertion_pass<T: Copy>(v: &mut [T], less: &mut impl FnMut(&T, &T) -> bool) {
    let head = v.len().min(SPLIT_ABOVE);
    for i in 1..v.len() {
        let value = v[i];
        if i < head && less(&value, &v[0]) {
            v.copy_within(0..i, 1);
            v[0] = value;
            continue;
        }
        let mut place = i;
        while less(&value, &v[place - 1]) {
            v[place] = v[place - 1];
            place -= 1;
        }
        v[place] = value;
    }
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
    /// `IsReverse` of the strand: the strand is used for nothing else.
    reverse: bool,
}

/// A Seq-loc as `CMaskedRangeSet` holds it: the one given, or the result of `CSeq_loc::Add`.
#[derive(Clone, Debug, Default)]
struct Loc {
    ids: Vec<LocId>,
    parts: Vec<Part>,
    /// Set when the ranges NCBI writes of this Seq-loc depend on memory addresses: the Seq-ids
    /// of the tied handles (`AsFastaString`, joined by " and ").
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
// NCBI reference (598d8ae6): c++/include/objects/seqloc/Na_strand.hpp:74-79
const NA_STRAND_MINUS: u8 = 2;
const NA_STRAND_BOTH_REV: u8 = 4;

fn is_reverse(strand: Option<u8>) -> bool {
    matches!(strand, Some(NA_STRAND_MINUS | NA_STRAND_BOTH_REV))
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
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4640-4693
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4546-4594
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1439-1471,1500-1531
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:4392-4437
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:335-338
// `CSeq_loc_CI` (`eEmpty_Allow`) first makes the handle of every range's Seq-id (`x_SetId`), then
// the loop takes `GetBestSynonym` = `GetHandle` of each again; both sets of handles go when
// `x_MergeAndSort` returns. A `null` Seq-loc has an unset Seq-id, whose handle is the null one:
// its group writes one `null` Seq-loc. A whole range becomes a `whole` Seq-loc and an empty one an
// `empty` Seq-loc, both without strand (the next iteration sees them on the plus strand).
/// `Add` of `L2` to `L1` with the flags `fMerge_All` and `fSort`, and no synonym mapper: the
/// Seq-loc `x_MergeAndSort` writes, with what in it depends on memory addresses (`Loc::written`).
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
    // `operator<`: `m_Packed` less 1, unsigned, then the info's address (the null handle's is 0).
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
                sort_ranges(ranges, reverse);
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
                out.written = Some(names());
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
    // NCBI reference (598d8ae6): c++/src/app/blastdb/masked_range_set.cpp:86-106
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
    /// `GetRanges`: the ranges of the sequence that has the Seq-ids `ids` (`bs->GetId()`), looked
    /// up in that order: one `MaskedRanges` for each algorithm that has any mask (in the order
    /// of its first `insert`), with no offsets where the sequence has none. A range whose end is
    /// before its start or beyond the sequence is passed on as it is (the writer checks it).
    /// Every algorithm is computed; where NCBI's order is the addresses', LOSAT's is written and
    /// the first such algorithm is named (`SequenceMasks::address_ordered`).
    pub(crate) fn get_ranges(&mut self, ids: &[ParsedSeqId]) -> SequenceMasks<'_> {
        let Self {
            ranges,
            values,
            giim,
        } = self;
        let mut address_ordered = None;
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
                if address_ordered.is_none() {
                    address_ordered.clone_from(&loc.written);
                }
                algo.offsets.extend(loc.offsets());
            }
        }
        SequenceMasks {
            ranges,
            address_ordered,
        }
    }
}

/// What `MaskedRangeSet::get_ranges` gives for one sequence.
#[derive(Debug)]
pub(crate) struct SequenceMasks<'a> {
    /// `GetRanges`' vector: the offsets of each algorithm that has any mask.
    pub(crate) ranges: &'a [MaskedRanges],
    /// The Seq-ids (`AsFastaString`, joined by " and ") of the first algorithm whose ranges NCBI
    /// writes in the order of the memory addresses of their handles; `ranges` has LOSAT's order
    /// there. `None` when NCBI's order is decided for every algorithm.
    pub(crate) address_ordered: Option<String>,
}

/// LOSAT's refusal of a build in which the masks of the Seq-ids `seq_ids` of the sequence
/// `sequence` (its first Seq-id) are written in an order that NCBI takes from memory addresses
/// (`SequenceMasks::address_ordered`): "losat makeblastdb: ...", exit status 1.
pub(crate) fn address_order_refusal(seq_ids: &str, sequence: &str) -> WriteDbError {
    WriteDbError::Unsupported(format!(
        "the masks of the Seq-ids {seq_ids} of the sequence {sequence} are written in an order \
         that NCBI takes from the memory addresses of their handles (equal packed numbers); LOSAT \
         does not write this database"
    ))
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
