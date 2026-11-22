//! No NCBI counterpart: `LOSAT_X_REDOCONTAINIDX`. The alignments one composition redo has kept for
//! a query (NCBI's list `alignments[query_index]` of `Blast_RedoOneMatch`, newest first) with an
//! index by subject start, so the two questions asked of the list for every HSP, `s_IsContained`
//! and `s_WithDistinctEnds`, visit only the alignments whose subject range can meet the new one.
//!
//! NCBI walks the whole list each time (`c++/src/algo/blast/composition_adjustment/redo_alignment.c` 385-450, 965-996); a (query, subject)
//! match with tens of thousands of HSPs (a protein against a human chromosome) spends most of its
//! redo there. The answers are the same: `x_contains` is true when some alignment contains the HSP,
//! `x_add_distinct` keeps the new alignment unless some alignment with a shared end scores at least
//! as well and then removes those with an equal start or end, and `x_finish` puts back the kept
//! alignments newest first, which is the order NCBI's list has (the new alignment goes to the head;
//! the others keep their order).
//!
//! The redo keeps NCBI's lists (`Vec<Option<Box<BlastCompoAlignment>>>`, one per query) and the
//! original functions answer while a list is shorter than `X_INDEX_FROM` alignments; a list that
//! reaches it moves, once, into an `XAlignIndex` of its own. The vector of indexes stays empty
//! (no allocation) in a redo where no list gets that long, which is most of them.
use std::collections::BTreeSet;
use std::sync::OnceLock;

use super::redo_alignment::{
    is_contained, is_similar_endpoint, with_distinct_ends, BlastCompoAlignment, KAPPA_BIT_TOL,
    LOCAL_LN2,
};

/// No NCBI counterpart: the switch, read once per process.
pub(crate) fn x_redo_contain_idx() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_REDOCONTAINIDX"))
}

/// The list length from which a query's alignments move into the index.
pub(crate) const X_INDEX_FROM: usize = 32;

/// The indexes of one redo, by query; empty or `None` where the query's list is NCBI's.
pub(crate) type XAlignIndexes = Vec<Option<Box<XAlignIndex>>>;

/// The alignments kept for one query once its list reached `X_INDEX_FROM`.
pub(crate) struct XAlignIndex {
    /// The alignments in the order they were kept; `None` where one was removed.
    slots: Vec<Option<Box<BlastCompoAlignment>>>,
    /// (frame, match start, slot) of each alignment in `slots`.
    by_start: BTreeSet<(i32, i32, u32)>,
    /// The frames seen.
    frames: Vec<i32>,
    /// The largest and smallest `match_end - match_start` of any alignment added.
    dmax: i64,
    dmin: i64,
}

fn clamp(v: i64) -> i32 {
    v.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

impl XAlignIndex {
    /// NCBI's list (newest first) moved into the slots (oldest first) and the index.
    fn from_list(list: Option<Box<BlastCompoAlignment>>) -> Box<Self> {
        let mut newest_first = Vec::new();
        let mut cur = list;
        while let Some(mut align) = cur {
            cur = align.next.take();
            newest_first.push(align);
        }
        let mut index = Box::new(Self {
            slots: Vec::with_capacity(newest_first.len() * 2),
            by_start: BTreeSet::new(),
            frames: Vec::new(),
            dmax: 0,
            dmin: 0,
        });
        for align in newest_first.into_iter().rev() {
            index.insert(align);
        }
        index
    }

    fn insert(&mut self, align: Box<BlastCompoAlignment>) {
        let slot = u32::try_from(self.slots.len()).expect("fewer than 2^32 alignments per query");
        let d = i64::from(align.match_end) - i64::from(align.match_start);
        if self.slots.is_empty() {
            self.dmax = d;
            self.dmin = d;
        } else {
            self.dmax = self.dmax.max(d);
            self.dmin = self.dmin.min(d);
        }
        if !self.frames.contains(&align.frame) {
            self.frames.push(align.frame);
        }
        self.by_start.insert((align.frame, align.match_start, slot));
        self.slots.push(Some(align));
    }

    /// The kept alignments, newest first (NCBI's list).
    fn into_list(self) -> Option<Box<BlastCompoAlignment>> {
        let mut head: Option<Box<BlastCompoAlignment>> = None;
        for mut align in self.slots.into_iter().flatten() {
            align.next = head;
            head = Some(align);
        }
        head
    }
}

/// `s_IsContained` (with `in_align`, `alignments[query_index]` and `lambda`): from the query's
/// index if it has one, else from NCBI's list `list`.
#[inline]
pub(crate) fn x_contains(
    indexes: &XAlignIndexes,
    query_index: usize,
    list: &Option<Box<BlastCompoAlignment>>,
    in_align: &BlastCompoAlignment,
    lambda: f64,
) -> bool {
    match indexes.get(query_index) {
        Some(Some(index)) => x_contains_indexed(index, in_align, lambda),
        _ => is_contained(in_align, list, lambda),
    }
}

/// `s_WithDistinctEnds` (with the new alignment `newAlign`, `alignments[query_index]` and
/// `is_same_adjustment`): into the query's index if it has one, else into NCBI's list `list`, which
/// moves into a new index (and `list` becomes `None`) once it holds `index_from` alignments
/// (`None`: never, the switch off).
#[inline]
pub(crate) fn x_add_distinct(
    indexes: &mut XAlignIndexes,
    query_index: usize,
    list: &mut Option<Box<BlastCompoAlignment>>,
    p_new_align: &mut Option<Box<BlastCompoAlignment>>,
    is_same_adjustment: bool,
    index_from: Option<usize>,
) {
    if let Some(Some(index)) = indexes.get_mut(query_index) {
        if let Some(new_align) = p_new_align.take() {
            x_add_distinct_indexed(index, new_align, is_same_adjustment);
        }
        return;
    }
    let head = list
        .as_deref()
        .map(|align| align as *const BlastCompoAlignment);
    with_distinct_ends(p_new_align, list, is_same_adjustment);
    let Some(from) = index_from else {
        return;
    };
    // The list grew only if the new alignment became its head.
    if list
        .as_deref()
        .map(|align| align as *const BlastCompoAlignment)
        != head
        && list_reaches(list, from)
    {
        x_start_index(indexes, query_index, list);
    }
}

#[inline(never)]
fn x_start_index(
    indexes: &mut XAlignIndexes,
    query_index: usize,
    list: &mut Option<Box<BlastCompoAlignment>>,
) {
    if indexes.len() <= query_index {
        indexes.resize_with(query_index + 1, || None);
    }
    indexes[query_index] = Some(XAlignIndex::from_list(list.take()));
}

/// Puts each index's alignments back into its query's list, newest first; returns whether any
/// query had an index.
pub(crate) fn x_finish(
    indexes: XAlignIndexes,
    lists: &mut [Option<Box<BlastCompoAlignment>>],
) -> bool {
    let indexed = !indexes.is_empty();
    for (query_index, index) in indexes.into_iter().enumerate() {
        if let Some(index) = index {
            lists[query_index] = index.into_list();
        }
    }
    indexed
}

/// Whether NCBI's list holds at least `n` alignments (visits at most `n` of them).
fn list_reaches(list: &Option<Box<BlastCompoAlignment>>, n: usize) -> bool {
    let mut seen = 0;
    let mut cur = list.as_deref();
    while let Some(align) = cur {
        seen += 1;
        if seen >= n {
            return true;
        }
        cur = align.next.as_deref();
    }
    seen >= n
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:980-995
// The same test, over the alignments of a frame of the same sign whose match start lies in
// [subject_end - dmax, subject_offset] only: an `align` that passes has matchStart <= subject_offset
// and matchEnd >= subject_end, so matchStart >= subject_end - (matchEnd - matchStart) >=
// subject_end - dmax. The answer does not depend on the order the alignments are met in.
fn x_contains_indexed(set: &XAlignIndex, in_align: &BlastCompoAlignment, lambda: f64) -> bool {
    let score_thresh = in_align.score as f64 + KAPPA_BIT_TOL * LOCAL_LN2 / lambda;
    let lo = clamp(i64::from(in_align.match_end) - set.dmax);
    let hi = in_align.match_start;
    if lo > hi {
        return false;
    }
    for &frame in &set.frames {
        if frame.signum() != in_align.frame.signum() {
            continue;
        }
        // Nearest start first: a containing alignment starts at or before the HSP, and the
        // answer does not depend on the order.
        for &(_, _, slot) in set
            .by_start
            .range((frame, lo, 0)..=(frame, hi, u32::MAX))
            .rev()
        {
            let align = set.slots[slot as usize]
                .as_deref()
                .expect("indexed slot holds an alignment");
            if align.query_start <= in_align.query_start
                && align.query_end >= in_align.query_start
                && align.match_start <= in_align.match_start
                && align.match_end >= in_align.match_start
                && align.query_start <= in_align.query_end
                && align.query_end >= in_align.query_end
                && align.match_start <= in_align.match_end
                && align.match_end >= in_align.match_end
                && score_thresh <= align.score as f64
            {
                return true;
            }
        }
    }
    false
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:401-446
// Both loops, over the alignments of the same frame whose match start lies in
// [min(ms, ms - dmax, me - dmax), max(ms, me, me - dmin)] only (ms, me: the new alignment's match
// start and end): a shared start has matchStart == ms; a shared end has matchStart == me -
// (matchEnd - matchStart), in [me - dmax, me - dmin]; a start contained in `align` (similar end
// point, `s_IsSimilarEndPoint`, 321-349) has matchStart in [ms - dmax, ms], a contained end in [me - dmax, me]. The first loop's
// answer does not depend on the order; the second removes the same alignments and the slots keep
// the others' order.
fn x_add_distinct_indexed(
    set: &mut XAlignIndex,
    new_align: Box<BlastCompoAlignment>,
    is_same_adjustment: bool,
) {
    let ms = i64::from(new_align.match_start);
    let me = i64::from(new_align.match_end);
    let lo = clamp(ms.min(ms - set.dmax).min(me - set.dmax));
    let hi = clamp(ms.max(me).max(me - set.dmin));
    let mut remove: Vec<u32> = Vec::new();
    for &(_, _, slot) in set
        .by_start
        .range((new_align.frame, lo, 0)..=(new_align.frame, hi, u32::MAX))
    {
        let align = set.slots[slot as usize]
            .as_deref()
            .expect("indexed slot holds an alignment");
        let same_start = align.query_start == new_align.query_start
            && align.match_start == new_align.match_start;
        let same_end =
            align.query_end == new_align.query_end && align.match_end == new_align.match_end;
        let shared = if is_same_adjustment {
            same_start || same_end
        } else {
            is_similar_endpoint(&new_align, align)
        };
        if shared && new_align.score <= align.score {
            return;
        }
        if same_start || same_end {
            remove.push(slot);
        }
    }
    for slot in remove {
        let align = set.slots[slot as usize]
            .take()
            .expect("indexed slot holds an alignment");
        set.by_start.remove(&(align.frame, align.match_start, slot));
    }
    set.insert(new_align);
}

#[cfg(test)]
#[path = "x_align_index_tests.rs"]
mod tests;
