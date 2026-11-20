//! No NCBI counterpart: `LOSAT_X_REDOCONTAINIDX`. The alignments one composition redo has kept for
//! a query (NCBI's list `alignments[query_index]` of `Blast_RedoOneMatch`, newest first) with an
//! index by subject start, so the two questions asked of the list for every HSP, `s_IsContained`
//! and `s_WithDistinctEnds`, visit only the alignments whose subject range can meet the new one.
//!
//! NCBI walks the whole list each time (`redo_alignment.c` 385-450, 965-996); a (query, subject)
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

/// `s_IsContained(in_align, alignments[query_index], lambda)`: from the query's index if it has
/// one, else from NCBI's list `list`.
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

/// `s_WithDistinctEnds(&newAlign, &alignments[query_index], ..., is_same_adjustment)`: into the
/// query's index if it has one, else into NCBI's list `list`, which moves into a new index (and
/// `list` becomes `None`) once it holds `index_from` alignments (`None`: never, the switch off).
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
// ```c
//     for (align = alignments;  align != NULL;  align = align->next ) {
//         /* for all elements of alignments */
//         if (KAPPA_SIGN(in_align->frame) == KAPPA_SIGN(align->frame)) {
//             /* hsp1 and hsp2 are in the same query/subject frame */
//             if (KAPPA_CONTAINED_IN_HSP
//                 (align->queryStart, align->queryEnd, query_offset,
//                  align->matchStart, align->matchEnd, subject_offset) &&
//                 KAPPA_CONTAINED_IN_HSP
//                 (align->queryStart, align->queryEnd, query_end,
//                  align->matchStart, align->matchEnd, subject_end) &&
//                 scoreThresh <= align->score) {
//                 return 1;
//             }
//         }
//     }
//     return 0;
// ```
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
// ```c
//     for (align = oldAlignments;  align != NULL;  align = align->next) {
//         ...
//         if (align->frame == newAlign->frame &&
//             ((is_same_adjustment && s_IsSameEndPoint(newAlign, align)) ||
//              (!is_same_adjustment && s_IsSimilarEndPoint(newAlign, align)))) {
//             ...
//             if (newAlign->score <= align->score) {
//                 /* newAlign cannot be added to the list. */
//                 include_new_align = 0;
//                 break;
//             }
//         }
//     }
//     if (include_new_align) {
//         ...
//         while (align != NULL) {
//             ...
//             if (align->frame == newAlign->frame &&
//                 ((align->queryStart == newAlign->queryStart &&
//                   align->matchStart == newAlign->matchStart)
//                  || (align->queryEnd == newAlign->queryEnd &&
//                      align->matchEnd == newAlign->matchEnd))) {
//                 /* The alignment shares an end with newAlign; */
//                 /* delete it. */
//                 BlastCompo_AlignmentsFree(&align, free_align_context);
//             } else { /* The alignment does not share an end with newAlign; */
//                 /* add it to the output list. */
//                 *tail =  align;
//                 tail  = &align->next;
//             }
//             align = align_next;
//         } /* end while align != NULL */
//         *p_oldAlignments = newAlign;
// ```
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
mod tests {
    use super::*;
    use crate::core::composition_adjustment::redo_alignment::{
        blast_compo_alignment_new, BlastCompoAlignmentContext::PreliminaryHspIndex,
        EMatrixAdjustRule,
    };

    /// A small linear congruential generator (the tests need no external crate).
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self, n: u32) -> i32 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 33) % u64::from(n)) as i32
        }
    }

    /// Small coordinate ranges so that shared ends, contained HSPs and equal scores are common;
    /// one alignment in 40 has its end before its start.
    fn random_align(rng: &mut Lcg, frames: &[i32], tag: usize) -> Box<BlastCompoAlignment> {
        let frame = frames[rng.next(frames.len() as u32) as usize];
        let qs = rng.next(40);
        let qe = qs + rng.next(20);
        let ms = rng.next(60) - 10;
        let me = if rng.next(40) == 0 {
            ms - 1 - rng.next(5)
        } else {
            ms + rng.next(25)
        };
        blast_compo_alignment_new(
            1 + rng.next(30),
            EMatrixAdjustRule::DontAdjustMatrix,
            qs,
            qe,
            0,
            ms,
            me,
            frame,
            Some(PreliminaryHspIndex(tag)),
        )
    }

    fn fields(
        list: &Option<Box<BlastCompoAlignment>>,
    ) -> Vec<(i32, i32, i32, i32, i32, i32, String)> {
        let mut out = Vec::new();
        let mut cur = list.as_deref();
        while let Some(a) = cur {
            out.push((
                a.score,
                a.query_start,
                a.query_end,
                a.match_start,
                a.match_end,
                a.frame,
                format!("{:?}", a.context),
            ));
            cur = a.next.as_deref();
        }
        out
    }

    /// A redo's per-query state as `blast_redo_one_match_...` holds it: NCBI's lists and the indexes.
    struct Redo {
        lists: Vec<Option<Box<BlastCompoAlignment>>>,
        indexes: XAlignIndexes,
        index_from: Option<usize>,
    }

    impl Redo {
        fn new(queries: usize, index_from: Option<usize>) -> Self {
            Self {
                lists: (0..queries).map(|_| None).collect(),
                indexes: Vec::new(),
                index_from,
            }
        }

        fn contains(&self, query: usize, align: &BlastCompoAlignment, lambda: f64) -> bool {
            x_contains(&self.indexes, query, &self.lists[query], align, lambda)
        }

        fn add(&mut self, query: usize, align: Box<BlastCompoAlignment>, same: bool) {
            let mut new_align = Some(align);
            x_add_distinct(
                &mut self.indexes,
                query,
                &mut self.lists[query],
                &mut new_align,
                same,
                self.index_from,
            );
            assert!(new_align.is_none());
        }

        fn indexed(&self, query: usize) -> bool {
            matches!(self.indexes.get(query), Some(Some(_)))
        }

        fn finish(self) -> (Vec<Option<Box<BlastCompoAlignment>>>, bool) {
            let mut lists = self.lists;
            let indexed = x_finish(self.indexes, &mut lists);
            (lists, indexed)
        }
    }

    /// One random sequence of calls on three queries, against NCBI's list functions: every
    /// `x_contains` answer and every final list (fields, context tags and order). `ends` puts the
    /// match coordinates near both ends of i32, with alignments that span from one end to the
    /// other (so `dmax` is near 2^32 and the ranges clamp). Returns the queries indexed.
    fn same_as_lists(
        trial: u64,
        frames: &[i32],
        index_from: Option<usize>,
        steps: usize,
        ends: bool,
    ) -> usize {
        const EDGE: i32 = i32::MAX - 1_000;
        let mut rng = Lcg(trial * 7 + 1);
        let mut reference: Vec<Option<Box<BlastCompoAlignment>>> = vec![None, None, None];
        let mut redo = Redo::new(3, index_from);
        let lambda = if rng.next(2) == 0 { 0.267 } else { 0.3176 };
        for step in 0..steps {
            let query = rng.next(3) as usize;
            let mut align = random_align(&mut rng, frames, step);
            if ends {
                let (ms, me) = match rng.next(6) {
                    0 => (-EDGE, EDGE),
                    1 => (EDGE, -EDGE),
                    2 | 3 => {
                        let ms = EDGE - rng.next(60);
                        (ms, (ms + rng.next(30)).min(EDGE))
                    }
                    _ => {
                        let ms = -EDGE + rng.next(60);
                        (ms, ms + rng.next(30))
                    }
                };
                align.match_start = ms;
                align.match_end = me;
            }
            if rng.next(2) == 0 {
                assert_eq!(
                    redo.contains(query, &align, lambda),
                    is_contained(&align, &reference[query], lambda),
                    "trial {trial} step {step}"
                );
            } else {
                let same = rng.next(2) == 0;
                with_distinct_ends(&mut Some(align.clone()), &mut reference[query], same);
                redo.add(query, align, same);
            }
        }
        let indexed = (0..3).filter(|&query| redo.indexed(query)).count();
        let (lists, any) = redo.finish();
        assert_eq!(any, indexed > 0, "trial {trial}");
        for query in 0..3 {
            assert_eq!(
                fields(&lists[query]),
                fields(&reference[query]),
                "trial {trial} query {query}"
            );
        }
        indexed
    }

    /// For blastp (frame 0), translated frames and two frames, and indexes started from the first
    /// alignment, from a few, from `X_INDEX_FROM` and never (`None`, the switch off).
    #[test]
    fn the_index_answers_as_ncbis_list() {
        let frame_sets: [&[i32]; 3] = [&[0], &[1, 2, 3, -1, -2, -3], &[1, -1]];
        let index_froms = [Some(0), Some(1), Some(4), Some(X_INDEX_FROM), None];
        let mut indexed = [0usize; 5];
        for trial in 0..600u64 {
            let which = trial as usize % 5;
            indexed[which] += same_as_lists(
                trial,
                frame_sets[trial as usize % 3],
                index_froms[which],
                1_200,
                false,
            );
        }
        // 120 runs (360 queries) per threshold; `None` never indexes.
        assert!(
            indexed[..4].iter().all(|&queries| queries >= 120) && indexed[4] == 0,
            "queries that reached the index per threshold: {indexed:?}"
        );
    }

    #[test]
    fn coordinates_at_the_ends_of_i32_give_the_same_answers() {
        for trial in 0..200u64 {
            let index_from = Some(if trial % 2 == 0 { 0 } else { 8 });
            same_as_lists(trial, &[1, -1], index_from, 600, true);
        }
    }

    #[test]
    fn an_empty_redo_contains_nothing_and_returns_no_list() {
        let probe = blast_compo_alignment_new(
            5,
            EMatrixAdjustRule::DontAdjustMatrix,
            0,
            10,
            0,
            0,
            10,
            1,
            None,
        );
        for index_from in [Some(0), None] {
            let redo = Redo::new(2, index_from);
            assert!(!redo.contains(1, &probe, 0.267));
            let (lists, indexed) = redo.finish();
            assert!(lists.iter().all(Option::is_none) && !indexed);
        }
    }

    /// A list moves into an index when it reaches the threshold, and only that query's; with the
    /// switch off (`None`) no index is made and the vector of indexes stays unallocated.
    #[test]
    fn a_list_moves_into_the_index_once_it_reaches_the_threshold() {
        let mut rng = Lcg(11);
        let mut redo = Redo::new(3, Some(4));
        let mut step = 0;
        while !redo.indexed(1) {
            redo.add(1, random_align(&mut rng, &[1], step), true);
            step += 1;
        }
        assert!(redo.lists[1].is_none());
        let index = redo.indexes[1].as_deref().expect("index of query 1");
        assert!(index.slots.iter().flatten().count() >= 4);
        assert!(!redo.indexed(0) && !redo.indexed(2));
        let mut off = Redo::new(3, None);
        for step in 0..200 {
            off.add(step % 3, random_align(&mut rng, &[1], step), true);
        }
        assert!(off.indexes.is_empty() && off.indexes.capacity() == 0);
    }
}
