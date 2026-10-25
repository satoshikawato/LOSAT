//! No NCBI counterpart: `LOSAT_X_REDOCONTAINIDX`. The alignments one composition redo has kept for
//! a query (NCBI's list `alignments[query_index]` of `Blast_RedoOneMatch`, newest first) with an
//! index by subject start, so the two questions asked of the list for every HSP, `s_IsContained`
//! and `s_WithDistinctEnds`, visit only the alignments whose subject range can meet the new one.
//!
//! NCBI walks the whole list each time (`redo_alignment.c` 386-451, 965-997); a (query, subject)
//! match with tens of thousands of HSPs (a protein against a human chromosome) spends most of its
//! redo there. The answers are the same: `contains` is true when some alignment contains the HSP,
//! `add_distinct` keeps the new alignment unless some alignment with a shared end scores at least as
//! well and then removes those with an equal start or end, and `into_list` returns the kept
//! alignments newest first, which is the order NCBI's list has (the new alignment goes to the head;
//! the others keep their order).
//!
//! While a list is short (fewer than `X_INDEX_FROM` alignments) it stays NCBI's list and the
//! original functions answer; past that it moves into the index once.
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

/// The alignments kept for one query.
pub(crate) struct XAlignSet {
    /// NCBI's list while the set is not indexed.
    list: Option<Box<BlastCompoAlignment>>,
    /// The list's length while not indexed (counted only when `index_from` is set).
    len: usize,
    /// `None` keeps NCBI's list for good (the switch is off).
    index_from: Option<usize>,
    indexed: bool,
    /// Indexed: the alignments in the order they were kept; `None` where one was removed.
    slots: Vec<Option<Box<BlastCompoAlignment>>>,
    /// Indexed: (frame, match start, slot) of each alignment in `slots`.
    by_start: BTreeSet<(i32, i32, u32)>,
    /// Indexed: the frames seen.
    frames: Vec<i32>,
    /// Indexed: the largest and smallest `match_end - match_start` of any alignment added.
    dmax: i64,
    dmin: i64,
}

fn clamp(v: i64) -> i32 {
    v.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

impl XAlignSet {
    pub(crate) fn new(index_from: Option<usize>) -> Self {
        Self {
            list: None,
            len: 0,
            index_from,
            indexed: false,
            slots: Vec::new(),
            by_start: BTreeSet::new(),
            frames: Vec::new(),
            dmax: 0,
            dmin: 0,
        }
    }

    /// Whether the set moved into the index.
    pub(crate) fn indexed(&self) -> bool {
        self.indexed
    }

    /// `s_IsContained(in_align, alignments[query_index], lambda)`.
    pub(crate) fn contains(&self, in_align: &BlastCompoAlignment, lambda: f64) -> bool {
        if !self.indexed {
            return is_contained(in_align, &self.list, lambda);
        }
        x_contains_indexed(self, in_align, lambda)
    }

    /// `s_WithDistinctEnds(&newAlign, &alignments[query_index], ..., is_same_adjustment)`.
    pub(crate) fn add_distinct(
        &mut self,
        p_new_align: &mut Option<Box<BlastCompoAlignment>>,
        is_same_adjustment: bool,
    ) {
        if !self.indexed {
            with_distinct_ends(p_new_align, &mut self.list, is_same_adjustment);
            if let Some(from) = self.index_from {
                self.len = list_len(&self.list);
                if self.len >= from {
                    self.build_index();
                }
            }
            return;
        }
        let Some(new_align) = p_new_align.take() else {
            return;
        };
        x_add_distinct_indexed(self, new_align, is_same_adjustment);
    }

    /// The kept alignments, newest first (NCBI's list).
    pub(crate) fn into_list(self) -> Option<Box<BlastCompoAlignment>> {
        if !self.indexed {
            return self.list;
        }
        let mut head: Option<Box<BlastCompoAlignment>> = None;
        for mut align in self.slots.into_iter().flatten() {
            align.next = head;
            head = Some(align);
        }
        head
    }

    /// Moves NCBI's list (newest first) into the slots (oldest first) and the index.
    fn build_index(&mut self) {
        let mut newest_first = Vec::with_capacity(self.len);
        let mut cur = self.list.take();
        while let Some(mut align) = cur {
            cur = align.next.take();
            newest_first.push(align);
        }
        self.indexed = true;
        self.slots.reserve(newest_first.len() * 2);
        for align in newest_first.into_iter().rev() {
            self.insert(align);
        }
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
}

fn list_len(list: &Option<Box<BlastCompoAlignment>>) -> usize {
    let mut n = 0;
    let mut cur = list.as_deref();
    while let Some(align) = cur {
        n += 1;
        cur = align.next.as_deref();
    }
    n
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
fn x_contains_indexed(set: &XAlignSet, in_align: &BlastCompoAlignment, lambda: f64) -> bool {
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
    set: &mut XAlignSet,
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

    /// Every `contains` answer and the final list (fields, context tags and order) equal those of
    /// NCBI's list functions on the same calls, for blastp (frame 0), translated frames, and
    /// indexes built from the first alignment, from a few, and from `X_INDEX_FROM`.
    #[test]
    fn the_index_answers_as_ncbis_list() {
        let frame_sets: [&[i32]; 3] = [&[0], &[1, 2, 3, -1, -2, -3], &[1, -1]];
        let index_froms = [Some(0), Some(1), Some(4), Some(X_INDEX_FROM), None];
        let mut indexed_runs = [0usize; 5];
        for trial in 0..600u64 {
            let mut rng = Lcg(trial * 7 + 1);
            let frames = frame_sets[trial as usize % frame_sets.len()];
            let which = trial as usize % 5;
            let index_from = index_froms[which];
            let mut reference: Option<Box<BlastCompoAlignment>> = None;
            let mut set = XAlignSet::new(index_from);
            let lambda = 0.267;
            for step in 0..400 {
                if rng.next(2) == 0 {
                    let probe = random_align(&mut rng, frames, usize::MAX);
                    assert_eq!(
                        set.contains(&probe, lambda),
                        is_contained(&probe, &reference, lambda),
                        "trial {trial} step {step}"
                    );
                } else {
                    let same = rng.next(2) == 0;
                    let new_align = random_align(&mut rng, frames, step);
                    let mut a = Some(new_align.clone());
                    let mut b = Some(new_align);
                    with_distinct_ends(&mut a, &mut reference, same);
                    set.add_distinct(&mut b, same);
                    assert!(a.is_none() && b.is_none());
                }
            }
            indexed_runs[which] += usize::from(set.indexed());
            assert_eq!(
                fields(&set.into_list()),
                fields(&reference),
                "trial {trial}"
            );
        }
        // 120 runs per threshold; `None` never indexes.
        assert!(
            indexed_runs[..4].iter().all(|&runs| runs >= 60) && indexed_runs[4] == 0,
            "runs that reached the index per threshold: {indexed_runs:?}"
        );
    }

    #[test]
    fn an_empty_set_contains_nothing_and_returns_no_list() {
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
            let set = XAlignSet::new(index_from);
            assert!(!set.contains(&probe, 0.267));
            assert!(set.into_list().is_none());
        }
    }

    /// Match coordinates near both ends of i32 (and alignments that span from one end to the
    /// other, so `dmax` is near 2^32 and the ranges clamp), against NCBI's list functions.
    #[test]
    fn coordinates_at_the_ends_of_i32_give_the_same_answers() {
        const EDGE: i32 = i32::MAX - 1_000;
        for trial in 0..200u64 {
            let mut rng = Lcg(trial * 13 + 5);
            let mut reference: Option<Box<BlastCompoAlignment>> = None;
            let mut set = XAlignSet::new(Some(if trial % 2 == 0 { 0 } else { 8 }));
            for step in 0..200 {
                let mut align = random_align(&mut rng, &[1, -1], step);
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
                if rng.next(2) == 0 {
                    assert_eq!(
                        set.contains(&align, 0.267),
                        is_contained(&align, &reference, 0.267),
                        "trial {trial} step {step}"
                    );
                } else {
                    let same = rng.next(2) == 0;
                    let mut a = Some(align.clone());
                    let mut b = Some(align);
                    with_distinct_ends(&mut a, &mut reference, same);
                    set.add_distinct(&mut b, same);
                }
            }
            assert_eq!(
                fields(&set.into_list()),
                fields(&reference),
                "trial {trial}"
            );
        }
    }

    #[test]
    fn a_list_moves_into_the_index_once_it_reaches_the_threshold() {
        let mut rng = Lcg(11);
        let mut set = XAlignSet::new(Some(4));
        let mut step = 0;
        while !set.indexed() {
            set.add_distinct(&mut Some(random_align(&mut rng, &[1], step)), true);
            step += 1;
        }
        assert!(set.list.is_none());
        assert!(set.slots.iter().flatten().count() >= 4);
        let mut off = XAlignSet::new(None);
        for step in 0..200 {
            off.add_distinct(&mut Some(random_align(&mut rng, &[1], step)), true);
        }
        assert!(!off.indexed());
    }
}
