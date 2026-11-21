use super::*;
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, BlastCompoAlignmentContext::PreliminaryHspIndex, EMatrixAdjustRule,
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

fn fields(list: &Option<Box<BlastCompoAlignment>>) -> Vec<(i32, i32, i32, i32, i32, i32, String)> {
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
