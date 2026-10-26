//! No NCBI counterpart: `LOSAT_X_TBNREAPIDX`. The reap of the HSPs that one composition redo
//! gives a match (`s_HitlistReapContained`, `blast_kappa.c` 224-273; `reap_contained_postredo_hsps`
//! in `kappa.rs`) answered from an index of the HSPs kept so far by subject start.
//!
//! NCBI tests every HSP against every earlier HSP still kept, so a match with tens of thousands of
//! HSPs (a protein against a human chromosome) costs the square of that: the B11h probe (35 human
//! proteins x the human genome `-db`, 8 threads) counted 21.11 G iterations in 73.6 thread-s, 10 %
//! of the search's user time, and the index of this file visited 6,326 candidates in 0.30
//! thread-s with the same answers. The answer for an HSP is whether some earlier HSP still kept,
//! of the same context and frame, contains both of its end points (`CONTAINED_IN_HSP`,
//! `blast_hits_priv.h` 65-70) and scores at least as much; it does not depend on the order the
//! earlier HSPs are met in, so the index may meet them in any order. Lists shorter than
//! `X_REAP_INDEX_FROM` HSPs keep NCBI's loop.
use std::collections::BTreeSet;
use std::sync::OnceLock;

use super::GappedHsp;

/// The list length from which the reap takes the index.
pub(super) const X_REAP_INDEX_FROM: usize = 32;

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// Whether the reap takes the index (`LOSAT_X_TBNREAPIDX`, read once; a unit test sets it per test
/// thread).
fn x_reap_idx_on() -> bool {
    #[cfg(test)]
    if let Some(on) = tests::X_ON.get() {
        return on;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNREAPIDX"))
}

/// Whether the reap of a list of `len` HSPs takes the index (the switch on and `len` at least
/// `X_REAP_INDEX_FROM`); records the side with `LOSAT_X_REACH`.
pub(super) fn x_indexed(len: usize) -> bool {
    if !x_reap_idx_on() {
        return false;
    }
    let indexed = len >= X_REAP_INDEX_FROM;
    crate::utils::xreach::x_reach!("LOSAT_X_TBNREAPIDX", indexed);
    indexed
}

/// `CONTAINED_IN_HSP` of both end points of `inner` in `outer`, and `inner.score <= outer.score`
/// (the test of `reap_keep` in `kappa.rs`, without the context and frame).
#[inline]
fn contains(outer: &GappedHsp, inner: &GappedHsp) -> bool {
    outer.q_start <= inner.q_start
        && outer.q_end >= inner.q_start
        && outer.s_start <= inner.s_start
        && outer.s_end >= inner.s_start
        && outer.q_start <= inner.q_end
        && outer.q_end >= inner.q_end
        && outer.s_start <= inner.s_end
        && outer.s_end >= inner.s_end
        && inner.score <= outer.score
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:224-273
// ```c
// for (iread = 1;  iread < hsp_list->hspcnt;  iread++) {
//     ...
//     for (ihsp = 0;  ihsp < iread;  ihsp++) {
//         ...
//         if (hsp2 == NULL) continue;
//         if (hsp2->query.frame == hsp1->query.frame &&
//             hsp2->subject.frame == hsp1->subject.frame &&
//             CONTAINED_IN_HSP(q2start, q2end, q1start, s2start, s2end, s1start) &&
//             CONTAINED_IN_HSP(q2start, q2end, q1end, s2start, s2end, s1end) &&
//             hsp1->score <= hsp2->score) {
//             hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
//             break;
// ```
// The same test, over the earlier HSPs still kept of the same context and frame whose subject start
// lies in [max(s_start, s_end) - dmax, min(s_start, s_end)] only (dmax: the largest `s_end -
// s_start` of an HSP kept so far, at least 0): an HSP `outer` that contains both end points has
// s_start <= both of them and s_end >= both, so 0 <= s_end - s_start <= dmax and s_start >=
// max(s_start, s_end) of the HSP - dmax. Every HSP in that range is tested with the whole
// condition, the score included.
/// Whether each HSP of `hsps` is kept by the reap (`reap_keep` in `kappa.rs`), from the index.
#[inline(never)]
pub(super) fn x_reap_keep(hsps: &[(usize, GappedHsp)]) -> Vec<bool> {
    let mut keep = vec![true; hsps.len()];
    // (context, frame, subject start, position in `hsps`) of each HSP kept so far.
    let mut by_start: BTreeSet<(usize, i8, i32, u32)> = BTreeSet::new();
    let mut dmax = 0i64;
    for (read, &(context, inner)) in hsps.iter().enumerate() {
        let lo = (i64::from(inner.s_start.max(inner.s_end)) - dmax)
            .clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
        let hi = inner.s_start.min(inner.s_end);
        // Nearest start first: a containing HSP starts at or before this one, and the answer does
        // not depend on the order.
        let contained = lo <= hi
            && by_start
                .range((context, inner.frame, lo, 0)..=(context, inner.frame, hi, u32::MAX))
                .rev()
                .any(|&(_, _, _, earlier)| contains(&hsps[earlier as usize].1, &inner));
        if contained {
            keep[read] = false;
        } else {
            let position = u32::try_from(read).expect("fewer than 2^32 HSPs per match");
            by_start.insert((context, inner.frame, inner.s_start, position));
            dmax = dmax.max(i64::from(inner.s_end) - i64::from(inner.s_start));
        }
    }
    keep
}

#[cfg(test)]
mod tests {
    use super::super::reap_keep;
    use super::*;
    use std::cell::Cell;

    thread_local! {
        /// `Some(on)`: `LOSAT_X_TBNREAPIDX` of the reaps of this test thread.
        pub(super) static X_ON: Cell<Option<bool>> = const { Cell::new(None) };
    }

    /// xorshift64*: the test's own random numbers.
    struct Random(u64);
    impl Random {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn range(&mut self, low: i32, high: i32) -> i32 {
            low + (self.next() % (high - low + 1) as u64) as i32
        }
    }

    fn hsp(frame: i8, score: i32, q: (i32, i32), s: (i32, i32)) -> GappedHsp {
        GappedHsp {
            frame,
            score,
            q_start: q.0,
            q_end: q.1,
            q_gapped_start: q.0,
            s_start: s.0,
            s_end: s.1,
            s_gapped_start: s.0,
        }
    }

    /// The shape of one random list.
    struct Shape {
        len: usize,
        contexts: usize,
        frames: &'static [i8],
        /// Query and subject coordinates lie in 0..=span.
        span: i32,
        /// The longest HSP, in either sequence.
        longest: i32,
        /// Scores lie in 1..=scores (few values: ties).
        scores: i32,
        /// NCBI's order (`Blast_HSPListSortByScore`: score down) or none.
        by_score: bool,
        /// Some HSPs are copies of an earlier one shrunk at either end (contained ones).
        nested: bool,
        /// Some HSPs have an end before their start.
        reversed: bool,
    }

    fn random_list(random: &mut Random, shape: &Shape) -> Vec<(usize, GappedHsp)> {
        let mut list: Vec<(usize, GappedHsp)> = Vec::with_capacity(shape.len);
        while list.len() < shape.len {
            let context = (random.next() % shape.contexts as u64) as usize;
            let frame = shape.frames[(random.next() % shape.frames.len() as u64) as usize];
            let score = random.range(1, shape.scores);
            let row = if shape.nested && !list.is_empty() && random.next() % 2 == 0 {
                // Inside (or at the edges of) an earlier HSP, any score, sometimes another frame
                // or context.
                let (outer_context, outer) = list[(random.next() % list.len() as u64) as usize];
                let shrink = |random: &mut Random, start: i32, end: i32| {
                    let a = random.range(start, end.max(start));
                    let b = random.range(a, end.max(a));
                    match random.next() % 4 {
                        0 => (start, end),
                        1 => (start, b),
                        2 => (a, end),
                        _ => (a, b),
                    }
                };
                let q = shrink(random, outer.q_start, outer.q_end);
                let s = shrink(random, outer.s_start, outer.s_end);
                let mut frame = outer.frame;
                let mut context = outer_context;
                match random.next() % 8 {
                    0 => frame = -frame,
                    1 => context = (context + 1) % shape.contexts,
                    _ => {}
                }
                // One past either end now and then: not contained.
                let s = match random.next() % 8 {
                    0 => (s.0 - 1, s.1),
                    1 => (s.0, s.1 + 1),
                    _ => s,
                };
                (context, hsp(frame, score, q, s))
            } else {
                let q0 = random.range(0, shape.span);
                let s0 = random.range(0, shape.span);
                let q = (q0, (q0 + random.range(0, shape.longest)).min(shape.span));
                let s = (s0, (s0 + random.range(0, shape.longest)).min(shape.span));
                (context, hsp(frame, score, q, s))
            };
            let row = if shape.reversed && random.next() % 6 == 0 {
                let (context, mut h) = row;
                std::mem::swap(&mut h.s_start, &mut h.s_end);
                if random.next() % 2 == 0 {
                    std::mem::swap(&mut h.q_start, &mut h.q_end);
                }
                (context, h)
            } else {
                row
            };
            list.push(row);
        }
        if shape.by_score {
            // Stable, as NCBI's comparison breaks ties by position (`score_compare`).
            list.sort_by(|a, b| b.1.score.cmp(&a.1.score));
        }
        list
    }

    fn shapes() -> Vec<Shape> {
        let mut shapes = Vec::new();
        for &len in &[2, 5, 31, 32, 33, 64, 300, 800] {
            for &(span, longest) in &[(40, 12), (400, 60), (20_000, 400)] {
                for &by_score in &[true, false] {
                    for &nested in &[true, false] {
                        shapes.push(Shape {
                            len,
                            contexts: if len % 2 == 0 { 1 } else { 3 },
                            frames: if span == 40 {
                                &[1, -1]
                            } else {
                                &[1, 2, 3, -1, -2, -3]
                            },
                            span,
                            longest,
                            scores: if span == 40 { 4 } else { 200 },
                            by_score,
                            nested,
                            reversed: !by_score,
                        });
                    }
                }
            }
        }
        shapes
    }

    /// The index keeps exactly the HSPs NCBI's loop keeps, on random lists (dense and sparse
    /// coordinates, score order and none, ties, nested HSPs, ends one past the edges, other frames
    /// and contexts, ends before starts) of every length around `X_REAP_INDEX_FROM`.
    #[test]
    fn x_reap_keep_keeps_what_ncbi_keeps() {
        let mut random = Random(0x9E37_79B9_7F4A_7C15);
        let mut reaped = 0usize;
        for shape in shapes() {
            for round in 0..12 {
                let list = random_list(&mut random, &shape);
                let expected = reap_keep(&list);
                assert_eq!(
                    x_reap_keep(&list),
                    expected,
                    "len {} span {} by_score {} nested {} round {round}",
                    shape.len,
                    shape.span,
                    shape.by_score,
                    shape.nested
                );
                reaped += expected.iter().filter(|&&kept| !kept).count();
            }
        }
        // The lists do reap.
        assert!(reaped > 10_000, "reaped {reaped}");
    }

    /// Edges: a list of one, equal HSPs (the first kept), equal scores, end points on the edges of
    /// the containing HSP, one past them, another frame or context, and an HSP longer than every
    /// other kept (dmax grows after a short one was kept).
    #[test]
    fn x_reap_keep_edges() {
        let one = vec![(0, hsp(1, 5, (0, 10), (0, 10)))];
        assert_eq!(x_reap_keep(&one), vec![true]);
        let rows = vec![
            (0, hsp(1, 50, (100, 120), (500, 520))),
            (0, hsp(1, 50, (100, 120), (500, 520))),
            (0, hsp(1, 50, (100, 120), (500, 521))),
            (0, hsp(1, 49, (100, 120), (499, 520))),
            (0, hsp(2, 40, (105, 110), (505, 510))),
            (1, hsp(1, 40, (105, 110), (505, 510))),
            (0, hsp(1, 51, (105, 110), (505, 510))),
            (0, hsp(1, 30, (120, 120), (520, 520))),
            (0, hsp(1, 30, (0, 400), (0, 2_000))),
            (0, hsp(1, 20, (10, 20), (1_000, 1_010))),
            (0, hsp(1, 20, (10, 20), (1_999, 2_001))),
        ];
        let expected = vec![
            true, false, true, true, true, true, true, false, true, false, true,
        ];
        assert_eq!(reap_keep(&rows), expected);
        assert_eq!(x_reap_keep(&rows), expected);
    }

    /// `reap_contained_postredo_hsps` gives the same list and survivors with the switch on and off,
    /// below and from `X_REAP_INDEX_FROM`, and takes the index only from there.
    #[test]
    fn the_reap_is_the_same_with_and_without_the_index() {
        let mut random = Random(0x2545_F491_4F6C_DD1D);
        for len in [0, 1, 2, X_REAP_INDEX_FROM - 1, X_REAP_INDEX_FROM, 200] {
            let shape = Shape {
                len,
                contexts: 1,
                frames: &[1, -1],
                span: 60,
                longest: 20,
                scores: 6,
                by_score: true,
                nested: true,
                reversed: false,
            };
            let list = random_list(&mut random, &shape);
            let mut results = Vec::new();
            for on in [false, true] {
                X_ON.set(Some(on));
                assert_eq!(x_indexed(len), on && len >= X_REAP_INDEX_FROM);
                let mut rows = list.clone();
                let survivors = super::super::reap_contained_postredo_hsps(&mut rows);
                results.push((rows, survivors));
            }
            X_ON.set(None);
            assert_eq!(results[0], results[1], "len {len}");
        }
    }
}
