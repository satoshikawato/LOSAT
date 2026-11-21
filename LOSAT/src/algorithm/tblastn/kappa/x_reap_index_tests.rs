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
        let row = if shape.nested && !list.is_empty() && random.next().is_multiple_of(2) {
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
        let row = if shape.reversed && random.next().is_multiple_of(6) {
            let (context, mut h) = row;
            let which = random.next() % 3;
            if which != 1 {
                std::mem::swap(&mut h.s_start, &mut h.s_end);
            }
            if which != 0 {
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
/// the containing HSP, one past them, another frame or context, an HSP longer than every
/// other kept (dmax grows after a short one was kept), and subject coordinates at the ends of
/// `i32` (the range start below `i32::MIN` is clamped).
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
    let ends = vec![
        (0, hsp(1, 90, (0, 100), (i32::MIN, i32::MAX))),
        (0, hsp(1, 80, (10, 20), (i32::MAX - 5, i32::MAX - 1))),
        (0, hsp(1, 80, (10, 20), (i32::MIN, i32::MIN + 3))),
        (0, hsp(1, 95, (10, 20), (0, 5))),
        (0, hsp(-1, 10, (10, 20), (0, 5))),
    ];
    let expected = vec![true, false, false, true, true];
    assert_eq!(reap_keep(&ends), expected);
    assert_eq!(x_reap_keep(&ends), expected);
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
