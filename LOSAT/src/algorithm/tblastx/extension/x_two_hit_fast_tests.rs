use super::super::query_batch::InBatch;
use super::*;
use crate::algorithm::tblastx::scoring::{Blosum62Score, MatrixScore, TbxScore};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// `len` residues from `lo..hi`.
fn random_seq(rng: &mut Rng, len: usize, lo: u8, hi: u8) -> Vec<u8> {
    (0..len)
        .map(|_| lo + rng.below((hi - lo) as u64) as u8)
        .collect()
}

/// A copy of `base` with substitutions of probability `1/rate`.
fn mutate(rng: &mut Rng, base: &[u8], rate: u64) -> Vec<u8> {
    base.iter()
        .map(|&c| {
            if rng.below(rate) != 0 {
                c
            } else {
                1 + rng.below(24) as u8
            }
        })
        .collect()
}

/// One call both ways; the fast path must handle it.
fn check_call(q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32) {
    let expect =
        extend_hit_two_hit_reference(Blosum62Score, q, InBatch::ALONE, s, sl, sr, qr, d, false);
    let got = two_hit(q, s, sl, sr, qr, d).expect("the fast path handles the call");
    assert_eq!(got, expect, "sl={sl} sr={sr} qr={qr} d={d}");
}

fn check_pair(q: &[u8], s: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
    if q.len() < WORD || s.len() < WORD {
        return 0;
    }
    for _ in 0..calls {
        let d = 1 + rng.below(dmax as u64) as i32;
        let qr = rng.below((q.len() - WORD + 1) as u64) as usize;
        let sr = rng.below((s.len() - WORD + 1) as u64) as usize;
        // mostly at or before R; sometimes after it, or wrapped (`last_hit + wordsize` as usize)
        let sl = match rng.below(16) {
            0 => sr + 1 + rng.below(3) as usize,
            1 => usize::MAX - rng.below(3) as usize,
            _ => sr.saturating_sub(rng.below(60) as usize),
        };
        check_call(q, s, sl, sr, qr, d);
    }
    calls as u64
}

// No NCBI counterpart (two LOSAT paths). Random and related sequences (long walks), every
// code 1..=27, short and long, X-drops up to the largest the table is exact for.
#[test]
fn fast_two_hit_equals_the_reference_on_random_inputs() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut compared = 0u64;
    for round in 0..400 {
        let len = match round % 4 {
            0 => WORD + rng.below(40) as usize,
            1 => 20 + rng.below(120) as usize,
            2 => 100 + rng.below(600) as usize,
            _ => 600 + rng.below(3000) as usize,
        };
        let q = random_seq(&mut rng, len, 1, 28);
        let rate = 2 + rng.below(12);
        let s_related = mutate(&mut rng, &q, rate);
        let extra = rng.below(30) as usize;
        let s_random = random_seq(&mut rng, len + extra, 1, 28);
        let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
        compared += check_pair(&q, &s_related, &mut rng, 60, dmax);
        compared += check_pair(&q, &s_random, &mut rng, 60, dmax);
        compared += check_pair(
            &s_related[..(len / 2 + WORD).min(len)],
            &q,
            &mut rng,
            20,
            dmax,
        );
    }
    assert!(compared > 50_000);
}

// No NCBI counterpart. Two-letter alphabets give ties for the best and drops equal to the
// X-drop; a gap or sentinel (code 0) inside a sequence stops a walk as BLAST_SCORE_MIN does;
// codes 28..=255 score as the reference scores them.
#[test]
fn fast_two_hit_equals_the_reference_at_ties_sentinels_and_odd_codes() {
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    for (a, b) in [
        (1u8, 2u8),
        (1, 20),
        (4, 5),
        (11, 12),
        (19, 25),
        (5, 40),
        (27, 255),
    ] {
        for _ in 0..300 {
            let len = WORD + rng.below(200) as usize;
            let pick = |rng: &mut Rng| if rng.below(2) == 0 { a } else { b };
            let q: Vec<u8> = (0..len).map(|_| pick(&mut rng)).collect();
            let s: Vec<u8> = (0..len + rng.below(10) as usize)
                .map(|_| pick(&mut rng))
                .collect();
            let d = 1 + rng.below(12) as i32;
            check_pair(&q, &s, &mut rng, 20, d);
            let mut inner = s.clone();
            let at = rng.below(inner.len() as u64) as usize;
            inner[at] = 0;
            check_pair(&q, &inner, &mut rng, 20, 127);
        }
    }
}

// No NCBI counterpart. The calls the fast path leaves to the reference: X-drops outside
// 1..=127 and a word that does not fit.
#[test]
fn fast_two_hit_refuses_what_the_table_cannot_do() {
    let q = vec![5u8; 50];
    let s = vec![5u8; 50];
    assert!(two_hit(&q, &s, 10, 20, 20, 0).is_none());
    assert!(two_hit(&q, &s, 10, 20, 20, 128).is_none());
    check_call(&q, &s, 10, 20, 20, 127);
    assert!(two_hit(&q, &s, 10, 20, 48, 16).is_none());
    assert!(two_hit(&q, &s, 10, 48, 20, 16).is_none());
    check_call(&q, &s, 10, 47, 47, 16);
}

// No NCBI counterpart. Every cell of the table is `get_score` of its two bytes, with
// BLAST_SCORE_MIN as -128.
#[test]
fn the_table_holds_get_score_for_every_two_bytes() {
    let tab = build_table().expect("BLOSUM62 fits the byte table");
    for q in 0..=255u8 {
        for s in 0..=255u8 {
            let v = get_score(q, s);
            let want = if v == BLAST_SCORE_MIN { -128 } else { v };
            assert_eq!(
                tab[(q as usize) << 8 | s as usize] as i32,
                want,
                "q={q} s={s}"
            );
        }
    }
}

// No NCBI counterpart. Every matrix keeps a byte table for the walks at every word length of
// NCBI's TBLASTX (2..=4): its scores other than BLAST_SCORE_MIN lie in -127..=MAX_MATRIX_SCORE,
// every cell of the bytes below 32 is the score the reference walks read (`TbxScore::score`),
// with BLAST_SCORE_MIN as -128, code 0 stops the sentinel walks, a matrix has one table, and the
// search's instance takes the walks at each of those word lengths (`x_sent_walks`,
// `table_two_hit`); another word length takes the reference.
#[test]
fn every_matrix_keeps_its_walk_table_at_every_word_length() {
    for matrix in ScoringMatrix::ALL {
        let scores = MatrixScore::new(matrix);
        let sc = scores;
        for q in 0..MATRIX_LETTERS as u8 {
            for s in 0..MATRIX_LETTERS as u8 {
                let v = sc.score(q, s);
                assert!(
                    v == BLAST_SCORE_MIN || (-127..=MAX_MATRIX_SCORE).contains(&v),
                    "{matrix:?} {q} {s} {v}"
                );
            }
        }
        let t = matrix_walk_table(matrix, |a, b| sc.score(a, b)).expect("the matrix fits");
        for q in 0..MATRIX_LETTERS as u8 {
            for s in 0..MATRIX_LETTERS as u8 {
                let v = sc.score(q, s);
                let want = if v == BLAST_SCORE_MIN { -128 } else { v };
                assert_eq!(t.at(q, s), want, "{matrix:?} {q} {s}");
            }
        }
        let again = matrix_walk_table(matrix, |_, _| unreachable!("built once"));
        assert!(std::ptr::eq(t, again.expect("the matrix fits")));
        assert!(t.x_sentinels_stop(), "{matrix:?}");
        let q = vec![5u8; 50];
        for word in 2..=MAX_MATRIX_WORD {
            let scores = MatrixScore::with_word(matrix, word);
            let sc = scores;
            assert!(sc.x_sent_walks(), "{matrix:?} {word}");
            assert!(sc.table_two_hit(&q, &q, 10, 20, 20, 16).is_some());
        }
        for word in [1, 5, 7] {
            let scores = MatrixScore::with_word(matrix, word);
            let sc = scores;
            assert!(!sc.x_sent_walks(), "{matrix:?} {word}");
            assert!(sc.table_two_hit(&q, &q, 10, 20, 20, 16).is_none());
        }
    }
    // A score above MAX_MATRIX_SCORE (or below -127) leaves the calls to the reference.
    let mut t = [0i8; MATRIX_TABLE_LEN];
    assert!(fill_matrix_table(&mut t, |a, b| if a == 3 && b == 3 {
        MAX_MATRIX_SCORE
    } else {
        1
    }));
    assert!(!fill_matrix_table(&mut t, |a, b| if a == 3 && b == 3 {
        MAX_MATRIX_SCORE + 1
    } else {
        1
    }));
    assert!(!fill_matrix_table(&mut t, |a, b| if a == 31 && b == 31 {
        -128
    } else {
        1
    }));
}

/// One call of the walks with words of `W` letters on the byte table `t` of `sc`'s scores
/// against the reference with `sc` (`TbxScore::word` `W`): the table walks and the sentinel
/// walks must handle it. Returns the reference's values.
#[expect(
    clippy::too_many_arguments,
    reason = "the extension's six arguments, the scores and their table"
)]
fn check_walks<const W: usize, S: TbxScore>(
    sc: S,
    t: &MatrixWalkTable,
    q: &[u8],
    s: &[u8],
    sl: usize,
    sr: usize,
    qr: usize,
    d: i32,
) -> TwoHit {
    assert_eq!(sc.word(), W);
    let expect = extend_hit_two_hit_reference(sc, q, InBatch::ALONE, s, sl, sr, qr, d, false);
    let got =
        two_hit_in_word::<W>(t, q, s, sl, sr, qr, d).expect("the table walks handle the call");
    assert_eq!(got, expect, "W={W} sl={sl} sr={sr} qr={qr} d={d}");
    assert!(t.x_sentinels_stop());
    let (qb, qb_range) = with_sentinels(q);
    let (sb, sb_range) = with_sentinels(s);
    // SAFETY: both buffers are [0, residues.., 0] and short; `t` stops at a sentinel.
    let sent = unsafe { two_hit_sent_in::<W>(t, &qb[qb_range], &sb[sb_range], sl, sr, qr, d) }
        .expect("the sentinel walks handle the call");
    assert_eq!(sent, expect, "sent W={W} sl={sl} sr={sr} qr={qr} d={d}");
    expect
}

/// A check of one call: (query, subject, s_left_off, s_right_off, q_right_off, x_drop).
type Check<'a> = dyn FnMut(&[u8], &[u8], usize, usize, usize, i32) + 'a;

/// `calls` calls of `check` with random X-drops up to `dmax`, words of `W` letters anywhere and
/// often at either end of either sequence, and first hits mostly at or before R, sometimes
/// after it, or wrapped (`last_hit + wordsize` as usize).
fn walk_calls<const W: usize>(
    q: &[u8],
    s: &[u8],
    rng: &mut Rng,
    calls: usize,
    dmax: i32,
    check: &mut Check<'_>,
) -> u64 {
    if q.len() < W || s.len() < W {
        return 0;
    }
    for _ in 0..calls {
        let d = 1 + rng.below(dmax as u64) as i32;
        let qr = match rng.below(8) {
            0 => 0,
            1 => q.len() - W,
            _ => rng.below((q.len() - W + 1) as u64) as usize,
        };
        let sr = match rng.below(8) {
            0 => 0,
            1 => s.len() - W,
            _ => rng.below((s.len() - W + 1) as u64) as usize,
        };
        let sl = match rng.below(16) {
            0 => sr + 1 + rng.below(3) as usize,
            1 => usize::MAX - rng.below(3) as usize,
            _ => sr.saturating_sub(rng.below(60) as usize),
        };
        check(q, s, sl, sr, qr, d);
    }
    calls as u64
}

/// `check` on random and related sequences (long walks; the subject shorter than the query and
/// the other way round), every code 1..=27, short and long, X-drops up to 127; then on the
/// two-letter alphabets `ties` (ties for the best, drops equal to the X-drop) with a code 0
/// inside a sequence. Returns the number of calls.
fn walk_rounds<const W: usize>(
    rng: &mut Rng,
    rounds: usize,
    ties: &[(u8, u8)],
    check: &mut Check<'_>,
) -> u64 {
    let mut n = 0u64;
    for round in 0..rounds {
        let len = match round % 4 {
            0 => W + rng.below(40) as usize,
            1 => 20 + rng.below(120) as usize,
            2 => 100 + rng.below(600) as usize,
            _ => 600 + rng.below(2000) as usize,
        };
        let q = random_seq(rng, len, 1, 28);
        let rate = 2 + rng.below(40);
        let s_related = mutate(rng, &q, rate);
        let extra = rng.below(30) as usize;
        let s_random = random_seq(rng, len + extra, 1, 28);
        let dmax = [4, 16, 40, 127][round / 4 % 4];
        n += walk_calls::<W>(&q, &s_related, rng, 30, dmax, check);
        n += walk_calls::<W>(&q, &s_random, rng, 30, dmax, check);
        let shorter = &s_related[..(len / 2 + W).min(len)];
        n += walk_calls::<W>(shorter, &q, rng, 10, dmax, check);
        let longer = &s_related[(len / 3).min(len - W)..];
        n += walk_calls::<W>(&q, longer, rng, 10, dmax, check);
    }
    for &(a, b) in ties {
        for _ in 0..40 {
            let len = W + rng.below(200) as usize;
            let pick = |rng: &mut Rng| if rng.below(2) == 0 { a } else { b };
            let q: Vec<u8> = (0..len).map(|_| pick(rng)).collect();
            let mut s: Vec<u8> = (0..len + rng.below(10) as usize)
                .map(|_| pick(rng))
                .collect();
            let d = 1 + rng.below(20) as i32;
            n += walk_calls::<W>(&q, &s, rng, 10, d, check);
            let at = rng.below(s.len() as u64) as usize;
            s[at] = 0;
            n += walk_calls::<W>(&q, &s, rng, 10, 127, check);
            n += walk_calls::<W>(&s, &q, rng, 10, 127, check);
        }
    }
    n
}

/// `check` on sequences of the letter `top` (the largest score against itself): a code 0 at
/// every position of the query and of the subject in turn, around a word of `W` letters (in the
/// word, at the first left step R - 1, at the first right step R, further out), so that the word
/// loop adds a -128 to up to `W - 1` of the largest scores; and identical sequences, where both
/// walks run out at the ends (the sentinels of the sentinel walks, the reference's bound), from
/// every word position, with the subject offset from the query on both sides. First hits at,
/// before and after R; X-drops 1, 11, 127.
fn walk_edges<const W: usize>(top: u8, check: &mut Check<'_>) {
    let len = 48;
    for at in 0..len {
        for in_query in [true, false] {
            let mut q = vec![top; len];
            let mut s = vec![top; len];
            if in_query {
                q[at] = 0;
            } else {
                s[at] = 0;
            }
            for (qr, sr) in [(20, 20), (20, 24), (0, 0), (len - W, len - W), (3, 40)] {
                for sl in [sr.saturating_sub(10), sr, sr + 2, usize::MAX] {
                    for d in [1, 11, 127] {
                        check(&q, &s, sl, sr, qr, d);
                    }
                }
            }
        }
    }
    for (q_len, s_len) in [(W, W), (W, 9), (9, W), (17, 17), (25, 14), (14, 25)] {
        let q = vec![top; q_len];
        let s = vec![top; s_len];
        for qr in 0..=q_len - W {
            for sr in 0..=s_len - W {
                for sl in [sr.saturating_sub(5), sr, sr + 1, usize::MAX] {
                    for d in [1, 11, 127] {
                        check(&q, &s, sl, sr, qr, d);
                    }
                }
            }
        }
    }
}

/// The table walks and the sentinel walks of `matrix` with words of `W` letters against the
/// reference with `MatrixScore::with_word(matrix, W)`, directly and through the search's
/// dispatch (`TbxScore::table_two_hit`, `TbxScore::sent_two_hit`); the calls they leave to the
/// reference (X-drops outside 1..=127, a word that does not fit) are left.
fn matrix_walks_equal_the_reference<const W: usize>(matrix: ScoringMatrix, rng: &mut Rng) {
    let scores = MatrixScore::with_word(matrix, W);
    let sc = scores;
    assert_eq!(sc.word(), W);
    assert!(sc.x_sent_walks(), "{matrix:?}");
    let t = matrix_walk_table(matrix, |a, b| sc.score(a, b)).expect("the matrix fits");
    let mut check = |q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32| {
        let expect = check_walks::<W, _>(sc, t, q, s, sl, sr, qr, d);
        let what = format!("{matrix:?} W={W} sl={sl} sr={sr} qr={qr} d={d}");
        assert_eq!(
            sc.table_two_hit(q, s, sl, sr, qr, d),
            Some(expect),
            "{what}"
        );
        let (qb, qb_range) = with_sentinels(q);
        let (sb, sb_range) = with_sentinels(s);
        // SAFETY: both buffers are [0, residues.., 0] and short; `x_sent_walks` holds.
        let (qf, sf) = (&qb[qb_range], &sb[sb_range]);
        let sent = unsafe { sc.sent_two_hit(qf, InBatch::ALONE, sf, sl, sr, qr, d) };
        assert_eq!(sent, expect, "sent {what}");
        // SAFETY: as above, and the words have `W` letters (the unit path's walks).
        let sent = unsafe { sc.sent_two_hit_of::<W>(qf, InBatch::ALONE, sf, sl, sr, qr, d) };
        assert_eq!(sent, expect, "sent_two_hit_of {what}");
    };
    let ties = [(1u8, 2u8), (4, 5), (19, 25), (5, 31), (27, 29)];
    let compared = walk_rounds::<W>(rng, 60, &ties, &mut check);
    assert!(compared > 10_000, "{matrix:?} W={W} {compared}");
    let top = (1..28u8)
        .max_by_key(|&c| (sc.score(c, c), std::cmp::Reverse(c)))
        .expect("letters");
    walk_edges::<W>(top, &mut check);

    // The calls the walks leave to the reference: X-drops outside 1..=127, a word that does
    // not fit (the last word that fits is taken).
    let (qb, qb_range) = with_sentinels(&[5u8; 50]);
    let (q, q_full) = (&qb[qb_range.start + 1..qb_range.end - 1], &qb[qb_range]);
    for (sr, qr, d) in [
        (20, 20, 0),
        (20, 20, 128),
        (20, 51 - W, 16),
        (51 - W, 20, 16),
    ] {
        assert!(two_hit_in_word::<W>(t, q, q, 10, sr, qr, d).is_none());
        assert!(sc.table_two_hit(q, q, 10, sr, qr, d).is_none());
        // SAFETY: the buffer is [0, residues.., 0].
        assert!(unsafe { two_hit_sent_in::<W>(t, q_full, q_full, 10, sr, qr, d) }.is_none());
    }
    for (sr, qr, d) in [(20, 20, 0), (20, 20, 128), (50 - W, 50 - W, 16)] {
        let expect = extend_hit_two_hit_reference(sc, q, InBatch::ALONE, q, 10, sr, qr, d, false);
        let with = super::super::two_hit::extend_hit_two_hit_with(
            sc,
            q,
            InBatch::ALONE,
            q,
            10,
            sr,
            qr,
            d,
            false,
        );
        assert_eq!(with, expect, "{matrix:?} W={W} sr={sr} qr={qr} d={d}");
        // SAFETY: the buffer is [0, residues.., 0]; `x_sent_walks` holds.
        let sent = unsafe { sc.sent_two_hit(q_full, InBatch::ALONE, q_full, 10, sr, qr, d) };
        assert_eq!(sent, expect, "sent {matrix:?} W={W} sr={sr} qr={qr} d={d}");
        // SAFETY: as above, and the words have `W` letters.
        let sent =
            unsafe { sc.sent_two_hit_of::<W>(q_full, InBatch::ALONE, q_full, 10, sr, qr, d) };
        assert_eq!(
            sent, expect,
            "sent_two_hit_of {matrix:?} W={W} sr={sr} qr={qr} d={d}"
        );
    }
}

// No NCBI counterpart (two LOSAT paths). The table walks and the sentinel walks of every
// matrix with words of 2 letters (`-word_size 2`) give the reference walks' values.
#[test]
fn matrix_walks_equal_the_reference_at_word_length_2() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    for matrix in ScoringMatrix::ALL {
        matrix_walks_equal_the_reference::<2>(matrix, &mut rng);
    }
}

// No NCBI counterpart (two LOSAT paths). The same with words of 3 letters (`-matrix` other
// than BLOSUM62; BLOSUM62's search takes `two_hit` and `two_hit_sent`).
#[test]
fn matrix_walks_equal_the_reference_at_word_length_3() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for matrix in ScoringMatrix::ALL {
        matrix_walks_equal_the_reference::<3>(matrix, &mut rng);
    }
}

// No NCBI counterpart (two LOSAT paths). The same with words of 4 letters (`-word_size 4`).
#[test]
fn matrix_walks_equal_the_reference_at_word_length_4() {
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    for matrix in ScoringMatrix::ALL {
        matrix_walks_equal_the_reference::<4>(matrix, &mut rng);
    }
}

/// The scores of a made-up matrix (`bound_cells`) for words of `W` letters.
#[derive(Clone, Copy)]
struct BoundScore<const W: usize>(&'static [[i32; MATRIX_LETTERS]; MATRIX_LETTERS]);

impl<const W: usize> TbxScore for BoundScore<W> {
    const BLOSUM62: bool = false;

    fn score(self, a: u8, b: u8) -> i32 {
        self.0[a as usize][b as usize]
    }

    // Not read by the walks (`TbxScore::score` only).
    fn matrix(self) -> ScoringMatrix {
        ScoringMatrix::Identity
    }

    fn word(self) -> usize {
        W
    }
}

/// The letter that scores `MAX_MATRIX_SCORE` against itself in `bound_cells`.
const BOUND_TOP: u8 = 20;

/// A made-up matrix at the walks' bounds: random scores in -MAX_MATRIX_SCORE..=MAX_MATRIX_SCORE,
/// `BOUND_TOP` against itself MAX_MATRIX_SCORE, -127 for 1 against 2, and BLAST_SCORE_MIN for
/// code 0 against anything.
fn bound_cells(rng: &mut Rng) -> &'static [[i32; MATRIX_LETTERS]; MATRIX_LETTERS] {
    let mut cells = [[0i32; MATRIX_LETTERS]; MATRIX_LETTERS];
    for (q, row) in cells.iter_mut().enumerate() {
        for (s, v) in row.iter_mut().enumerate() {
            *v = if q == 0 || s == 0 {
                BLAST_SCORE_MIN
            } else if q == BOUND_TOP as usize && s == q {
                MAX_MATRIX_SCORE
            } else if (q, s) == (1, 2) {
                -127
            } else {
                rng.below(2 * MAX_MATRIX_SCORE as u64 + 1) as i32 - MAX_MATRIX_SCORE
            };
        }
    }
    Box::leak(Box::new(cells))
}

/// `check_walks` of the made-up matrix with words of `W` letters.
fn bound_walks_equal_the_reference<const W: usize>(rng: &mut Rng) {
    let cells = bound_cells(rng);
    let mut t = Box::new(MatrixWalkTable([0; MATRIX_TABLE_LEN], [0; KEY_TABLE_LEN]));
    assert!(fill_matrix_table(&mut t.0, |a, b| cells[a as usize][b as usize]));
    let t0 = t.0;
    fill_keys::<5>(&mut t.1, &t0);
    let sc = BoundScore::<W>(cells);
    let mut check = |q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32| {
        check_walks::<W, _>(sc, &t, q, s, sl, sr, qr, d);
    };
    let ties = [(BOUND_TOP, 3u8), (1, 2), (BOUND_TOP, 1), (5, 31)];
    let compared = walk_rounds::<W>(rng, 80, &ties, &mut check);
    assert!(compared > 10_000, "W={W} {compared}");
    walk_edges::<W>(BOUND_TOP, &mut check);
}

// No NCBI counterpart (two LOSAT paths). The walks' argument at its bound (module comment):
// a matrix whose scores reach MAX_MATRIX_SCORE (42) and -127, with words of 2, 3 and 4
// letters, gives the reference walks' values, also where the word loop adds a code 0
// (BLAST_SCORE_MIN, read as -128) to `W - 1` scores of 42 (126 at 4 letters).
#[test]
fn walks_at_the_score_bound_equal_the_reference_at_every_word_length() {
    let mut rng = Rng(0x9E6C_63D0_676A_9A99);
    bound_walks_equal_the_reference::<2>(&mut rng);
    bound_walks_equal_the_reference::<3>(&mut rng);
    bound_walks_equal_the_reference::<4>(&mut rng);
}

/// `[W.., 0, seq.., 0, W..]` with `PAD` letters W around, and the range of `[0, seq.., 0]`:
/// a step of the sentinel walks past a sentinel would read a W (+11) and change the result.
const PAD: usize = 8;
fn with_sentinels(seq: &[u8]) -> (Vec<u8>, std::ops::Range<usize>) {
    let mut v = vec![20u8; PAD];
    v.push(0);
    v.extend_from_slice(seq);
    v.push(0);
    let r = PAD..v.len();
    v.extend(std::iter::repeat_n(20u8, PAD));
    (v, r)
}

/// One call of the sentinel walks against the reference; they must handle it.
fn check_sent(q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32) {
    let (qb, qr_range) = with_sentinels(q);
    let (sb, sr_range) = with_sentinels(s);
    let expect =
        extend_hit_two_hit_reference(Blosum62Score, q, InBatch::ALONE, s, sl, sr, qr, d, false);
    // SAFETY: both buffers are [0, residues.., 0] and short.
    let got = unsafe { two_hit_sent(&qb[qr_range], &sb[sr_range], sl, sr, qr, d) }
        .expect("the sentinel walks handle the call");
    assert_eq!(got, expect, "sl={sl} sr={sr} qr={qr} d={d}");
}

fn check_sent_pair(q: &[u8], s: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
    if q.len() < WORD || s.len() < WORD {
        return 0;
    }
    for _ in 0..calls {
        let d = 1 + rng.below(dmax as u64) as i32;
        // the word anywhere, and often at either end of either sequence
        let qr = match rng.below(8) {
            0 => 0,
            1 => q.len() - WORD,
            _ => rng.below((q.len() - WORD + 1) as u64) as usize,
        };
        let sr = match rng.below(8) {
            0 => 0,
            1 => s.len() - WORD,
            _ => rng.below((s.len() - WORD + 1) as u64) as usize,
        };
        let sl = match rng.below(16) {
            0 => sr + 1 + rng.below(3) as usize,
            1 => usize::MAX - rng.below(3) as usize,
            _ => sr.saturating_sub(rng.below(60) as usize),
        };
        check_sent(q, s, sl, sr, qr, d);
    }
    calls as u64
}

// No NCBI counterpart (two LOSAT paths): LOSAT_X_TBXUNGAPSENT against the reference on random
// and related sequences (walks that run to the sentinels at both ends), short and long, every
// code 1..=27, X-drops up to 127.
#[test]
fn sentinel_walks_equal_the_reference_on_random_inputs() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut compared = 0u64;
    for round in 0..400 {
        let len = match round % 4 {
            0 => WORD + rng.below(40) as usize,
            1 => 20 + rng.below(120) as usize,
            2 => 100 + rng.below(600) as usize,
            _ => 600 + rng.below(3000) as usize,
        };
        let q = random_seq(&mut rng, len, 1, 28);
        let rate = 2 + rng.below(40);
        let s_related = mutate(&mut rng, &q, rate);
        let extra = rng.below(30) as usize;
        let s_random = random_seq(&mut rng, len + extra, 1, 28);
        let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
        compared += check_sent_pair(&q, &s_related, &mut rng, 60, dmax);
        compared += check_sent_pair(&q, &s_random, &mut rng, 60, dmax);
        // the subject shorter than the query, and the other way round
        compared += check_sent_pair(
            &s_related[..(len / 2 + WORD).min(len)],
            &q,
            &mut rng,
            20,
            dmax,
        );
        compared += check_sent_pair(
            &q,
            &s_related[(len / 3).min(len - WORD)..],
            &mut rng,
            20,
            dmax,
        );
    }
    assert!(compared > 60_000);
}

// No NCBI counterpart. Identical sequences of W (+11 a step): both walks run out at the
// sentinels (the reference's bound), from every word position, with the subject offset from
// the query on both sides.
#[test]
fn sentinel_walks_stop_at_both_ends_as_the_bound_does() {
    for (q_len, s_len) in [(3, 3), (3, 9), (9, 3), (17, 17), (40, 25), (25, 40)] {
        let q = vec![20u8; q_len];
        let s = vec![20u8; s_len];
        for qr in 0..=q_len - WORD {
            for sr in 0..=s_len - WORD {
                for sl in [sr.saturating_sub(5), sr, sr + 1, usize::MAX] {
                    for d in [1, 11, 127] {
                        check_sent(&q, &s, sl, sr, qr, d);
                    }
                }
            }
        }
    }
}

// No NCBI counterpart. Two-letter alphabets (ties for the best, drops equal to the X-drop)
// and a code 0 inside a sequence (it stops a walk as BLAST_SCORE_MIN does), codes 28..=31 (the
// keys' bytes past the letters); a buffer with a byte of 32 or more does not take the path.
#[test]
fn sentinel_walks_equal_the_reference_at_ties_inner_sentinels_and_odd_codes() {
    let mut rng = Rng(0x9E6C_63D0_676A_9A99);
    for (a, b) in [(5u8, 40u8), (27, 255), (31, 32)] {
        let frame = [0u8, a, b, a, b, 0];
        assert!(!sent_frame(true, &frame), "{a} {b}");
        assert!(!sent_contexts(std::iter::once(&frame[..])), "{a} {b}");
        assert!(sent_frame(true, &[0u8, a, a, a, 0]) == (a < 32), "{a}");
    }
    for (a, b) in [
        (1u8, 2u8),
        (1, 20),
        (4, 5),
        (11, 12),
        (19, 25),
        (27, 31),
        (28, 30),
    ] {
        for _ in 0..300 {
            let len = WORD + rng.below(200) as usize;
            let pick = |rng: &mut Rng| if rng.below(2) == 0 { a } else { b };
            let q: Vec<u8> = (0..len).map(|_| pick(&mut rng)).collect();
            let s: Vec<u8> = (0..len + rng.below(10) as usize)
                .map(|_| pick(&mut rng))
                .collect();
            let d = 1 + rng.below(12) as i32;
            check_sent_pair(&q, &s, &mut rng, 20, d);
            let mut inner = s.clone();
            let at = rng.below(inner.len() as u64) as usize;
            inner[at] = 0;
            check_sent_pair(&q, &inner, &mut rng, 20, 127);
            check_sent_pair(&inner, &q, &mut rng, 20, 127);
        }
    }
}

// No NCBI counterpart. The calls the sentinel walks leave to the reference (as `two_hit`):
// X-drops outside 1..=127 and a word that does not fit.
#[test]
fn sentinel_walks_refuse_what_the_table_cannot_do() {
    let (qb, qr) = with_sentinels(&[5u8; 50]);
    let (sb, sr) = with_sentinels(&[5u8; 50]);
    let (q, s) = (&qb[qr.clone()], &sb[sr.clone()]);
    // SAFETY: both buffers are [0, residues.., 0].
    unsafe {
        assert!(two_hit_sent(q, s, 10, 20, 20, 0).is_none());
        assert!(two_hit_sent(q, s, 10, 20, 20, 128).is_none());
        assert!(two_hit_sent(q, s, 10, 20, 48, 16).is_none());
        assert!(two_hit_sent(q, s, 10, 48, 20, 16).is_none());
    }
    check_sent(&[5u8; 50], &[5u8; 50], 10, 20, 20, 127);
    check_sent(&[5u8; 50], &[5u8; 50], 10, 47, 47, 16);
    // `extend_hit_two_hit_sent` runs the reference for the X-drops (the callers' words fit).
    let (qv, sv) = (&qb[qr.start + 1..qr.end - 1], &sb[sr.start + 1..sr.end - 1]);
    for (sl, srr, qrr, d) in [(10, 20, 20, 0), (10, 20, 20, 128), (10, 47, 47, 0)] {
        let expect = extend_hit_two_hit_reference(
            Blosum62Score,
            qv,
            InBatch::ALONE,
            sv,
            sl,
            srr,
            qrr,
            d,
            false,
        );
        // SAFETY: both buffers are [0, residues.., 0].
        let got = unsafe {
            super::super::two_hit::extend_hit_two_hit_sent(q, InBatch::ALONE, s, sl, srr, qrr, d)
        };
        assert_eq!(got, expect, "sl={sl} sr={srr} qr={qrr} d={d}");
    }
}

// No NCBI counterpart. A chunk takes the sentinel walks only when its frame may (`sent_frame`),
// the chunk is the whole frame and the buffer has its two sentinels; otherwise the walks of
// `two_hit`.
#[test]
fn only_a_whole_frame_with_its_sentinels_takes_the_sentinel_walks() {
    let frame = [0u8, 4, 5, 6, 7, 0];
    assert!(sent_subject(true, &frame, 0, 4));
    assert!(!sent_subject(false, &frame, 0, 4));
    // a chunk of the frame (MAX_DBSEQ_LEN splits frames longer than 5 M letters)
    assert!(!sent_subject(true, &frame, 1, 3));
    assert!(!sent_subject(true, &frame, 0, 3));
    assert!(!sent_subject(true, &frame, 0, usize::MAX));
    assert!(!sent_subject(true, &[0u8, 4, 5, 6, 7, 8], 0, 4));
    assert!(!sent_subject(true, &[1u8, 4, 5, 6, 7, 0], 0, 4));
    // contexts: every buffer needs both sentinels (and the switch)
    let good: [&[u8]; 2] = [&[0, 4, 5, 6, 0], &[0, 0]];
    assert_eq!(sent_contexts(good.into_iter()), state().sent);
    let bad: [&[u8]; 2] = [&[0, 4, 5, 6, 0], &[0, 4, 5]];
    assert!(!sent_contexts(bad.into_iter()));
    assert!(!sent_contexts([&[0u8][..]].into_iter()));
}

// No NCBI counterpart. A gap or sentinel code (0) at every position of the query and of the
// subject in turn, around a word of W (the largest score, 11, which bounds the word-loop
// argument): in the word, at the first left step (R - 1), at the first right step (R), and
// further out; with the first hit at, before and after R, and X-drops 1, 11, 127.
#[test]
fn a_sentinel_anywhere_stops_the_walks_as_the_reference_does() {
    const W: u8 = 20;
    assert_eq!(get_score(W, W), MAX_WORD_SCORE);
    let len = 48;
    for at in 0..len {
        for in_query in [true, false] {
            let mut q = vec![W; len];
            let mut s = vec![W; len];
            if in_query {
                q[at] = 0;
            } else {
                s[at] = 0;
            }
            for (qr, sr) in [
                (20, 20),
                (20, 24),
                (0, 0),
                (len - WORD, len - WORD),
                (3, 40),
            ] {
                for sl in [sr.saturating_sub(10), sr, sr + 2, usize::MAX] {
                    for d in [1, 11, 127] {
                        check_call(&q, &s, sl, sr, qr, d);
                    }
                }
            }
        }
    }
}
