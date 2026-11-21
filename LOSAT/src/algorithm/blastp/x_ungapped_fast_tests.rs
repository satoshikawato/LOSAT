use super::*;
use crate::algorithm::blastp::extension::{extend_one_hit, extend_two_hit};
use crate::config::ScoringMatrix;

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

/// A sequence of `len` residues from `lo..hi` with the sentinel 0 at both ends.
fn random_seq(rng: &mut Rng, len: usize, lo: u8, hi: u8) -> Vec<u8> {
    let mut v = vec![0u8];
    for _ in 0..len {
        v.push(lo + rng.below((hi - lo) as u64) as u8);
    }
    v.push(0);
    v
}

/// A mutated copy of `base` (substitutions with probability `1/rate`).
fn mutate(rng: &mut Rng, base: &[u8], rate: u64) -> Vec<u8> {
    base.iter()
        .enumerate()
        .map(|(i, &c)| {
            if i == 0 || i + 1 == base.len() || rng.below(rate) != 0 {
                c
            } else {
                1 + rng.below(24) as u8
            }
        })
        .collect()
}

/// Both paths on one call (the reference through the matrix-generic entry points, which the
/// switch does not touch); the fast path must handle it.
fn check_call(q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32, w: usize) {
    let expect = extend_two_hit(ScoringMatrix::Blosum62, q, s, sl, sr, qr, d, w);
    let got = two_hit(q, s, sl, sr, qr, d, w).expect("fast path handles the call");
    assert_eq!(got, expect, "two_hit sl={sl} sr={sr} qr={qr} d={d} w={w}");
    let expect = extend_one_hit(ScoringMatrix::Blosum62, q, s, qr, sr, d, w);
    let got = one_hit(q, s, qr, sr, d, w).expect("fast path handles the call");
    assert_eq!(got, expect, "one_hit qr={qr} sr={sr} d={d} w={w}");
}

fn check_pair(q: &[u8], s: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
    for _ in 0..calls {
        let w = 1 + rng.below(4) as usize;
        let d = 1 + rng.below(dmax as u64) as i32;
        let qr = rng.below(q.len() as u64) as usize;
        let sr = rng.below(s.len() as u64) as usize;
        let sl = if sr > 0 {
            rng.below(sr as u64) as usize
        } else {
            0
        };
        check_call(q, s, sl, sr, qr, d, w);
    }
    2 * calls as u64
}

/// Random sequences, related sequences (long walks), all codes, short and long.
#[test]
fn fast_extensions_equal_the_reference_on_random_inputs() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut compared = 0u64;
    for round in 0..400 {
        let len = match round % 4 {
            0 => 1 + rng.below(40) as usize,
            1 => 20 + rng.below(120) as usize,
            2 => 100 + rng.below(600) as usize,
            _ => 600 + rng.below(3000) as usize,
        };
        // residues 1..=27, so every code except the sentinel appears
        let q = random_seq(&mut rng, len, 1, 28);
        let rate = 2 + rng.below(12);
        let s_related = mutate(&mut rng, &q, rate);
        let extra = rng.below(30) as usize;
        let s_random = random_seq(&mut rng, len + extra, 1, 28);
        let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
        compared += check_pair(&q, &s_related, &mut rng, 60, dmax);
        compared += check_pair(&q, &s_random, &mut rng, 60, dmax);
        // a sentinel inside the sequence (between contexts)
        let mut inner = q.clone();
        let at = 1 + rng.below(len as u64) as usize;
        inner[at] = 0;
        compared += check_pair(&inner, &s_related, &mut rng, 40, dmax);
    }
    assert!(compared > 100_000);
}

/// X-drop at the boundary and ties for the best position: walks whose score returns to the
/// best value exactly, drops that equal the dropoff, repeated patterns.
#[test]
fn fast_extensions_equal_the_reference_at_boundaries_and_ties() {
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let mut compared = 0u64;
    // Two-letter alphabets give plus, zero and minus steps in equal measure, so the score
    // returns to the same value over and over (ties) and meets best - score == dropoff.
    for pairs in [
        (1u8, 2u8),
        (1, 20),
        (20, 23),
        (4, 5),
        (13, 14),
        (1, 4),
        (17, 20),
    ] {
        for round in 0..300 {
            let len = 1 + rng.below(if round % 3 == 0 { 50 } else { 400 }) as usize;
            let mut draw = |n: usize| -> Vec<u8> {
                std::iter::once(0)
                    .chain((0..n).map(|_| if rng.below(2) == 0 { pairs.0 } else { pairs.1 }))
                    .chain(std::iter::once(0))
                    .collect()
            };
            let q = draw(len);
            let s = draw(len + 3);
            compared += check_pair(&q, &s, &mut rng, 30, 25);
        }
    }
    // Every dropoff 1..=60 on a periodic pattern: +11 (W) among +4 (A-A) steps.
    let w_code = 20u8;
    let q: Vec<u8> = [0u8]
        .into_iter()
        .chain((0..300).map(|i| if i % 7 == 3 { w_code } else { 1 }))
        .chain([0u8])
        .collect();
    let s: Vec<u8> = [0u8]
        .into_iter()
        .chain((0..300).map(|i| if i % 5 == 1 { w_code } else { 1 }))
        .chain([0u8])
        .collect();
    for d in 1..=60 {
        for sr in 1..q.len() - 4 {
            for w in [1usize, 3] {
                check_call(&q, &s, sr.saturating_sub(5), sr, sr, d, w);
                compared += 2;
            }
        }
    }
    assert!(compared > 50_000);
}

/// Very long walks: identical sequences whose score passes 16 bits, with the sentinels at both
/// ends, and walks that start next to either end.
#[test]
fn fast_extensions_equal_the_reference_on_very_long_walks() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    // 9000 tryptophans: 99,000 > 2^16 (W-W = 11)
    let mut q = vec![0u8];
    q.extend(std::iter::repeat_n(20u8, 9000));
    q.push(0);
    let mut s = q.clone();
    for d in [1, 5, 16, 100, 127] {
        for (qr, sr) in [
            (1usize, 1usize),
            (4000, 4000),
            (8990, 8990),
            (5000, 20),
            (20, 5000),
        ] {
            check_call(&q, &s, sr.saturating_sub(10), sr, qr, d, 3);
        }
    }
    // a long walk with mismatches that never reach the dropoff
    for i in (0..s.len()).step_by(9) {
        if s[i] != 0 {
            s[i] = 1;
        }
    }
    for _ in 0..200 {
        let qr = 1 + rng.below(8990) as usize;
        let sr = 1 + rng.below(8990) as usize;
        let sl = sr.saturating_sub(rng.below(30) as usize);
        let d = 20 + rng.below(100) as i32;
        check_call(&q, &s, sl, sr, qr, d, 3);
    }
}

/// Words of 1..=11 residues with sentinels and the extreme scores inside them, and words next
/// to either end of a sequence.
#[test]
fn fast_extensions_equal_the_reference_for_sentinels_in_words() {
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    for _ in 0..20_000 {
        let len = 12 + rng.below(40) as usize;
        let mut q = random_seq(&mut rng, len, 1, 28);
        let mut s = random_seq(&mut rng, len, 1, 28);
        // plenty of W/W (11) and sentinels (0) in the sequences
        for v in [&mut q, &mut s] {
            for _ in 0..rng.below(6) {
                let at = rng.below(v.len() as u64) as usize;
                v[at] = if rng.below(2) == 0 { 20 } else { 0 };
            }
        }
        let w = 1 + rng.below(11) as usize;
        let d = 1 + rng.below(127) as i32;
        let qr = rng.below(q.len() as u64) as usize;
        let sr = rng.below(s.len() as u64) as usize;
        let sl = sr.saturating_sub(rng.below(8) as usize);
        check_call(&q, &s, sl, sr, qr, d, w);
    }
}

/// Parameters outside the proven range, and a code the table does not cover, go to the
/// reference.
#[test]
fn fast_path_hands_unsupported_inputs_to_the_reference() {
    let q = vec![
        0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0,
    ];
    let s = q.clone();
    for d in [i32::MIN, -5, 0, 128, 8193, i32::MAX] {
        assert_eq!(two_hit(&q, &s, 2, 6, 6, d, 3), None, "d={d}");
        assert_eq!(one_hit(&q, &s, 6, 6, d, 3), None, "d={d}");
    }
    for w in [0usize, 12, 20] {
        assert_eq!(two_hit(&q, &s, 0, 2, 2, 16, w), None, "w={w}");
        assert_eq!(one_hit(&q, &s, 2, 2, 16, w), None, "w={w}");
    }
    let mut wide = vec![0u8; 80];
    let mut other = vec![0u8; 80];
    for i in 1..79 {
        wide[i] = 1 + (i % 20) as u8;
        other[i] = 1 + ((i * 3) % 20) as u8;
    }
    wide[36] = 200;
    assert_eq!(two_hit(&wide, &other, 30, 35, 35, 16, 3), None);
    assert_eq!(one_hit(&wide, &other, 35, 35, 16, 3), None);
    // a code >= 32 only in the walk, not in the word
    let mut deep = other.clone();
    deep[20] = 77;
    assert_eq!(two_hit(&other, &deep, 28, 30, 30, 100, 3), None);
    // the early return for a word that does not fit is the reference's
    assert_eq!(two_hit(&q, &s, 0, 21, 3, 16, 3), Some(None));
    assert_eq!(one_hit(&q, &s, 3, 21, 16, 3), Some(None));
}

/// The table the fast path scores with is the matrix function the reference calls (the
/// sentinel is stored as -128, a stop in the same step).
#[test]
fn score_table_equals_the_reference_matrix_function() {
    let t = build_table().expect("BLOSUM62 fits bytes");
    for q in 0..32usize {
        for s in 0..32usize {
            let v = blosum62_score_ncbistdaa_direct(q as u8, s as u8);
            let expect = if v == BLAST_SCORE_MIN { -128 } else { v };
            assert_eq!(t[q << 5 | s] as i32, expect, "q={q} s={s}");
        }
    }
    assert!(t[CODES * CODES..].iter().all(|&b| b == 0));
}

/// Every matrix's table holds the scores of the function the reference calls for that matrix
/// (the sentinel stored as -128), and its longest word is `127 / M`; BLOSUM62's is the table
/// and the longest word of [`two_hit`].
#[test]
fn matrix_tables_hold_the_reference_scores() {
    for matrix in ScoringMatrix::ALL {
        let m = matrix_table(matrix).expect("every matrix fits bytes");
        let mut max = 1;
        for q in 0..32usize {
            for s in 0..32usize {
                let v = protein_score(matrix, q as u8, s as u8);
                let expect = if v == BLAST_SCORE_MIN {
                    -128
                } else {
                    max = max.max(v);
                    v
                };
                assert_eq!(m.table[q << 5 | s] as i32, expect, "{matrix:?} q={q} s={s}");
            }
        }
        assert!(m.table[CODES * CODES..].iter().all(|&b| b == 0));
        assert_eq!(
            m.max_word_size,
            (WORD_SCORE_ROOM / max) as usize,
            "{matrix:?}"
        );
        assert!(
            m.max_word_size >= 7,
            "{matrix:?}: words of up to 7 residues (tblastn)"
        );
    }
    let b62 = matrix_table(ScoringMatrix::Blosum62).unwrap();
    assert_eq!(&b62.table[..], &build_table().unwrap()[..]);
    assert_eq!(b62.max_word_size, MAX_WORD_SIZE);
}

/// `two_hit_of` with each matrix's table gives the reference's result for that matrix: random
/// and related sequences with sentinels and the matrix's best pair inside them, words of every
/// size the table takes with a sentinel at each place among best pairs, and a longer word goes
/// to the reference.
#[test]
fn matrix_walks_equal_the_reference_for_every_matrix() {
    let mut rng = Rng(0x5851_F42D_4C95_7F2D);
    for matrix in ScoringMatrix::ALL {
        let m = matrix_table(matrix).unwrap();
        let (bq, bs) = (1..28u8)
            .flat_map(|q| (1..28u8).map(move |s| (q, s)))
            .max_by_key(|&(q, s)| protein_score(matrix, q, s))
            .unwrap();
        let mut compared = 0u64;
        for round in 0..300 {
            let len = 12 + rng.below(if round % 2 == 0 { 40 } else { 400 }) as usize;
            // every fifth round also codes 28..=31, which the table holds (`CODES`)
            let hi = if round % 5 == 4 { 32 } else { 28 };
            let mut q = random_seq(&mut rng, len, 1, hi);
            let mut s = if round % 3 == 0 {
                random_seq(&mut rng, len + 3, 1, hi)
            } else {
                let rate = 2 + rng.below(10);
                mutate(&mut rng, &q, rate)
            };
            for _ in 0..rng.below(10) {
                let at = rng.below(q.len().min(s.len()) as u64) as usize;
                match rng.below(3) {
                    0 => q[at] = 0,
                    1 => (q[at], s[at]) = (bq, bs),
                    _ => s[at] = 0,
                }
            }
            for _ in 0..40 {
                let w = 1 + rng.below(m.max_word_size as u64) as usize;
                let d = 1 + rng.below(127) as i32;
                let qr = rng.below(q.len() as u64) as usize;
                let sr = rng.below(s.len() as u64) as usize;
                let sl = sr.saturating_sub(rng.below(8) as usize);
                let expect = extend_two_hit(matrix, &q, &s, sl, sr, qr, d, w);
                let got =
                    two_hit_of(m, &q, &s, sl, sr, qr, d, w).expect("fast path handles the call");
                assert_eq!(
                    got, expect,
                    "{matrix:?} sl={sl} sr={sr} qr={qr} d={d} w={w}"
                );
                compared += 1;
            }
        }
        // Words of best pairs with one sentinel at each place, at the longest word too.
        for w in 1..=m.max_word_size {
            for k in 0..w {
                let n = 2 * w + 12;
                let mut q = vec![bq; n];
                let s = vec![bs; n];
                q[0] = 0;
                q[n - 1] = 0;
                let p = 5;
                q[p + k] = 0;
                for d in [1, 30, 127] {
                    let expect = extend_two_hit(matrix, &q, &s, p - 3, p, p, d, w);
                    let got = two_hit_of(m, &q, &s, p - 3, p, p, d, w)
                        .expect("fast path handles the call");
                    assert_eq!(got, expect, "{matrix:?} w={w} k={k} d={d}");
                    compared += 1;
                }
            }
        }
        assert!(compared >= 12_000, "{matrix:?}: {compared}");
        let q = random_seq(&mut rng, 40, 1, 28);
        assert_eq!(
            two_hit_of(m, &q, &q, 2, 6, 6, 16, m.max_word_size + 1),
            None
        );
        assert_eq!(two_hit_of(m, &q, &q, 2, 6, 6, 128, 3), None);
        // the other inputs `fast_path_hands_unsupported_inputs_to_the_reference` gives
        // `two_hit`: X-drops outside 1..=127, no word, and a code >= 32 in the word or only
        // in the walk; the early return of a word that does not fit is the reference's
        for d in [i32::MIN, -5, 0, 8193, i32::MAX] {
            assert_eq!(
                two_hit_of(m, &q, &q, 2, 6, 6, d, 3),
                None,
                "{matrix:?} d={d}"
            );
        }
        assert_eq!(two_hit_of(m, &q, &q, 0, 2, 2, 16, 0), None, "{matrix:?}");
        let mut wide = vec![0u8; 80];
        let mut other = vec![0u8; 80];
        for i in 1..79 {
            wide[i] = 1 + (i % 20) as u8;
            other[i] = 1 + ((i * 3) % 20) as u8;
        }
        wide[36] = 200;
        assert_eq!(
            two_hit_of(m, &wide, &other, 30, 35, 35, 16, 3),
            None,
            "{matrix:?}"
        );
        let mut deep = other.clone();
        deep[20] = 77;
        assert_eq!(
            two_hit_of(m, &other, &deep, 28, 30, 30, 100, 3),
            None,
            "{matrix:?}"
        );
        let short: Vec<u8> = [0].into_iter().chain(1..=20).chain([0]).collect();
        assert_eq!(
            two_hit_of(m, &short, &short, 0, 21, 3, 16, 3),
            Some(extend_two_hit(matrix, &short, &short, 0, 21, 3, 16, 3)),
            "{matrix:?}"
        );
        assert_eq!(
            two_hit_of(m, &short, &short, 0, 21, 3, 16, 3),
            Some(None),
            "{matrix:?}"
        );
    }
}
