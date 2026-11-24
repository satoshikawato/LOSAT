use super::super::x_ungapped_fast::build_table;
use super::*;
use crate::algorithm::blastp::extension::extend_two_hit;
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

/// A buffer of `len` residues from `lo..hi` with the sentinel 0 at both ends.
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

/// Another matrix's byte table for the walks (`MatrixTable::sent_table`, item SC2c).
fn table_of(matrix: ScoringMatrix) -> &'static Table {
    super::super::x_ungapped_fast::matrix_table(matrix)
        .and_then(|table| table.sent_table())
        .expect("the matrix's table is on the path")
}

/// The keys the walks take for `matrix`: BLOSUM62's, or another matrix's (item SC2c).
fn keys_for(matrix: ScoringMatrix) -> &'static KeyTable {
    if matrix == ScoringMatrix::Blosum62 {
        state().keys.as_deref().expect("BLOSUM62 fits bytes")
    } else {
        keys_of(matrix).expect("the matrix's table is on the path")
    }
}

/// The path on its own (independent of the switch) with the table of `matrix`: the query
/// block as blastp builds it, without its leading sentinel.
fn query_of(matrix: ScoringMatrix, block: &[u8], x_drops: &[i32]) -> Option<XSentQuery> {
    if x_drops.iter().any(|d| !(1..=MAX_DROPOFF).contains(d))
        || block.iter().fold(0u8, |a, &c| a | c) >= CODES
    {
        return None;
    }
    Some(XSentQuery {
        padded: padded_rows(block),
        keys: keys_for(matrix),
    })
}

/// One call against the reference (the matrix-generic entry point, which no switch
/// touches). `block` is the query block (blastp's `concatenated_query`, which has no leading
/// sentinel), `subject` a subject buffer with both sentinels.
fn check_call(block: &[u8], subject: &[u8], sl: usize, sr: usize, qr: usize, d: i32) {
    check_call_of(ScoringMatrix::Blosum62, block, subject, sl, sr, qr, d);
}

/// [`check_call`] with the table of `matrix` against the reference of `matrix`.
fn check_call_of(
    matrix: ScoringMatrix,
    block: &[u8],
    subject: &[u8],
    sl: usize,
    sr: usize,
    qr: usize,
    d: i32,
) {
    let q = query_of(matrix, block, &[d]).expect("block on the path");
    let s = XSentSubject::new(subject).expect("subject on the path");
    let raw = &subject[1..subject.len() - 1];
    let expect = extend_two_hit(matrix, block, raw, sl as i32, sr, qr, d, 3);
    // SAFETY: d is in 1..=127 (checked by `query_of`).
    let got = unsafe { two_hit(&q, s, sl as i32, sr, qr, d) };
    assert_eq!(
        got, expect,
        "{matrix:?} two_hit sl={sl} sr={sr} qr={qr} d={d}"
    );
}

/// Random calls on one pair; `block` and `subject` as in [`check_call`].
fn check_pair(block: &[u8], subject: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
    check_pair_of(ScoringMatrix::Blosum62, block, subject, rng, calls, dmax)
}

/// [`check_pair`] with the table of `matrix`.
fn check_pair_of(
    matrix: ScoringMatrix,
    block: &[u8],
    subject: &[u8],
    rng: &mut Rng,
    calls: usize,
    dmax: i32,
) -> u64 {
    let slen = subject.len() - 2;
    for _ in 0..calls {
        let d = 1 + rng.below(dmax as u64) as i32;
        // offsets up to the lengths, so that some words do not fit (the early return)
        let qr = rng.below(block.len() as u64 + 1) as usize;
        let sr = rng.below(slen as u64 + 1) as usize;
        let sl = if sr > 0 {
            rng.below(sr as u64 + 4) as usize
        } else {
            0
        };
        check_call_of(matrix, block, subject, sl, sr, qr, d);
    }
    calls as u64
}

/// A query block of contexts as blastp builds it: each context's residues and its trailing
/// sentinel.
fn block_of(rng: &mut Rng, contexts: usize, len: usize) -> Vec<u8> {
    let mut block = Vec::new();
    for _ in 0..contexts {
        let n = 1 + rng.below(len as u64) as usize;
        let seq = random_seq(rng, n, 1, 28);
        block.extend_from_slice(&seq[1..]);
    }
    block
}

/// Random sequences, related sequences (long walks), all codes, short and long, blocks of
/// several contexts.
#[test]
fn sentinel_extension_equals_the_reference_on_random_inputs() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut compared = 0u64;
    for round in 0..400 {
        let len = match round % 4 {
            0 => 1 + rng.below(40) as usize,
            1 => 20 + rng.below(120) as usize,
            2 => 100 + rng.below(600) as usize,
            _ => 600 + rng.below(3000) as usize,
        };
        let dmax = [4, 16, 40, 127][round / 4 % 4];
        // residues 1..=27, so every code except the sentinel appears
        let q = random_seq(&mut rng, len, 1, 28);
        let rate = 2 + rng.below(12);
        let s_related = mutate(&mut rng, &q, rate);
        let extra = rng.below(30) as usize;
        let s_random = random_seq(&mut rng, len + extra, 1, 28);
        // one context: the block is the residues and the trailing sentinel
        compared += check_pair(&q[1..], &s_related, &mut rng, 60, dmax);
        compared += check_pair(&q[1..], &s_random, &mut rng, 60, dmax);
        // a block without the trailing sentinel (the right bound is the block's end)
        compared += check_pair(&q[1..q.len() - 1], &s_related, &mut rng, 30, dmax);
        // several contexts, related to the subject in the first
        let mut block = q[1..].to_vec();
        block.extend(block_of(&mut rng, 1 + round % 3, len));
        compared += check_pair(&block, &s_related, &mut rng, 60, dmax);
        // a sentinel inside the subject
        let mut inner = s_related.clone();
        let at = 1 + rng.below(len as u64) as usize;
        inner[at] = 0;
        compared += check_pair(&q[1..], &inner, &mut rng, 40, dmax);
    }
    assert!(compared >= 100_000);
}

/// X-drop at the boundary and ties for the best position: walks whose score returns to the
/// best value exactly, drops that equal the dropoff, repeated patterns.
#[test]
fn sentinel_extension_equals_the_reference_at_boundaries_and_ties() {
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
            compared += check_pair(&q[1..], &s, &mut rng, 30, 25);
        }
    }
    // Every dropoff 1..=127 on a periodic pattern: +11 (W) among +4 (A-A) steps.
    let w_code = 20u8;
    let q: Vec<u8> = (0..300)
        .map(|i| if i % 7 == 3 { w_code } else { 1 })
        .chain([0u8])
        .collect();
    let s: Vec<u8> = [0u8]
        .into_iter()
        .chain((0..300).map(|i| if i % 5 == 1 { w_code } else { 1 }))
        .chain([0u8])
        .collect();
    for d in 1..=MAX_DROPOFF {
        for sr in 0..q.len() - 3 {
            check_call(&q, &s, sr.saturating_sub(5), sr, sr, d);
            compared += 1;
        }
    }
    assert!(compared > 50_000);
}

/// Very long walks: identical sequences whose score passes 16 bits, and walks that start next
/// to either end of either sequence (the sentinel stops at offset -1 and at the length).
#[test]
fn sentinel_extension_equals_the_reference_on_very_long_walks() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    // 9000 tryptophans: 99,000 > 2^16 (W-W = 11)
    let mut s = vec![0u8];
    s.extend(std::iter::repeat_n(20u8, 9000));
    s.push(0);
    let q = s[1..].to_vec();
    for d in [1, 5, 16, 100, 127] {
        for (qr, sr) in [
            (0usize, 0usize),
            (1, 1),
            (4000, 4000),
            (8990, 8990),
            (8997, 8997),
            (5000, 20),
            (20, 5000),
            (0, 8000),
            (8000, 0),
        ] {
            for back in [0usize, 1, 10] {
                check_call(&q, &s, sr.saturating_sub(back), sr, qr, d);
            }
        }
    }
    // a long walk with mismatches that never reach the dropoff
    for i in (0..s.len()).step_by(9) {
        if s[i] != 0 {
            s[i] = 1;
        }
    }
    for _ in 0..200 {
        let qr = rng.below(8998) as usize;
        let sr = rng.below(8998) as usize;
        let sl = sr.saturating_sub(rng.below(30) as usize);
        let d = 20 + rng.below(100) as i32;
        check_call(&q, &s, sl, sr, qr, d);
    }
}

/// Words with sentinels and the extreme scores inside them, and words next to either end of a
/// sequence.
#[test]
fn sentinel_extension_equals_the_reference_for_sentinels_in_words() {
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    for _ in 0..20_000 {
        let len = 4 + rng.below(40) as usize;
        let q = random_seq(&mut rng, len, 1, 28);
        let mut block = q[1..].to_vec();
        let mut s = random_seq(&mut rng, len, 1, 28);
        // plenty of W/W (11) and sentinels (0) inside the sequences
        for _ in 0..rng.below(6) {
            let at = rng.below(block.len() as u64) as usize;
            block[at] = if rng.below(2) == 0 { 20 } else { 0 };
        }
        for _ in 0..rng.below(6) {
            let at = 1 + rng.below(s.len() as u64 - 2) as usize;
            s[at] = if rng.below(2) == 0 { 20 } else { 0 };
        }
        let d = 1 + rng.below(127) as i32;
        let qr = rng.below(block.len() as u64 + 1) as usize;
        let sr = rng.below(s.len() as u64 - 1) as usize;
        let sl = sr.saturating_sub(rng.below(8) as usize);
        check_call(&block, &s, sl, sr, qr, d);
    }
}

/// The keys of the module comment: a best key decodes to its score and length, a tie keeps the
/// first step, and the keys of a matrix are `(t << 32) - 1`.
#[test]
fn keys_decode_and_order_as_the_reference_compares() {
    // no best yet: the length is 0
    for maxscore in [0i32, 1, 77, 1 << 20] {
        assert_eq!(best_of((maxscore as i64) << 32), (maxscore, 0));
    }
    // a best at step i: (score << 32) - (i + 1)
    for (score, i) in [
        (0i32, 0i64),
        (5, 0),
        (-3, 7),
        (100_000, 9000),
        (11, (1 << 31) - 2),
    ] {
        let key = ((score as i64) << 32) - (i + 1);
        assert_eq!(
            best_of(key),
            (score, (i + 1) as u32),
            "score {score} step {i}"
        );
    }
    // the same score later is not above the first (the reference's `>`), a higher one is
    let first = (40i64 << 32) - 3;
    assert!(((40i64 << 32) - 9) < first);
    assert!(((41i64 << 32) - 9_000_000) > first);
    // the X-drop test: best - key >= dropoff << 32 exactly when best - score >= dropoff
    for (best, bi, score, si, d) in [
        (30i64, 2i64, 14i64, 9i64, 16i64),
        (30, 2, 15, 9, 16),
        (30, 2, 15, 3, 15),
    ] {
        let bk = (best << 32) - (bi + 1);
        let k = (score << 32) - (si + 1);
        assert_eq!(k + (d << 32) <= bk, best - score >= d, "{best} {score} {d}");
    }
    let table = build_table().expect("BLOSUM62 fits bytes");
    let keys = key_table(&table);
    assert_eq!(blosum62_keys().as_deref(), Some(&*keys));
    for q in 0..CODES as usize {
        for c in 0..CODES as usize {
            assert_eq!(keys[q << 5 | c], ((table[q << 5 | c] as i64) << 32) - 1);
        }
    }
    assert_eq!(
        padded_rows(&[1, 20, 31]),
        vec![0, 1 << 5, 20 << 5, 31 << 5, 0]
    );
}

/// The checks keep everything the sentinel stop is not exact for on the reference closure.
#[test]
fn checks_refuse_inputs_outside_the_path() {
    let block: Vec<u8> = (1..=20).chain([0]).collect();
    let x_drops = [16, 127, 1];
    if state().on {
        assert!(XSentQuery::new(&block, &x_drops, 3, 3).is_some());
    }
    // X-drops outside 1..=127, other word sizes, a code >= 32: refused whatever the switch says
    for d in [i32::MIN, -5, 0, 128, 8193, i32::MAX] {
        assert!(XSentQuery::new(&block, &[16, d], 3, 3).is_none(), "d={d}");
    }
    for (lw, ew) in [(2, 2), (2, 3), (3, 2), (5, 5), (6, 6), (3, 11)] {
        assert!(
            XSentQuery::new(&block, &x_drops, lw, ew).is_none(),
            "w={lw}/{ew}"
        );
    }
    let mut wide = block.clone();
    wide[7] = 32;
    assert!(XSentQuery::new(&wide, &x_drops, 3, 3).is_none());
    // subjects: a sentinel at each end, codes below 32, two bytes at least
    let s: Vec<u8> = [0u8].into_iter().chain(1..=20).chain([0]).collect();
    assert!(XSentSubject::new(&s).is_some());
    assert!(XSentSubject::new(&[0, 0]).is_some());
    assert!(XSentSubject::new(&[0]).is_none());
    assert!(XSentSubject::new(&[]).is_none());
    assert!(XSentSubject::new(&s[1..]).is_none());
    assert!(XSentSubject::new(&s[..s.len() - 1]).is_none());
    let mut deep = s.clone();
    deep[9] = 77;
    assert!(XSentSubject::new(&deep).is_none());
    deep[9] = 31;
    assert!(XSentSubject::new(&deep).is_some());
    // the early return for a word that does not fit is the reference's
    check_call(&block, &s, 0, 18, 3, 16);
    check_call(&block, &s, 0, 3, 19, 16);
    check_call(&block, &s, 0, 20, 21, 16);
}

/// The matrices other than BLOSUM62 that BLASTP reads (item SC2c).
const OTHER_MATRICES: [ScoringMatrix; 7] = [
    ScoringMatrix::Blosum45,
    ScoringMatrix::Blosum50,
    ScoringMatrix::Blosum80,
    ScoringMatrix::Blosum90,
    ScoringMatrix::Pam30,
    ScoringMatrix::Pam70,
    ScoringMatrix::Pam250,
];

/// Every other matrix that BLASTP reads takes the sentinel walks: its table gives the sentinel
/// -128 both ways and its largest score is at most 31, and `new_of` takes a block on the path
/// (with the switch on) and refuses the inputs that `new` refuses.
#[test]
fn other_matrices_take_the_sentinel_walks() {
    let block: Vec<u8> = (1..=20).chain([0]).collect();
    let x_drops = [16, 127, 1];
    for matrix in OTHER_MATRICES {
        let table = table_of(matrix);
        for c in 0..CODES as usize {
            assert_eq!(table[c], -128, "{matrix:?} sentinel against {c}");
            assert_eq!(table[c << 5], -128, "{matrix:?} {c} against the sentinel");
        }
        if state().on {
            assert!(XSentQuery::new_of(matrix, &block, &x_drops, 3, 3).is_some());
        }
        for d in [0, 128] {
            assert!(XSentQuery::new_of(matrix, &block, &[16, d], 3, 3).is_none());
        }
        for (lw, ew) in [(2, 3), (3, 2), (5, 5), (6, 6)] {
            assert!(XSentQuery::new_of(matrix, &block, &x_drops, lw, ew).is_none());
        }
        let mut wide = block.clone();
        wide[7] = 32;
        assert!(XSentQuery::new_of(matrix, &wide, &x_drops, 3, 3).is_none());
    }
    // The length bound (`MAX_BLOCK_LEN_OF` in `new_of`), here on a small bound.
    assert!(fits_of(&block, &x_drops, 3, 3, block.len()));
    assert!(!fits_of(&block, &x_drops, 3, 3, block.len() - 1));
}

/// The other matrices' tables on random, related and boundary inputs, as for BLOSUM62.
#[test]
fn sentinel_extension_of_other_matrices_equals_the_reference() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut compared = 0u64;
    for matrix in OTHER_MATRICES {
        for round in 0..80 {
            let len = match round % 4 {
                0 => 1 + rng.below(40) as usize,
                1 => 20 + rng.below(120) as usize,
                2 => 100 + rng.below(600) as usize,
                _ => 600 + rng.below(2000) as usize,
            };
            let dmax = [4, 16, 40, 127][round / 4 % 4];
            let q = random_seq(&mut rng, len, 1, 28);
            let rate = 2 + rng.below(12);
            let s_related = mutate(&mut rng, &q, rate);
            let extra = rng.below(30) as usize;
            let s_random = random_seq(&mut rng, len + extra, 1, 28);
            compared += check_pair_of(matrix, &q[1..], &s_related, &mut rng, 40, dmax);
            compared += check_pair_of(matrix, &q[1..], &s_random, &mut rng, 40, dmax);
            compared += check_pair_of(matrix, &q[1..q.len() - 1], &s_related, &mut rng, 20, dmax);
            let mut block = q[1..].to_vec();
            block.extend(block_of(&mut rng, 1 + round % 3, len));
            compared += check_pair_of(matrix, &block, &s_related, &mut rng, 40, dmax);
            let mut inner = s_related.clone();
            inner[1 + rng.below(len as u64) as usize] = 0;
            compared += check_pair_of(matrix, &q[1..], &inner, &mut rng, 20, dmax);
        }
    }
    assert!(compared >= 7 * 80 * 160);
}
