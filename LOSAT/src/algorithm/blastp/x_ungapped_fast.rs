//! EXPERIMENT (LOSAT_X_BPUNGAP / LOSAT_X_BPUNGAPSHADOW, off unless set): the BLOSUM62 protein
//! ungapped extension (`s_BlastAaExtendOneHit`, `s_BlastAaExtendTwoHit` and the
//! `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` walks they share) with a tighter inner loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883
//! ```c
//! for (i = 0; i < n; i++) {
//!     score += matrix[q[i]][s[i]];
//!     if (score > maxscore) { maxscore = score; best_i = i; }
//!     if (score <= 0 || (maxscore - score) >= dropoff) break;
//! }
//! *length = best_i + 1;
//! *s_last_off = s_off + i;
//! ```
//! (right walk; the left walk counts `i` down from `n`, has no `score <= 0` test, and returns
//! `n - best_i + 1`.)
//!
//! The reference walks cost about 22 instructions a step: two bounds checks, a range test and a
//! branch for each of the two residue codes, a multiply for the matrix row, and the update. A
//! step here is a load of each residue, `q << 5 | s`, one signed-byte load from a flat table, the
//! add, and the two compare-and-select updates of NCBI's loop; the statements, their order and
//! the exit tests are those of NCBI. (A 16-lane AVX2 version of the walk, with the scores
//! gathered from the same table and the stop found by prefix sums and prefix maxima, was
//! measured too and was not faster: the two `vpgatherdd` of a block and the dependent scans cost
//! as much as 16 scalar steps, whose exit branch is predictable.)
//!
//! Why the result is the same.
//! * Bounds. The reference indexes with `[]`. Here each walk reads `n` positions, with `n`
//!   derived from the slice lengths exactly as the reference derives it (right: the room left in
//!   both slices; left: `min(q_off, s_off) + 1`, after an `assert!` that the start is inside),
//!   so every position is inside both slices.
//! * Table. `T[q << 5 | s]` is `blosum62_score_ncbistdaa_direct(q, s)` (the function the
//!   reference calls) for `q, s < 32`, with the one change that `BLAST_SCORE_MIN` (-32768, the
//!   score of the sentinel code 0 against anything) is stored as -128 so that an entry is one
//!   byte; the table is used only if every other entry is within -127..=127. A residue code
//!   >= 32 would index a part of the 8 KiB table that holds zeros, so each extension ORs
//!   together every code it reads and is handed to the reference when the OR is >= 32.
//! * Sentinel. Take `1 <= dropoff <= 127` (otherwise the reference runs). The score before a
//!   step is never above the best (`best` is raised whenever `score` exceeds it, and both start
//!   equal), so a sentinel step makes `best - score >= 128 > dropoff` with the stored -128 as it
//!   does with -32768: the walk stops at that step in both, the step does not raise the best in
//!   either, and the rest of the result (`best`, `best_i`, the position of the stop) is the
//!   same. The word loops add scores without a stop: for a word of at most 11 residues a
//!   sentinel (-128) is not made positive again by the others (<= 11 each), as -32768 is not,
//!   so every comparison with a score >= 0 has the same outcome; longer words run the
//!   reference.

use super::extension::{
    BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult, BlastpUngappedData,
};
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, BLAST_SCORE_MIN};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// Codes below this index the table; `q << 5 | s` is the index.
const CODES: usize = 32;
/// `q` and `s` are bytes, so every `q << 5 | s` is below this: a code >= 32 reads zeros.
const TABLE_LEN: usize = 255 * CODES + 256;
/// Largest X-drop the byte table is exact for (see the module comment).
const MAX_DROPOFF: i32 = 127;
/// Longest word whose sentinel sum stays below zero (11 * 11 < 128).
const MAX_WORD_SIZE: usize = 11;

type Table = [i8; TABLE_LEN];

struct State {
    /// 0 = off, 1 = LOSAT_X_BPUNGAP, 2 = LOSAT_X_BPUNGAPSHADOW (both paths, compared).
    mode: u8,
    table: Box<Table>,
}

/// The table of the walks, `None` if a score does not fit a byte.
// NCBI reference: aa_ungapped.c:804-883 `matrix[q[i]][s[i]]`; the scores are those of the
// function the reference walks call.
fn build_table() -> Option<Box<Table>> {
    let mut t = Box::new([0i8; TABLE_LEN]);
    for q in 0..CODES {
        for s in 0..CODES {
            let v = blosum62_score_ncbistdaa_direct(q as u8, s as u8);
            t[q << 5 | s] = if v == BLAST_SCORE_MIN {
                -128
            } else if (-127..=127).contains(&v) {
                v as i8
            } else {
                return None;
            };
        }
    }
    Some(t)
}

fn state() -> &'static State {
    static STATE: OnceLock<State> = OnceLock::new();
    STATE.get_or_init(|| {
        // Not in `xswitch::STRICT_SET`, so not read with `xswitch::on`.
        let requested = if std::env::var_os("LOSAT_X_BPUNGAPSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_BPUNGAP").is_some() {
            1
        } else {
            0
        };
        // The table is built in every mode (the tests use it); a score that does not fit a byte
        // keeps the switch off.
        match build_table() {
            Some(table) => State {
                mode: requested,
                table,
            },
            None => State {
                mode: 0,
                table: Box::new([0i8; TABLE_LEN]),
            },
        }
    })
}

/// 0 = off, 1 = LOSAT_X_BPUNGAP, 2 = LOSAT_X_BPUNGAPSHADOW.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline]
pub(crate) fn mode() -> u8 {
    state().mode
}

static SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// Counts one compared extension (shadow mode).
// No NCBI counterpart: a counter.
#[inline]
pub(crate) fn count_shadow() {
    SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
}

/// The summary line (printed by `main` at exit), only in shadow mode.
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] BPUNGAP compared={} (ungapped extensions identical)",
            SHADOW_COMPARED.load(Ordering::Relaxed)
        );
    }
}

/// What one walk returns: best score, length, step at which the loop left, and the OR of every
/// residue code it read.
struct Walk {
    best: i32,
    len: i32,
    end: usize,
    codes: u8,
}

/// One walk of `s_BlastAaExtendRight` (`RIGHT`) or `s_BlastAaExtendLeft` over `n` steps: step `i`
/// reads `qp[i]`, `sp[i]` (right) or `qp[-i]`, `sp[-i]` (left).
///
/// SAFETY: the `n` positions of both sequences are readable.
// NCBI reference: aa_ungapped.c:804-883 (see the module comment).
#[inline(always)]
unsafe fn walk<const RIGHT: bool>(
    tab: &Table,
    qp: *const u8,
    sp: *const u8,
    n: usize,
    dropoff: i32,
    maxscore: i32,
) -> Walk {
    let mut score = maxscore;
    let mut best = maxscore;
    let mut best_len = 0i32;
    let mut codes = 0u8;
    let mut i = 0usize;
    while i < n {
        let (qc, sc) = if RIGHT {
            (*qp.add(i), *sp.add(i))
        } else {
            (*qp.sub(i), *sp.sub(i))
        };
        codes |= qc | sc;
        // SAFETY: q << 5 | s < 255 * 32 + 256 = TABLE_LEN for any two bytes.
        score += *tab.get_unchecked((qc as usize) << 5 | sc as usize) as i32;
        if score > best {
            best = score;
            best_len = i as i32 + 1;
        }
        let stop = if RIGHT {
            score <= 0 || (best - score) >= dropoff
        } else {
            (best - score) >= dropoff
        };
        if stop {
            break;
        }
        i += 1;
    }
    Walk {
        best,
        len: best_len,
        end: i,
        codes,
    }
}

/// `s_BlastAaExtendLeft(.., s_off, q_off, ..)`: (score, length, codes read).
#[inline(always)]
fn extend_left(
    tab: &Table,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, u8) {
    let n = s_off.min(q_off);
    if n < 0 {
        // NCBI reference: aa_ungapped.c:853-883, `best_i = n + 1` with no step taken.
        return (maxscore, 0, 0);
    }
    // n + 1 steps down from (q_off, s_off); both offsets are >= n.
    assert!((q_off as usize) < query.len() && (s_off as usize) < subject.len());
    // SAFETY: the assert and n <= min(q_off, s_off) keep the n + 1 positions inside.
    let w = unsafe {
        walk::<false>(
            tab,
            query.as_ptr().add(q_off as usize),
            subject.as_ptr().add(s_off as usize),
            n as usize + 1,
            dropoff,
            maxscore,
        )
    };
    (w.best, w.len, w.codes)
}

/// `s_BlastAaExtendRight(.., s_off, q_off, ..)`: (score, length, s_last_off, codes read).
#[inline(always)]
fn extend_right(
    tab: &Table,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, i32, u8) {
    let n = (subject.len() as i32 - s_off).min(query.len() as i32 - q_off);
    if n <= 0 {
        // NCBI reference: aa_ungapped.c:804-835, the loop does not run: i = 0.
        return (maxscore, 0, s_off, 0);
    }
    // n steps up from (q_off, s_off); n is the room left in both slices. A negative offset was
    // an index panic in the reference.
    assert!(q_off >= 0 && s_off >= 0);
    // SAFETY: the n positions are inside both slices.
    let w = unsafe {
        walk::<true>(
            tab,
            query.as_ptr().add(q_off as usize),
            subject.as_ptr().add(s_off as usize),
            n as usize,
            dropoff,
            maxscore,
        )
    };
    (w.best, w.len, s_off + w.end as i32, w.codes)
}

/// Whether the fast path covers these parameters.
#[inline(always)]
fn supported(query: &[u8], subject: &[u8], dropoff: i32, word_size: usize) -> bool {
    (1..=MAX_DROPOFF).contains(&dropoff)
        && (1..=MAX_WORD_SIZE).contains(&word_size)
        && query.len() <= i32::MAX as usize
        && subject.len() <= i32::MAX as usize
}

/// `Some(result)` when the fast path handled the call (`result` may be `None`, the reference's
/// early return), `None` when the reference has to run.
// NCBI reference: aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
pub(crate) fn two_hit(
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpTwoHitUngappedResult>> {
    // The early return of the reference wrapper (a word that does not fit).
    if q_right_off + word_size > query.len() || s_right_off + word_size > subject.len() {
        return Some(None);
    }
    if !supported(query, subject, dropoff, word_size) {
        return None;
    }
    let tab: &Table = &state().table;
    let mut q_right_off = q_right_off as i32;
    let s_left_off = s_left_off as i32;
    let mut s_right_off = s_right_off as i32;
    let word_size = word_size as i32;

    let mut right_d = 0i32;
    let mut left_score = 0i32;
    let mut score = 0i32;
    let mut codes = 0u8;
    for i in 0..word_size {
        // the word fits in both sequences (checked above)
        let qc = query[(q_right_off + i) as usize];
        let sc = subject[(s_right_off + i) as usize];
        codes |= qc | sc;
        score += tab[(qc as usize) << 5 | sc as usize] as i32;
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;

    let (left_score, left_d, c) = extend_left(
        tab,
        query,
        subject,
        q_right_off - 1,
        s_right_off - 1,
        dropoff,
        0,
    );
    codes |= c;

    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    let mut right_score = 0i32;
    let mut right_d = 0i32;
    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        let (rs, rd, last, c) = extend_right(
            tab,
            query,
            subject,
            q_right_off,
            s_right_off,
            dropoff,
            left_score,
        );
        (right_score, right_d, s_last_off) = (rs, rd, last);
        codes |= c;
    }
    if codes >= CODES as u8 {
        return None;
    }

    Some(Some(BlastpTwoHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_d,
            score: left_score.max(right_score),
        },
        right_extend,
        s_last_off,
    }))
}

/// As [`two_hit`], for `s_BlastAaExtendOneHit`.
// NCBI reference: aa_ungapped.c:1020-1083 `s_BlastAaExtendOneHit`.
pub(crate) fn one_hit(
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpOneHitUngappedResult>> {
    if q_off + word_size > query.len() || s_off + word_size > subject.len() {
        return Some(None);
    }
    if !supported(query, subject, dropoff, word_size) {
        return None;
    }
    let tab: &Table = &state().table;
    let q_off = q_off as i32;
    let s_off = s_off as i32;
    let word_size = word_size as i32;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut q_left_off = q_off;
    let mut q_right_off = q_off + word_size;
    let mut q_best_left_off = q_off;
    let mut codes = 0u8;
    for i in 0..word_size {
        // the word fits in both sequences (checked above)
        let qc = query[(q_off + i) as usize];
        let sc = subject[(s_off + i) as usize];
        codes |= qc | sc;
        sum += tab[(qc as usize) << 5 | sc as usize] as i32;
        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = q_off + i;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = q_off + i + 1;
        }
    }

    let init_hit_width = q_right_off - q_left_off + 1;
    q_left_off = q_best_left_off;
    let s_left_off = q_left_off + (s_off - q_off);
    let s_right_off = q_right_off + (s_off - q_off);

    let (left_score, left_disp, c) = extend_left(
        tab,
        query,
        subject,
        q_left_off - 1,
        s_left_off - 1,
        dropoff,
        score,
    );
    codes |= c;
    let (total_score, right_disp, s_last_off, c) = extend_right(
        tab,
        query,
        subject,
        q_right_off + 1,
        s_right_off + 1,
        dropoff,
        left_score,
    );
    codes |= c;
    if codes >= CODES as u8 {
        return None;
    }

    Some(Some(BlastpOneHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_left_off - left_disp,
            s_start: s_left_off - left_disp,
            length: left_disp + right_disp + init_hit_width,
            score: total_score,
        },
        s_last_off,
    }))
}

#[cfg(test)]
mod tests {
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
}
