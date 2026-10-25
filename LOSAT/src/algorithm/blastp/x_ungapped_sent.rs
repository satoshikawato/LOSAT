//! EXPERIMENT (LOSAT_X_BPUNGAPSENT): blastp's two-hit ungapped extension (`s_BlastAaExtendTwoHit`
//! with the `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` walks) for word size 3, inlined into the
//! candidate closure of the fused seed loop (`blast_engine.rs`), with walks that the sentinels stop
//! instead of a bound.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 (the walks) and
//! 1089-1155 (`s_BlastAaExtendTwoHit`), as in `x_ungapped_fast.rs` (LOSAT_X_BPUNGAP), whose byte
//! table this path scores with and whose statements it keeps in the same order. What changes:
//! * The checks the fast path makes on every call are made once: residue codes below 32 and the
//!   X-drops of every context in `1..=127` for a query block ([`XSentQuery::new`]), codes below 32
//!   and a sentinel at each end for a subject ([`XSentSubject::new`]), word size 3. A block or a
//!   subject that fails a check runs the reference closure.
//! * The walks have no bound test. The reference walks `n` steps at most, where `n` is the room
//!   left in both slices (right) or `min(q_off, s_off) + 1` (left); the position `n` steps away
//!   is just past the end of the slice that sets `n`, and that position holds a sentinel (code 0):
//!   the subject buffer carries one at each end (`aa_seq`), and the query block is copied with one
//!   added at its start and its end. A step on a sentinel scores -128 in the byte table; the score
//!   before a step is never above the best, so the step leaves the best as it is and makes
//!   `best - score >= 128 > dropoff`: the walk stops at step `n` with the best, its length and the
//!   stop position of the reference, which left its loop at `i = n`. Steps before `n` read the
//!   same bytes as the reference (a sentinel among them stops both, as in `x_ungapped_fast.rs`).
//! * The word loop is unrolled for three letters (the same statements and order) and the walks
//!   take two steps per loop.
//! * No score wraps (a wrapped score could step over a sentinel without stopping): a block on the
//!   path has at most `i32::MAX / 16` letters, the word and the two walks of one extension read
//!   distinct query positions at +11 a step at most, so a score stays below `11 * (len + 3)`, and
//!   a score is never below `best - 255` (a step costs at most 128 and the walk stops once
//!   `best - score` reaches `dropoff <= 127`).
//!
//! BLASTX's seed loops call the same extension, [`two_hit_in`] with `TAIL` 1, on buffers laid out
//! differently (`blastx/x_ungap_sent.rs`, whose header gives the argument for them).
//!
//! In a replay of the closure on 4.3 million captured candidates of an OrthoFinder job (line F,
//! `v4-F/results/s1/p2`, arm V6) every result equalled the reference's, at 204.7 instead of 264.7
//! cycles per candidate.

use super::extension::{BlastpTwoHitUngappedResult, BlastpUngappedData};
use super::x_ungapped_fast::{build_table, Table};
use std::sync::OnceLock;

/// Codes below this index the table (`q << 5 | s`).
const CODES: u8 = 32;
/// The only word size this path takes.
const WORD_SIZE: usize = 3;
/// Largest X-drop the sentinel stop is exact for (see the module comment).
const MAX_DROPOFF: i32 = 127;
/// Longest query block on the path: no score of a walk wraps (see the module comment).
const MAX_BLOCK_LEN: usize = (i32::MAX / 16) as usize;
const _: () = assert!(11 * (MAX_BLOCK_LEN as i64 + 3) < i32::MAX as i64);

struct State {
    /// `LOSAT_X_BPUNGAPSENT`.
    on: bool,
    table: Option<Box<Table>>,
}

fn state() -> &'static State {
    static STATE: OnceLock<State> = OnceLock::new();
    STATE.get_or_init(|| State {
        // `LOSAT_X_BPUNGAPSENT` is a strict-set switch.
        on: crate::utils::xswitch::on("LOSAT_X_BPUNGAPSENT"),
        table: build_table(),
    })
}

/// Whether the switch is on.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn on() -> bool {
    state().on
}

/// A query block made ready for [`two_hit`]: a copy with a sentinel added at each end.
pub(crate) struct XSentQuery {
    /// `[0, query.., 0]`.
    padded: Vec<u8>,
    table: &'static Table,
}

impl XSentQuery {
    /// `None` when the switch is off or the block is outside the path: a code >= 32, an X-drop
    /// outside `1..=127`, a word size other than 3 (`lookup_word_size` is the lookup table's,
    /// `extension_word_size` the one the extension is called with), or more than `i32::MAX / 16`
    /// letters (no score wraps; module comment).
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    pub(crate) fn new(
        query: &[u8],
        x_drops: &[i32],
        lookup_word_size: i32,
        extension_word_size: usize,
    ) -> Option<Self> {
        let state = state();
        if !state.on {
            return None;
        }
        let table = state.table.as_deref()?;
        let fits = lookup_word_size == WORD_SIZE as i32
            && extension_word_size == WORD_SIZE
            && x_drops.iter().all(|d| (1..=MAX_DROPOFF).contains(d))
            && query.len() <= MAX_BLOCK_LEN
            && query.iter().fold(0u8, |a, &c| a | c) < CODES;
        if !fits {
            return None;
        }
        let mut padded = Vec::with_capacity(query.len() + 2);
        padded.push(0);
        padded.extend_from_slice(query);
        padded.push(0);
        Some(Self { padded, table })
    }
}

/// A subject made ready for [`two_hit`]: its buffer with the sentinels, checked.
#[derive(Clone, Copy)]
pub(crate) struct XSentSubject<'a> {
    /// `[0, subject.., 0]`.
    aa_seq: &'a [u8],
}

impl<'a> XSentSubject<'a> {
    /// `None` unless `aa_seq` is a subject buffer with a sentinel at each end, no code >= 32
    /// inside and a length within `i32`.
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    #[inline]
    pub(crate) fn new(aa_seq: &'a [u8]) -> Option<Self> {
        let ok = aa_seq.len() >= 2
            && aa_seq.len() < i32::MAX as usize
            && aa_seq[0] == 0
            && aa_seq[aa_seq.len() - 1] == 0
            && aa_seq.iter().fold(0u8, |a, &c| a | c) < CODES;
        ok.then_some(Self { aa_seq })
    }

    /// The buffer, `[0, subject.., 0]`.
    #[inline]
    pub(crate) fn with_sentinels(self) -> &'a [u8] {
        self.aa_seq
    }
}

/// One walk of `s_BlastAaExtendRight` (`RIGHT`) or `s_BlastAaExtendLeft` from `qp`, `sp`, until
/// the X-drop test (or, right, `score <= 0`) holds: (best, length, steps before the stop).
///
/// SAFETY: walking from `qp`, `sp` in the direction, both sequences reach a sentinel (code 0)
/// with every position up to it readable, all codes up to it are below 32, `tab` is the byte
/// table and `1 <= dropoff <= 127`.
// NCBI reference: aa_ungapped.c:804-883 (see the module comment).
#[inline(always)]
unsafe fn walk<const RIGHT: bool>(
    tab: *const i8,
    qp: *const u8,
    sp: *const u8,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, usize) {
    let mut score = maxscore;
    let mut best = maxscore;
    // the offset of the best step; the start value means "no best yet" (length 0)
    let none: isize = if RIGHT { -1 } else { 1 };
    let mut best_off = none;
    let mut off: isize = 0;
    let dir: isize = if RIGHT { 1 } else { -1 };
    macro_rules! step {
        ($l:lifetime) => {{
            let qc = *qp.offset(off) as usize;
            let sc = *sp.offset(off) as usize;
            score += *tab.add(qc << 5 | sc) as i32;
            if score > best {
                best = score;
                best_off = off;
            }
            let stop = if RIGHT {
                score <= 0 || (best - score) >= dropoff
            } else {
                (best - score) >= dropoff
            };
            if stop {
                break $l;
            }
            off += dir;
        }};
    }
    'w: loop {
        step!('w);
        step!('w);
    }
    let len = if best_off == none {
        0
    } else {
        (best_off.unsigned_abs() + 1) as i32
    };
    (best, len, off.unsigned_abs())
}

/// `s_BlastAaExtendTwoHit` for word size 3 on a checked query block and subject (`subject_raw`
/// offsets, as `extend_two_hit_blosum62` takes them): the reference's result, `None` for its
/// early return (a word that does not fit).
///
/// SAFETY: `1 <= dropoff <= 127` (an X-drop that [`XSentQuery::new`] checked).
// NCBI reference: aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
#[inline(always)]
pub(crate) unsafe fn two_hit(
    query: &XSentQuery,
    subject: XSentSubject<'_>,
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
) -> Option<BlastpTwoHitUngappedResult> {
    two_hit_in::<2>(
        query.table,
        &query.padded,
        subject.aa_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        dropoff,
    )
}

/// [`two_hit`] on two buffers that start with a sentinel and end with one: the reference's
/// sequences are `query[1..query.len() - TAIL + 1]` and the same of `subject` (`TAIL` 2: blastp's
/// query block and `subject_raw`, without the last sentinel; 1: BLASTX's six-frame query and
/// subject, which keep it). A const per caller keeps each caller's instance its own.
///
/// SAFETY: `1 <= dropoff <= 127`; both buffers have a sentinel (code 0) at index 0 and at their
/// last index, at most `i32::MAX / 16` letters in `query`, fewer than `i32::MAX` bytes in
/// `subject`, and `TAIL` is 1 or 2.
// NCBI reference: aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
#[inline(always)]
pub(crate) unsafe fn two_hit_in<const TAIL: usize>(
    table: &Table,
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
) -> Option<BlastpTwoHitUngappedResult> {
    debug_assert!((1..=MAX_DROPOFF).contains(&dropoff));
    debug_assert!(TAIL == 1 || TAIL == 2);
    let qlen = query.len() - TAIL;
    let slen = subject.len() - TAIL;
    // The early return of the reference wrapper (a word that does not fit).
    if q_right_off + WORD_SIZE > qlen || s_right_off + WORD_SIZE > slen {
        return None;
    }
    let tab = table.as_ptr();
    // offset 0 of the reference's sequences; one sentinel before each
    let qb = query.as_ptr().add(1);
    let sb = subject.as_ptr().add(1);
    let mut q_right_off = q_right_off as i32;
    let s_left_off = s_left_off as i32;
    let mut s_right_off = s_right_off as i32;

    // the word loop of three letters (inside both sequences: checked above)
    let mut right_d = 0i32;
    let mut left_score = 0i32;
    let mut score = 0i32;
    let qw = qb.add(q_right_off as usize);
    let sw = sb.add(s_right_off as usize);
    for i in 0..WORD_SIZE {
        score += *tab.add((*qw.add(i) as usize) << 5 | *sw.add(i) as usize) as i32;
        if score > left_score {
            left_score = score;
            right_d = i as i32 + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;

    // left: from (q_right_off - 1, s_right_off - 1) down; no step when either is -1
    let (left_score, left_d) = if q_right_off == 0 || s_right_off == 0 {
        (0, 0)
    } else {
        // SAFETY: the query block and the subject have a sentinel before offset 0.
        let (best, len, _) = walk::<false>(
            tab,
            qb.offset(q_right_off as isize - 1),
            sb.offset(s_right_off as isize - 1),
            dropoff,
            0,
        );
        (best, len)
    };

    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    let mut right_score = 0i32;
    let mut right_d = 0i32;
    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        // SAFETY: with TAIL 2, q_right_off <= qlen and s_right_off <= slen, and both have a
        // sentinel at offset qlen / slen. With TAIL 1 the last sentinel is at offset qlen - 1 /
        // slen - 1, and the walk starts at or before it: the word ends at that offset at most, and
        // a word whose last letter is a sentinel has right_d < 3 (that letter scores -128, the
        // other two 11 each at most), so the start is past it only when it is not in the word.
        let (best, len, end) = walk::<true>(
            tab,
            qb.add(q_right_off as usize),
            sb.add(s_right_off as usize),
            dropoff,
            left_score,
        );
        (right_score, right_d, s_last_off) = (best, len, s_right_off + end as i32);
    }

    Some(BlastpTwoHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_d,
            score: left_score.max(right_score),
        },
        right_extend,
        s_last_off,
    })
}

#[cfg(test)]
mod tests {
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

    /// The path on its own (independent of the switch): the query block as blastp builds it,
    /// without its leading sentinel.
    fn query(block: &[u8], x_drops: &[i32]) -> Option<XSentQuery> {
        if x_drops.iter().any(|d| !(1..=MAX_DROPOFF).contains(d))
            || block.iter().fold(0u8, |a, &c| a | c) >= CODES
        {
            return None;
        }
        let mut padded = vec![0u8];
        padded.extend_from_slice(block);
        padded.push(0);
        Some(XSentQuery {
            padded,
            table: state().table.as_deref().expect("BLOSUM62 fits bytes"),
        })
    }

    /// One call against the reference (the matrix-generic entry point, which no switch
    /// touches). `block` is the query block (blastp's `concatenated_query`, which has no leading
    /// sentinel), `subject` a subject buffer with both sentinels.
    fn check_call(block: &[u8], subject: &[u8], sl: usize, sr: usize, qr: usize, d: i32) {
        let q = query(block, &[d]).expect("block on the path");
        let s = XSentSubject::new(subject).expect("subject on the path");
        let raw = &subject[1..subject.len() - 1];
        let expect = extend_two_hit(ScoringMatrix::Blosum62, block, raw, sl, sr, qr, d, 3);
        // SAFETY: d is in 1..=127 (checked by `query`).
        let got = unsafe { two_hit(&q, s, sl, sr, qr, d) };
        assert_eq!(got, expect, "two_hit sl={sl} sr={sr} qr={qr} d={d}");
    }

    /// Random calls on one pair; `block` and `subject` as in [`check_call`].
    fn check_pair(block: &[u8], subject: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
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
            check_call(block, subject, sl, sr, qr, d);
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
            let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
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
}
