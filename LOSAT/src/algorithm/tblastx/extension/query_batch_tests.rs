//! Tests of the walks past a context's ends (query_batch.rs): the extensions of `two_hit.rs` and
//! `one_hit.rs` on a context of a batch against NCBI's functions on the batch's buffer.

use super::super::one_hit::{extend_hit_one_hit, extend_hit_one_hit_alone, OneHitAlone};
use super::super::two_hit::{extend_hit_two_hit_reference, two_hit_walks};
use super::*;
use crate::algorithm::tblastx::scoring::Blosum62Score;
use crate::config::ScoringMatrix;
use crate::stats::KarlinParams;
use crate::utils::matrix::BLAST_SCORE_MIN;

/// Scores large enough that a few letters make up for a NULLB: `MATCH` for two equal letters,
/// `MISMATCH` otherwise, and BLAST_SCORE_MIN against the NULLB (code 0), as NCBI's matrix has it.
#[derive(Clone, Copy)]
struct Big;

const MATCH: i32 = 3000;
const MISMATCH: i32 = -1000;

impl TbxScore for Big {
    const BLOSUM62: bool = false;

    fn score(self, a: u8, b: u8) -> i32 {
        if a == NULLB || b == NULLB {
            BLAST_SCORE_MIN
        } else if a == b {
            MATCH
        } else {
            MISMATCH
        }
    }

    fn matrix(self) -> ScoringMatrix {
        ScoringMatrix::Blosum62
    }
}

/// Contexts with the letters `letters`, with `gaps[k]` NULLBs between context `k` and `k + 1`
/// (one in every batch LOSAT and NCBI make; more for the general layout of `walk_left` and
/// `walk_right`).
fn contexts_of(letters: &[Vec<u8>], gaps: &[usize]) -> Vec<QueryContext> {
    let mut frame_base = 0usize;
    letters
        .iter()
        .enumerate()
        .map(|(k, residues)| {
            let mut aa_seq = vec![NULLB];
            aa_seq.extend_from_slice(residues);
            aa_seq.push(NULLB);
            let context = QueryContext {
                q_idx: 0,
                frame: 1,
                aa_seq,
                aa_seq_nomask: None,
                aa_len: residues.len(),
                orig_len: 3 * residues.len() + 2,
                frame_base: frame_base as i32,
                is_valid: true,
                karlin_params: KarlinParams::default(),
            };
            frame_base += residues.len() + gaps.get(k).copied().unwrap_or(1);
            context
        })
        .collect()
}

/// NCBI's `query->sequence` of the batch (NULLB between the contexts) and `query->length`.
fn batch_buffer(contexts: &[QueryContext]) -> (Vec<u8>, i64) {
    let last = contexts.last().unwrap();
    let length = last.frame_base as usize + last.aa_len;
    let mut buffer = vec![NULLB; length + 1];
    for c in contexts {
        let at = c.frame_base as usize;
        buffer[at..at + c.aa_len].copy_from_slice(&c.aa_seq[1..1 + c.aa_len]);
    }
    (buffer, length as i64)
}

/// Whether batch offset `p` lies outside the letters of context `ctx`.
fn outside(contexts: &[QueryContext], ctx: usize, p: i64) -> bool {
    let c = &contexts[ctx];
    p < i64::from(c.frame_base) || p >= i64::from(c.frame_base) + c.aa_len as i64
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:886-920, as written (on the
// batch's buffer), with the best's position kept to tell a best outside the context.
/// `s_BlastAaExtendLeft`: (maxscore, length, whether a best lies outside context `ctx`).
#[allow(clippy::too_many_arguments)]
fn ncbi_left<S: TbxScore>(
    sc: S,
    q: &[u8],
    s: &[u8],
    s_off: i64,
    q_off: i64,
    dropoff: i64,
    maxscore: i64,
    contexts: &[QueryContext],
    ctx: usize,
) -> (i64, i64, bool) {
    let n = s_off.min(q_off);
    let mut best_i = n + 1;
    let (mut score, mut maxscore) = (maxscore, maxscore);
    let mut crossed = false;
    let mut i = n;
    while i >= 0 {
        let (qp, sp) = (q_off - n + i, s_off - n + i);
        score += i64::from(sc.score(q[qp as usize], s[sp as usize]));
        if score > maxscore {
            maxscore = score;
            best_i = i;
            crossed |= outside(contexts, ctx, qp);
        }
        if (maxscore - score) >= dropoff {
            break;
        }
        i -= 1;
    }
    (maxscore, n - best_i + 1, crossed)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:831-866, as written.
/// `s_BlastAaExtendRight`: (maxscore, length, s_last_off, whether a best lies outside context
/// `ctx`).
#[allow(clippy::too_many_arguments)]
fn ncbi_right<S: TbxScore>(
    sc: S,
    q: &[u8],
    q_length: i64,
    s: &[u8],
    s_off: i64,
    q_off: i64,
    dropoff: i64,
    maxscore: i64,
    contexts: &[QueryContext],
    ctx: usize,
) -> (i64, i64, i64, bool) {
    let n = (s.len() as i64 - s_off).min(q_length - q_off);
    let mut best_i = -1;
    let (mut score, mut maxscore) = (maxscore, maxscore);
    let mut crossed = false;
    let mut i = 0;
    while i < n {
        score += i64::from(sc.score(q[(q_off + i) as usize], s[(s_off + i) as usize]));
        if score > maxscore {
            maxscore = score;
            best_i = i;
            crossed |= outside(contexts, ctx, q_off + i);
        }
        if score <= 0 || (maxscore - score) >= dropoff {
            break;
        }
        i += 1;
    }
    (maxscore, best_i + 1, s_off + i, crossed)
}

/// NCBI's two-hit result in the context's offsets: (q_start, q_end, s_start, s_end, score,
/// right_extend, s_last_off), or `None` where its HSP has a best outside the context.
type Model = Option<(usize, usize, usize, usize, i32, bool, usize)>;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158, as written.
/// `s_BlastAaExtendTwoHit` on the batch from the word at context `ctx`'s offset `q_right`.
#[allow(clippy::too_many_arguments)]
fn ncbi_two_hit<S: TbxScore>(
    sc: S,
    contexts: &[QueryContext],
    ctx: usize,
    s: &[u8],
    s_left: usize,
    s_right: usize,
    q_right: usize,
    dropoff: i32,
) -> Model {
    let (q, q_length) = batch_buffer(contexts);
    let fb = i64::from(contexts[ctx].frame_base);
    let (mut q_right_off, mut s_right_off) = (fb + q_right as i64, s_right as i64);
    let (mut score, mut left_score, mut right_d) = (0i64, 0i64, 0i64);
    for i in 0..sc.word() as i64 {
        score += i64::from(sc.score(q[(q_right_off + i) as usize], s[(s_right_off + i) as usize]));
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;
    let mut right_d = 0;
    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    let dropoff = i64::from(dropoff);
    let (left_score, left_d, crossed_left) = ncbi_left(
        sc,
        &q,
        s,
        s_right_off - 1,
        q_right_off - 1,
        dropoff,
        0,
        contexts,
        ctx,
    );
    let mut right_score = 0;
    let mut crossed_right = false;
    if left_d >= s_right_off - s_left as i64 {
        right_extend = true;
        let (score, length, last, crossed) = ncbi_right(
            sc,
            &q,
            q_length,
            s,
            s_right_off,
            q_right_off,
            dropoff,
            left_score,
            contexts,
            ctx,
        );
        right_score = score;
        right_d = length;
        s_last_off = last;
        crossed_right = crossed;
    }
    if crossed_left || crossed_right {
        return None;
    }
    let q_start = q_right_off - left_d - fb;
    let s_start = s_right_off - left_d;
    let len = left_d + right_d;
    Some((
        q_start as usize,
        (q_start + len) as usize,
        s_start as usize,
        (s_start + len) as usize,
        left_score.max(right_score) as i32,
        right_extend,
        s_last_off as usize,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1020-1085, as written.
/// `s_BlastAaExtendOneHit` on the batch from the word at context `ctx`'s offset `q_word`:
/// (q_start in the context, s_start, length, score, s_last_off), or `None` where its HSP has a
/// best outside the context.
fn ncbi_one_hit<S: TbxScore>(
    sc: S,
    contexts: &[QueryContext],
    ctx: usize,
    s: &[u8],
    q_word: usize,
    s_off: usize,
    dropoff: i32,
) -> Option<(i64, i64, i64, i64, i64)> {
    let (q, q_length) = batch_buffer(contexts);
    let fb = i64::from(contexts[ctx].frame_base);
    let (q_off, s_off, word) = (fb + q_word as i64, s_off as i64, sc.word() as i64);
    let (mut score, mut sum) = (0i64, 0i64);
    let (mut q_left_off, mut q_right_off, mut q_best_left_off) = (q_off, q_off + word, q_off);
    for i in 0..word {
        sum += i64::from(sc.score(q[(q_off + i) as usize], s[(s_off + i) as usize]));
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
    let q_left_off = q_best_left_off;
    let s_left_off = q_left_off + (s_off - q_off);
    let s_right_off = q_right_off + (s_off - q_off);
    let dropoff = i64::from(dropoff);
    let (left_score, left_disp, crossed_left) = ncbi_left(
        sc,
        &q,
        s,
        s_left_off - 1,
        q_left_off - 1,
        dropoff,
        score,
        contexts,
        ctx,
    );
    let (total_score, right_disp, s_last_off, crossed_right) = ncbi_right(
        sc,
        &q,
        q_length,
        s,
        s_right_off + 1,
        q_right_off + 1,
        dropoff,
        left_score,
        contexts,
        ctx,
    );
    if crossed_left || crossed_right {
        return None;
    }
    Some((
        q_left_off - left_disp - fb,
        s_left_off - left_disp,
        left_disp + right_disp + init_hit_width,
        total_score,
        s_last_off,
    ))
}

/// A small xorshift generator for the random batches.
struct Rng(u64);

impl Rng {
    fn below(&mut self, n: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % n
    }
}

/// A random batch of 1 to 5 contexts over `alphabet` letters (from 1), with `gap_max` NULLBs at
/// most between two of them, a subject that copies the batch's buffer from a random offset (so
/// that diagonals run on across the NULLBs), with random letters changed, and that offset.
fn random_case(rng: &mut Rng, alphabet: u64, gap_max: u64) -> (Vec<QueryContext>, Vec<u8>, usize) {
    let count = 1 + rng.below(5) as usize;
    let letters: Vec<Vec<u8>> = (0..count)
        .map(|_| {
            let len = 3 + rng.below(40) as usize;
            (0..len).map(|_| 1 + rng.below(alphabet) as u8).collect()
        })
        .collect();
    let gaps: Vec<usize> = (0..count)
        .map(|_| 1 + rng.below(gap_max) as usize)
        .collect();
    let contexts = contexts_of(&letters, &gaps);
    let (buffer, length) = batch_buffer(&contexts);
    let from = rng.below(length as u64) as usize;
    let mut subject: Vec<u8> = buffer[from..length as usize]
        .iter()
        .map(|&c| {
            if c == NULLB {
                1 + rng.below(alphabet) as u8
            } else {
                c
            }
        })
        .collect();
    let extra = rng.below(30) as usize;
    subject.extend((0..extra).map(|_| 1 + rng.below(alphabet) as u8));
    for _ in 0..rng.below(1 + subject.len() as u64 / 4) {
        let at = rng.below(subject.len() as u64) as usize;
        subject[at] = 1 + rng.below(alphabet) as u8;
    }
    if subject.len() < 3 {
        subject.extend([1, 2, 3]);
    }
    (contexts, subject, from)
}

/// A subject offset for a word at batch offset `q`: half of the time on the copy's diagonal.
fn subject_offset(rng: &mut Rng, subject: &[u8], q: usize, from: usize) -> usize {
    match q.checked_sub(from) {
        Some(s) if s + 3 <= subject.len() && rng.below(2) == 0 => s,
        _ => rng.below((subject.len() - 3) as u64 + 1) as usize,
    }
}

const X_DROPS: [i32; 9] = [1, 16, 5000, 32767, 32768, 32769, 40000, 100_000, 50_000_000];

// No NCBI counterpart: the two-hit extension of a context of a batch equals NCBI's on the batch's
// buffer where NCBI's HSP stays in the context (the left walk's continuation changes nothing
// there; the right walk's moves `s_last_off`), and the search is refused exactly where NCBI's best
// lies in another context; for random batches with up to three NULLBs between contexts (the
// general layout; LOSAT's and NCBI's batches have one), X-drops below, at and above 32,768, and
// BLOSUM62's scores beside the large ones of `Big`.
#[test]
fn the_two_hit_walks_go_on_into_the_batch_as_ncbis_do() {
    let mut rng = Rng(0x2c2c_9e37_79b9_7f4a);
    let (mut moved, mut refused, mut equal, mut again) = (0usize, 0usize, 0usize, 0usize);
    for round in 0..30_000 {
        let blosum = round % 3 == 0;
        let (contexts, subject, from) = random_case(&mut rng, if blosum { 20 } else { 3 }, 3);
        let ctx = rng.below(contexts.len() as u64) as usize;
        let len = contexts[ctx].aa_len;
        let q_right = rng.below((len - 3) as u64 + 1) as usize;
        let fb = contexts[ctx].frame_base as usize;
        let s_right = subject_offset(&mut rng, &subject, fb + q_right, from);
        let s_left = s_right.saturating_sub(rng.below(20) as usize);
        let x_drop = X_DROPS[rng.below(X_DROPS.len() as u64) as usize];
        let failure = XReadFailure::new();
        let batch = QueryBatch::new(&contexts, &failure);
        let letters = &contexts[ctx].aa_seq[1..1 + len];
        macro_rules! run {
            ($sc:expr) => {{
                let got = extend_hit_two_hit_reference(
                    $sc,
                    letters,
                    batch.at(ctx, 7),
                    &subject,
                    s_left,
                    s_right,
                    q_right,
                    x_drop,
                    false,
                );
                let alone = extend_hit_two_hit_reference(
                    $sc,
                    letters,
                    InBatch::ALONE,
                    &subject,
                    s_left,
                    s_right,
                    q_right,
                    x_drop,
                    false,
                );
                let model = ncbi_two_hit(
                    $sc, &contexts, ctx, &subject, s_left, s_right, q_right, x_drop,
                );
                // The search's first pass, on the context alone.
                let first = two_hit_walks(
                    $sc,
                    letters,
                    InBatch::ALONE,
                    &subject,
                    s_left,
                    s_right,
                    q_right,
                    x_drop,
                    false,
                    true,
                );
                (got, alone, model, first)
            }};
        }
        let (got, alone, model, first) = if blosum {
            run!(Blosum62Score)
        } else {
            run!(Big)
        };
        // A refusal fails the search of the subject the extension reads (7).
        if failure.failed() {
            assert!(failure.skips(8) && !failure.skips(7), "round {round}");
        }
        let was_refused = failure.take().is_some();
        // The first pass gives the batch's result, or asks for the batch only with an X-drop
        // above 32,768.
        match first {
            (r, false) => {
                assert_eq!(r, got, "round {round}: first pass");
                assert!(!was_refused, "round {round}: first pass");
            }
            (_, true) => {
                assert!(x_drop > -BLAST_SCORE_MIN, "round {round}: again");
                again += 1;
            }
        }
        match model {
            Some(expect) => {
                assert!(!was_refused, "round {round}: refused, NCBI {expect:?}");
                assert_eq!(got, expect, "round {round}");
                if got == alone {
                    equal += 1;
                } else {
                    // Only the stop of the right walk moves.
                    assert_eq!(
                        (got.0, got.1, got.2, got.3, got.4, got.5),
                        (alone.0, alone.1, alone.2, alone.3, alone.4, alone.5),
                        "round {round}"
                    );
                    assert!(got.6 > alone.6, "round {round}");
                    moved += 1;
                }
            }
            None => {
                assert!(
                    was_refused,
                    "round {round}: NCBI's best lies in another context"
                );
                refused += 1;
            }
        }
        // With 32,768 or less the NULLB breaks every walk (`maxscore - score >= 32768` there).
        if x_drop <= -BLAST_SCORE_MIN {
            assert_eq!(got, alone, "round {round}: x_drop {x_drop}");
            assert!(!was_refused, "round {round}");
        }
    }
    assert!(
        moved > 100 && refused > 100 && equal > 1000 && again > 100,
        "{moved} {refused} {equal} {again}"
    );
}

// No NCBI counterpart: the one-hit extension (`-window_size 0`) of a context of a batch against
// NCBI's `s_BlastAaExtendOneHit` on the batch's buffer, as above; with BLOSUM62 and every byte at
// a low threshold also the words whose running sum never rises above 0 at a context's end, whose
// right walk NCBI starts on the next context's first letter (past the NULLB) at any X-drop.
#[test]
fn the_one_hit_walks_go_on_into_the_batch_as_ncbis_do() {
    let mut rng = Rng(0x0e17_1c0d_e5a5_5a5a);
    let (mut moved, mut refused, mut past_end, mut again) = (0usize, 0usize, 0usize, 0usize);
    for round in 0..30_000 {
        let blosum = round % 2 == 0;
        let (contexts, subject, from) = random_case(&mut rng, if blosum { 25 } else { 3 }, 3);
        let ctx = rng.below(contexts.len() as u64) as usize;
        let len = contexts[ctx].aa_len;
        // Words at a context's end half of the time.
        let q_word = if rng.below(2) == 0 {
            len - 3
        } else {
            rng.below((len - 3) as u64 + 1) as usize
        };
        let fb = contexts[ctx].frame_base as usize;
        let s_off = subject_offset(&mut rng, &subject, fb + q_word, from);
        let x_drop = X_DROPS[rng.below(X_DROPS.len() as u64) as usize];
        let failure = XReadFailure::new();
        let batch = QueryBatch::new(&contexts, &failure);
        let letters = &contexts[ctx].aa_seq[1..1 + len];
        macro_rules! run {
            ($sc:expr) => {{
                let got = extend_hit_one_hit(
                    $sc,
                    letters,
                    batch.at(ctx, 7),
                    &subject,
                    q_word,
                    s_off,
                    x_drop,
                )
                .expect("the word fits");
                let alone = extend_hit_one_hit(
                    $sc,
                    letters,
                    InBatch::ALONE,
                    &subject,
                    q_word,
                    s_off,
                    x_drop,
                )
                .expect("the word fits");
                let model = ncbi_one_hit($sc, &contexts, ctx, &subject, q_word, s_off, x_drop);
                // The search's first pass, on the context alone.
                let first = extend_hit_one_hit_alone($sc, letters, &subject, q_word, s_off, x_drop)
                    .expect("the word fits");
                (got, alone, model, first)
            }};
        }
        let (got, alone, model, first) = if blosum {
            run!(Blosum62Score)
        } else {
            run!(Big)
        };
        // A refusal fails the search of the subject the extension reads (7).
        if failure.failed() {
            assert!(failure.skips(8) && !failure.skips(7), "round {round}");
        }
        let was_refused = failure.take().is_some();
        // The first pass gives the batch's result, or asks for the batch only where a walk ends
        // at the context's end with an X-drop above 32,768 or starts past it (a word at the
        // context's end that never scores above 0).
        match first {
            OneHitAlone::Hit(hit) => {
                assert_eq!(hit, got, "round {round}: first pass");
                assert!(!was_refused, "round {round}: first pass");
            }
            OneHitAlone::Again(_) => {
                assert!(
                    x_drop > -BLAST_SCORE_MIN || q_word + 3 == len,
                    "round {round}: again"
                );
                again += 1;
            }
        }
        match model {
            Some((q_start, s_start, length, score, s_last_off)) => {
                assert!(!was_refused, "round {round}: refused");
                assert_eq!(
                    (
                        i64::from(got.q_start),
                        i64::from(got.s_start),
                        i64::from(got.length),
                        i64::from(got.score),
                        i64::from(got.s_last_off)
                    ),
                    (q_start, s_start, length, score, s_last_off),
                    "round {round}"
                );
                if got != alone {
                    assert!(got.s_last_off > alone.s_last_off, "round {round}");
                    moved += 1;
                    if x_drop <= -BLAST_SCORE_MIN {
                        past_end += 1;
                    }
                }
            }
            None => {
                assert!(
                    was_refused,
                    "round {round}: NCBI's best lies in another context"
                );
                refused += 1;
            }
        }
    }
    assert!(
        moved > 100 && refused > 100 && past_end > 0 && again > 100,
        "{moved} {refused} {past_end} {again}"
    );
}

// No NCBI counterpart: `NULLB_SCORE` is the NULLB's score against every byte in BLOSUM62's
// instance and in every matrix's (`nullb_goes_on_left` and `nullb_goes_on_right` take the NULLB's
// step from it).
#[test]
fn the_nullb_scores_blast_score_min_against_every_byte() {
    for b in 0..=u8::MAX {
        assert_eq!(Blosum62Score.score(NULLB, b), NULLB_SCORE, "{b}");
    }
    for matrix in ScoringMatrix::ALL {
        let scores = crate::algorithm::tblastx::scoring::MatrixScore::new(matrix);
        for b in 0..=u8::MAX {
            assert_eq!(scores.score(NULLB, b), NULLB_SCORE, "{matrix:?} {b}");
            assert_eq!(matrix.score(NULLB, b), NULLB_SCORE, "{matrix:?} {b}");
        }
    }
}

// No NCBI counterpart: the NULLB's step of NCBI's walks against the tests of `left_at_bound` and
// `right_at_bound` on a context alone, at every distance below the best that a walk without a
// break can be at, around the X-drops where it changes.
#[test]
fn the_nullb_step_breaks_the_walks_as_ncbis() {
    for x_drop in [1, 2, 100, 32767, 32768, 32769, 32770, 40000, i32::MAX] {
        for max_score in [1, 5, 32768, 32769, 40000, 1 << 30] {
            for below in [0, 1, 2, x_drop / 2, x_drop - 2, x_drop - 1] {
                if below < 0 || below >= x_drop {
                    continue;
                }
                let score = max_score - below;
                let next = i64::from(score) + i64::from(BLAST_SCORE_MIN);
                let left = i64::from(max_score) - next >= i64::from(x_drop);
                let right = next <= 0 || left;
                assert_eq!(nullb_goes_on_left(score, max_score, x_drop), !left);
                assert_eq!(nullb_goes_on_right(score, max_score, x_drop), !right);
            }
        }
    }
}

// No NCBI counterpart: the cases of the brief, one by one (Big's scores: a NULLB costs 32,768,
// eleven equal letters make it up).
#[test]
fn the_walks_past_a_context_end_case_by_case() {
    let failure = XReadFailure::new();
    // Context 0 holds 12 letters 1, context 1 starts with 2 2 2 2; the subject's diagonal through
    // context 1's start meets five 1s before the NULLB's place (3) and mismatches before them.
    let contexts = contexts_of(&[vec![1; 12], vec![2, 2, 2, 2, 3, 3, 3]], &[1]);
    let subject: Vec<u8> = [vec![4; 7], vec![1; 5], vec![3], vec![2, 2, 2, 2]].concat();
    let batch = QueryBatch::new(&contexts, &failure);
    let letters1 = &contexts[1].aa_seq[1..8];
    // From the word 2 2 2 at context 1's start the left walk (9,000) passes the NULLB with an
    // X-drop of 40,000 and of 1,000,000 and runs to the subject's start without a new best: the
    // result of the context alone.
    for x_drop in [40_000, 1_000_000] {
        let alone = extend_hit_two_hit_reference(
            Big,
            letters1,
            InBatch::ALONE,
            &subject,
            13,
            13,
            0,
            x_drop,
            false,
        );
        let got = extend_hit_two_hit_reference(
            Big,
            letters1,
            batch.at(1, 7),
            &subject,
            13,
            13,
            0,
            x_drop,
            false,
        );
        assert_eq!(got, alone, "{x_drop}");
        assert!(failure.take().is_none());
    }
    // With twelve 1s before the NULLB's place the walk gains past it at the eleventh (9,000 -
    // 32,768 + 33,000 > 9,000): refused.
    let subject: Vec<u8> = [vec![1; 12], vec![3], vec![2, 2, 2, 2]].concat();
    let _ = extend_hit_two_hit_reference(
        Big,
        letters1,
        batch.at(1, 7),
        &subject,
        13,
        13,
        0,
        40_000,
        false,
    );
    let error = failure.take().expect("refused");
    assert!(error
        .to_string()
        .ends_with("is not supported by LOSAT's TBLASTX"));
    // From context 0's word at 0: the right walk reaches the NULLB with 36,000, passes it with an
    // X-drop above 32,768 (score 3,232 > 0) and runs over 2 2 2 2 (no new best) to the subject's
    // end: `s_last_off` 17 in place of 12.
    let letters0 = &contexts[0].aa_seq[1..13];
    let alone = extend_hit_two_hit_reference(
        Big,
        letters0,
        InBatch::ALONE,
        &subject,
        0,
        0,
        0,
        40_000,
        false,
    );
    let got = extend_hit_two_hit_reference(
        Big,
        letters0,
        batch.at(0, 7),
        &subject,
        0,
        0,
        0,
        40_000,
        false,
    );
    assert_eq!(alone.6, 12);
    assert_eq!(
        got,
        (alone.0, alone.1, alone.2, alone.3, alone.4, alone.5, 17)
    );
    assert!(failure.take().is_none());
    // An X-drop of 32,768 stops it at the NULLB, as alone.
    let got = extend_hit_two_hit_reference(
        Big,
        letters0,
        batch.at(0, 7),
        &subject,
        0,
        0,
        0,
        32_768,
        false,
    );
    let alone = extend_hit_two_hit_reference(
        Big,
        letters0,
        InBatch::ALONE,
        &subject,
        0,
        0,
        0,
        32_768,
        false,
    );
    assert_eq!(got, alone);
    // Twelve letters 2 in context 1 make the right walk score higher past the NULLB: refused.
    let contexts = contexts_of(&[vec![1; 12], vec![2; 12]], &[1]);
    let subject: Vec<u8> = [vec![1; 12], vec![3], vec![2; 12]].concat();
    let batch = QueryBatch::new(&contexts, &failure);
    let letters0 = &contexts[0].aa_seq[1..13];
    let _ = extend_hit_two_hit_reference(
        Big,
        letters0,
        batch.at(0, 7),
        &subject,
        0,
        0,
        0,
        40_000,
        false,
    );
    assert!(failure.take().is_some());
    // Three NULLBs (a gap the batches never have): the walk pays three times, and 12 letters of
    // context 0 do not make up for them from context 2's start (36,000 < 3 x 32,768).
    let contexts = contexts_of(&[vec![1; 12], vec![2, 2, 2, 2]], &[3]);
    let subject: Vec<u8> = [vec![1; 12], vec![3, 3, 3], vec![2, 2, 2, 2]].concat();
    let batch = QueryBatch::new(&contexts, &failure);
    let letters1 = &contexts[1].aa_seq[1..5];
    let _ = extend_hit_two_hit_reference(
        Big,
        letters1,
        batch.at(1, 7),
        &subject,
        15,
        15,
        0,
        1_000_000_000,
        false,
    );
    assert!(failure.take().is_none());
    // ... while 40 letters do.
    let contexts = contexts_of(&[vec![1; 40], vec![2, 2, 2, 2]], &[3]);
    let subject: Vec<u8> = [vec![1; 40], vec![3, 3, 3], vec![2, 2, 2, 2]].concat();
    let batch = QueryBatch::new(&contexts, &failure);
    let letters1 = &contexts[1].aa_seq[1..5];
    let _ = extend_hit_two_hit_reference(
        Big,
        letters1,
        batch.at(1, 7),
        &subject,
        43,
        43,
        0,
        1_000_000_000,
        false,
    );
    assert!(failure.take().is_some());
    // The batch's first context's left end and its last context's right end are the batch's ends.
    let contexts = contexts_of(&[vec![2; 12]], &[1]);
    let subject: Vec<u8> = [vec![2; 30], vec![2; 12], vec![2; 30]].concat();
    let batch = QueryBatch::new(&contexts, &failure);
    let letters0 = &contexts[0].aa_seq[1..13];
    let got = extend_hit_two_hit_reference(
        Big,
        letters0,
        batch.at(0, 7),
        &subject,
        30,
        30,
        0,
        1_000_000_000,
        false,
    );
    let alone = extend_hit_two_hit_reference(
        Big,
        letters0,
        InBatch::ALONE,
        &subject,
        30,
        30,
        0,
        1_000_000_000,
        false,
    );
    assert_eq!(got, alone);
    assert!(failure.take().is_none());
}
