//! The query batch around the context whose letters an ungapped extension walks.
//!
//! NCBI's protein walks run on the query batch, not on a context:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:841 (`s_BlastAaExtendRight`)
//! and :897 (`s_BlastAaExtendLeft`)
//! `query` is the batch's sequence block: `query->sequence` holds the letters of every context of
//! the batch, `q_off` is an offset into it and `query->length` its length (two less than `buflen`,
//! api/blast_setup_cxx.cpp:654, with `QueryInfo_GetSeqBufLen` of blast_query_info.c:246-250: the
//! end of the last context's letters). The contexts lie one after the other with one NULLB between
//! two of them (api/blast_setup_cxx.cpp:65-88: `shift` is `prev_len` plus 1 when `prev_len` is not
//! 0, else 0, so a context without letters takes no room; LOSAT makes no such context,
//! `translation::generate_frames`, and `QueryContext::frame_base` is NCBI's `query_offset`,
//! lookup/backbone.rs). The NULLB's row scores `BLAST_SCORE_MIN` (-32768) against every letter
//! (blast_stat.c:1559-1576, kept for `-`).
//!
//! So a walk that reaches the end of its context without a break reads the NULLB next to it with
//! the same tests (left: `(maxscore - score) >= dropoff`; right: `score <= 0 ||
//! (maxscore - score) >= dropoff`), and while those do not break, the next context's letters (and
//! beyond it the next NULLB and context), until NCBI's `n` runs out: at the subject's end or start
//! (a subject frame, or a chunk of one, is the whole `subject`), at the batch's first letter on the
//! left and at the end of its last context's letters on the right.
//!
//! LOSAT's walks (`two_hit.rs`) read the context's letters (`&aa_seq[1..len - 1]`); a walk that
//! reads the last one on its side without a break, where the subject has letters left, goes on here
//! (`InBatch::left_at_bound`, `InBatch::right_at_bound`). With a raw X-drop of 32,768 or less the
//! NULLB always breaks the walk (`maxscore - score >= 32768`), on the left with the same best and
//! length, on the right also with the same stop: `*s_last_off = s_off + i` with `i` the NULLB's
//! step, the count of the walk that ended at the context's end.
//!
//! The searches so walk a context alone first (`InBatch::ALONE`, `Pass::First`): a walk that ends at
//! the context's end tests, from its own sum, best and X-drop, whether the NULLB's step would break
//! it (`nullb_goes_on_left`, `nullb_goes_on_right`), and only where it would not does the search
//! take the batch; the walks' loops keep no other values. The one-hit extension goes on from the
//! walk's end (`Pass::Resume` with `WalkEnd`; the scan calls its second pass itself,
//! `word_finder::one_hit_in_batch`, so that it records the first pass's result from its registers:
//! with one result for both passes the first pass's went through memory at every extension, X2c,
//! `-window_size 0`, +1.4 % instructions); the two-hit extension walks both again on the batch
//! (`extend_hit_two_hit_reference`).
//!
//! A walk past the NULLB that scores no higher than its best leaves the extension's result as it
//! was, except that the right walk's stop (`s_last_off`, and through it the diagonal's `last_hit`
//! for the later hits of that diagonal, aa_ungapped.c:593-597) moves past the NULLB: that is NCBI's
//! output, which LOSAT reproduces (PD-LOSAT-NCBI-DEFECTS, rule 2). A walk that scores higher there
//! moves the HSP's start into the previous context or its end into the next: NCBI keeps the HSP in
//! the seed's context with offsets outside it (a start before the context fails the assertion of
//! `s_AdjustInitialHSPOffsets`, which is compiled out of release builds, blast_gapalign.c:2390-2391;
//! an end past it is not tested), and its later steps write query coordinates past the query or
//! drop the alignment from every report. LOSAT refuses such a search (PD-LOSAT-NCBI-DEFECTS, rule
//! 3, the Owner's decision of 2026-11-16: `refusal`), after the reports of the query batches before.
//!
//! The one-hit extension starts its right walk one letter past a word whose running score never
//! rises above 0 (aa_ungapped.c:1030-1031, 1073-1076), past the context's end for such a word at
//! the end of its context; the walk then starts on the next context's first letter
//! (`InBatch::right_at_bound` with `q_at` one past the NULLB). Such a word scores below every
//! `-threshold` and seeds only as an exact match (blast_aalookup.c:498-504), of a subject word of
//! the same letters, and only `X` scores below 0 against itself; the seed stage's subject frames
//! are translated from the subject's 2-bit bases (blast_engine.c:1403-1405, `subjects.rs`) and
//! hold no `X`, so no such seed is found.

use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::scoring::TbxScore;
use crate::algorithm::tblastx::subjects::XReadFailure;
use crate::utils::matrix::BLAST_SCORE_MIN;
use std::sync::atomic::{AtomicBool, Ordering};

/// The byte between two contexts of NCBI's query batch (`NULLB`, ncbi_std.h), which is also each
/// context buffer's sentinel in LOSAT.
const NULLB: u8 = 0;

/// The NULLB's score against every byte in every matrix of the search: `BLAST_SCORE_MIN`, the row
/// and the column of `-` (blast_stat.c:1559-1576; `utils::matrix::protein_score`, and BLOSUM62's
/// direct table; test `the_nullb_scores_blast_score_min_against_every_byte`).
const NULLB_SCORE: i32 = BLAST_SCORE_MIN;

/// The query batch of a search, as its two-hit and one-hit extensions read it: the batch's
/// contexts, and where the search reports an extension that LOSAT refuses. The frame loop makes it
/// once per batch (`WordFinderCx::query_batch`), the unit path once per hit context (`XHitCtx`,
/// `noting`); an extension reads it through `QueryBatch::at`, with its subject. `Copy`: it holds references only (a subject's
/// `WordFinderCx` copies the batch's).
#[derive(Clone, Copy)]
pub(crate) struct QueryBatch<'a> {
    contexts: &'a [QueryContext],
    failure: &'a XReadFailure,
    /// The unit path's note of a right walk that went on past a NULLB (`QueryBatch::noting`;
    /// x_pair_par.rs, `x_units_in_order`), or None.
    passed: Option<&'a AtomicBool>,
}

impl<'a> QueryBatch<'a> {
    /// The batch's `contexts` (in the order of their `frame_base`), whose searches refuse into
    /// `failure`.
    pub(crate) fn new(contexts: &'a [QueryContext], failure: &'a XReadFailure) -> Self {
        Self::noting(contexts, failure, None)
    }

    /// As `new`, and where `passed` is given, every right walk that goes on past a NULLB sets it.
    pub(crate) fn noting(
        contexts: &'a [QueryContext],
        failure: &'a XReadFailure,
        passed: Option<&'a AtomicBool>,
    ) -> Self {
        Self {
            contexts,
            failure,
            passed,
        }
    }

    /// Context `ctx` of the batch, as an extension in it against subject `s_idx` reads the batch.
    #[inline(always)]
    pub(crate) fn at(&'a self, ctx: usize, s_idx: usize) -> InBatch<'a> {
        debug_assert!(ctx <= u32::MAX as usize && s_idx <= u32::MAX as usize);
        InBatch {
            batch: Some(self),
            // A batch's contexts and a search's subjects are counted in `u32`, as `InitHSP`
            // counts the subject.
            at: (ctx as u32 as u64) | ((s_idx as u32 as u64) << 32),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:897-916
    // The steps of NCBI's left walk after the context's first letter: `s_next` subject letters are
    // left (the next one is `s_seq[s_next - 1]`), and the batch's letters before the context are its
    // `frame_base`. The sums are `i64` (NCBI's `Int4` sum of many NULLBs would overflow).
    /// The left walk of context `ctx` past its start; refuses the search where it scores higher.
    #[cold]
    #[inline(never)]
    fn walk_left<S: TbxScore>(
        &self,
        ctx: usize,
        s_idx: usize,
        sc: S,
        s_seq: &[u8],
        s_next: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) {
        let contexts = self.contexts;
        let start = contexts[ctx].frame_base as usize;
        let n = s_next.min(start);
        let (mut score, max_score, x_drop) =
            (i64::from(score), i64::from(max_score), i64::from(x_drop));
        // Context by context down from `ctx`: the NULLBs between the walk's position and the
        // context's letters, then its letters from the last one (each run in a loop of its own).
        let (mut step, mut k) = (0usize, ctx);
        while step < n && k > 0 {
            k -= 1;
            let context = &contexts[k];
            let len = context.aa_len;
            let nullbs = (start - step)
                .saturating_sub(context.frame_base as usize + len)
                .min(n - step);
            for _ in 0..nullbs {
                score += i64::from(sc.score(NULLB, s_seq[s_next - 1 - step]));
                step += 1;
                if score > max_score {
                    self.refuse(s_idx);
                    return;
                }
                if max_score - score >= x_drop {
                    return;
                }
            }
            let take = len.min(n - step);
            let letters = &context.aa_seq[1 + len - take..1 + len];
            let subject = &s_seq[s_next - step - take..s_next - step];
            for (&q, &s) in letters.iter().rev().zip(subject.iter().rev()) {
                score += i64::from(sc.score(q, s));
                if score > max_score {
                    self.refuse(s_idx);
                    return;
                }
                if max_score - score >= x_drop {
                    return;
                }
            }
            step += take;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:841-866
    // The steps of NCBI's right walk after the context's letters: the walk has taken `j` steps, the
    // next subject letter is `s_seq[s_at]` and the next query byte is at offset `q_at` of the
    // context (its letters' count, the NULLB after them; or one more, where the one-hit extension
    // starts its walk past its context's end, `extend_hit_one_hit`), and the batch ends at its last
    // context's letters (`query->length`). Returns NCBI's `i` at the end of the walk.
    /// The right walk of context `ctx` past its end, and its stop; refuses the search where it
    /// scores higher. A walk that goes on past the NULLB after the context's letters sets
    /// `passed` (its stop then lies past the context, where the cell of its diagonal masks the
    /// later hits of the diagonal: x_pair_par.rs, module doc item 3).
    #[cold]
    #[inline(never)]
    fn walk_right<S: TbxScore>(
        &self,
        ctx: usize,
        s_idx: usize,
        sc: S,
        s_seq: &[u8],
        s_at: usize,
        q_at: usize,
        j: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) -> usize {
        let i = self.walk_right_steps(
            ctx, s_idx, sc, s_seq, s_at, q_at, j, score, max_score, x_drop,
        );
        if let (Some(passed), true) = (self.passed, i > j) {
            passed.store(true, Ordering::Relaxed);
        }
        i
    }

    /// The steps of `walk_right`.
    #[inline(always)]
    fn walk_right_steps<S: TbxScore>(
        &self,
        ctx: usize,
        s_idx: usize,
        sc: S,
        s_seq: &[u8],
        s_at: usize,
        q_at: usize,
        j: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) -> usize {
        let contexts = self.contexts;
        let start = contexts[ctx].frame_base as usize + q_at;
        let batch_end = contexts
            .last()
            .map_or(start, |last| last.frame_base as usize + last.aa_len);
        // NCBI's `n` of a walk that starts past the batch's end is negative: no step.
        let n = (s_seq.len() - s_at).min(batch_end.saturating_sub(start));
        let (mut score, max_score, x_drop) =
            (i64::from(score), i64::from(max_score), i64::from(x_drop));
        // Context by context up from `ctx`: the NULLBs before the next context's letters, then its
        // letters (from the walk's position in it; each run in a loop of its own).
        let (mut step, mut k) = (0usize, ctx);
        while step < n && k + 1 < contexts.len() {
            k += 1;
            let context = &contexts[k];
            let first = context.frame_base as usize;
            if start + step >= first + context.aa_len {
                continue;
            }
            let nullbs = first.saturating_sub(start + step).min(n - step);
            for _ in 0..nullbs {
                score += i64::from(sc.score(NULLB, s_seq[s_at + step]));
                if score > max_score {
                    self.refuse(s_idx);
                    return j + step;
                }
                if score <= 0 || max_score - score >= x_drop {
                    return j + step;
                }
                step += 1;
            }
            let off = (start + step).saturating_sub(first);
            let take = (context.aa_len - off).min(n - step);
            let letters = &context.aa_seq[1 + off..1 + off + take];
            let subject = &s_seq[s_at + step..s_at + step + take];
            for (i, (&q, &s)) in letters.iter().zip(subject).enumerate() {
                score += i64::from(sc.score(q, s));
                if score > max_score {
                    self.refuse(s_idx);
                    return j + step + i;
                }
                if score <= 0 || max_score - score >= x_drop {
                    return j + step + i;
                }
            }
            step += take;
        }
        j + n
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2383-2392
    // An HSP that a walk past a NULLB makes span two contexts keeps the seed's context with offsets
    // outside it (a start before it fails the debug-only assertion above; an end past it is not
    // tested), and NCBI's later steps report query coordinates past the query or no alignment
    // (PD-LOSAT-NCBI-DEFECTS, rule 3, the Owner's decision of 2026-11-16).
    /// Refuses the search of subject `s_idx` (and so the query batch's search).
    fn refuse(&self, s_idx: usize) {
        self.failure.fail(s_idx, refusal());
    }
}

/// The error of a search that `QueryBatch` refuses: an ungapped extension that scores higher past
/// the end of its query frame, in the frame next to it in the query batch. It comes after the
/// reports of the query batches before (`ErrorAfterBatchReports`; the Owner's decision of
/// 2026-11-16 on line RG's Q-X2c).
pub(crate) fn refusal() -> anyhow::Error {
    anyhow::Error::new(
        crate::algorithm::tblastx::blast_engine::ErrorAfterBatchReports(anyhow::anyhow!(
            "an ungapped extension that passes the end of a query frame into the frame next to it \
             in the query batch and scores higher there (an -xdrop_ungap of a raw X-drop above \
             32,768, with BLOSUM62 about 15,000 bits or more; NCBI BLAST+ 2.17.0 then keeps an \
             alignment of two frames with query coordinates outside the query's frame and writes \
             them past the query or drops the alignment) is not supported by LOSAT's TBLASTX"
        )),
    )
}

/// Where the letters an ungapped extension walks lie: context `ctx` of a query batch
/// (`QueryBatch::at`), or a context alone (`ALONE`: the walks stop at its ends, as NCBI's do when
/// the NULLB breaks them).
#[derive(Clone, Copy)]
pub(crate) struct InBatch<'a> {
    batch: Option<&'a QueryBatch<'a>>,
    /// The context (the low 32 bits) and the subject, whose search a refusal fails (the high 32
    /// bits): one word, so that the struct is a pair of scalars, which a call passes in two
    /// registers (with three fields it goes through memory at every out-of-line call of the unit
    /// path's extensions).
    at: u64,
}

impl InBatch<'_> {
    /// The context.
    #[inline(always)]
    fn ctx(self) -> usize {
        self.at as u32 as usize
    }

    /// The subject.
    #[inline(always)]
    fn s_idx(self) -> usize {
        (self.at >> 32) as usize
    }
}

impl InBatch<'static> {
    /// A context alone: the walks stop at its ends. The extensions walk so first and on the batch
    /// only where NCBI's walk could go on past an end (`InBatch::left_at_bound`,
    /// `InBatch::right_at_bound`; `extend_hit_one_hit_alone`, `extend_hit_two_hit_reference`);
    /// the tests, and `x_two_hit_cells`, which counts the cells of the same walks.
    pub(crate) const ALONE: Self = Self { batch: None, at: 0 };
}

impl InBatch<'_> {
    /// The left walk read the letters before its start without a break (`q_off` of the context's
    /// letters before it, `s_off` of the subject's: it stopped at the first one of the two), at
    /// `score` below its best `max_score`. On a batch NCBI's walk goes on into it where the
    /// subject has letters more (`QueryBatch::walk_left`; the batch's first context ends here),
    /// and the walk's result stands: false. On a context alone: whether NCBI's walk could go on,
    /// its next byte (the NULLB before the context) not breaking it (`nullb_goes_on_left`): the
    /// search then runs the extension again on the batch; elsewhere the walk's result is NCBI's.
    #[inline(always)]
    pub(super) fn left_at_bound<S: TbxScore>(
        self,
        sc: S,
        s_seq: &[u8],
        q_off: usize,
        s_off: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) -> bool {
        match self.batch {
            Some(batch) => {
                if q_off < s_off && self.ctx() > 0 {
                    batch.walk_left(
                        self.ctx(),
                        self.s_idx(),
                        sc,
                        s_seq,
                        s_off - q_off,
                        score,
                        max_score,
                        x_drop,
                    );
                }
                false
            }
            None => nullb_goes_on_left(score, max_score, x_drop),
        }
    }

    /// The right walk read the letters after its start without a break in `j` steps (it stopped
    /// at the context's end or the subject's), at `score` below its best `max_score`; its next
    /// query byte is at offset `q_at` of the context and the subject's next letter is
    /// `s_seq[s_at]`; `starts_past`: the walk started past the context's letters (the one-hit
    /// extension's walk from a word at the context's end that never scored above 0,
    /// `extend_hit_one_hit`; it took no step). On a batch NCBI's walk goes on into it where the
    /// subject has letters more (`QueryBatch::walk_right`; the batch's last context ends here): its
    /// count at its end, and false. On a context alone: `j`, and whether NCBI's walk could go on,
    /// its next byte (the NULLB after the context) not breaking it (`nullb_goes_on_right`), or
    /// the walk starting past that NULLB (`starts_past`): the search then runs the extension again
    /// on the batch; elsewhere the walk's result is NCBI's.
    #[inline(always)]
    pub(super) fn right_at_bound<S: TbxScore>(
        self,
        sc: S,
        s_seq: &[u8],
        s_at: usize,
        q_at: usize,
        starts_past: bool,
        j: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) -> (usize, bool) {
        match self.batch {
            Some(batch) if s_at < s_seq.len() => (
                batch.walk_right(
                    self.ctx(),
                    self.s_idx(),
                    sc,
                    s_seq,
                    s_at,
                    q_at,
                    j,
                    score,
                    max_score,
                    x_drop,
                ),
                false,
            ),
            Some(_) => (j, false),
            None => (
                j,
                starts_past || nullb_goes_on_right(score, max_score, x_drop),
            ),
        }
    }
}

/// Where a walk of a search's first pass, on a context alone, ended at the context's end and NCBI's
/// could go on (`InBatch::left_at_bound`, `InBatch::right_at_bound`): its best sum and length, its
/// count (the right walk's; 0 for the left walk) and its sum there. The one-hit extension's second
/// pass goes on from there on the query batch (`Pass::Resume`): a walk past a NULLB that scores
/// higher refuses the search, so the best and the length stand, and only the right walk's count
/// can change.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct WalkEnd {
    pub(crate) max_score: i32,
    pub(crate) disp: usize,
    pub(crate) steps: usize,
    pub(crate) score: i32,
}

/// Where a search's first pass stopped: at the end of its left walk, or of its right walk (with
/// the left walk's best and length).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Again {
    Left(WalkEnd),
    Right(WalkEnd, WalkEnd),
}

/// The pass of the one-hit extension (`one_hit::one_hit_walks`): the first, on a context alone,
/// which stops at a walk that asks for the batch; or the second, which goes on from there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pass {
    First,
    Resume(Again),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:904-916 (`s_BlastAaExtendLeft`)
/// Whether the NULLB's step leaves NCBI's left walk going, the walk being at `score` below its
/// best `max_score` without a break (`max_score - score < x_drop` after a step; 0 with no step):
/// the NULLB takes `-NULLB_SCORE` (32,768) off, so the sum stays below the best and the walk breaks
/// unless the X-drop is larger than what it is then below the best. In `i32`:
/// `x_drop - (max_score - score)` lies in `1..=x_drop` after a step and is `x_drop` with none (an
/// X-drop of 0 or less then gives false, as NCBI's step breaks there), so it does not overflow.
#[inline(always)]
fn nullb_goes_on_left(score: i32, max_score: i32, x_drop: i32) -> bool {
    x_drop - (max_score - score) > -NULLB_SCORE
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:847-862 (`s_BlastAaExtendRight`)
/// `nullb_goes_on_left` for the right walk, which also breaks where the sum falls to 0 or less.
#[inline(always)]
fn nullb_goes_on_right(score: i32, max_score: i32, x_drop: i32) -> bool {
    score > -NULLB_SCORE && nullb_goes_on_left(score, max_score, x_drop)
}

#[cfg(test)]
#[path = "query_batch_tests.rs"]
mod tests;
