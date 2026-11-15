//! The query batch around the context whose letters an ungapped extension walks.
//!
//! NCBI's protein walks run on the query batch, not on a context:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:841 (`s_BlastAaExtendRight`)
//! ```c
//!     n = MIN(subject->length - s_off, query->length - q_off);
//! ```
//! and :897 (`s_BlastAaExtendLeft`)
//! ```c
//!     n = MIN(s_off, q_off);
//! ```
//! `query` is the batch's sequence block: `query->sequence` holds the letters of every context of
//! the batch, `q_off` is an offset into it and `query->length` its length (`buflen - 2`,
//! api/blast_setup_cxx.cpp:654, with `QueryInfo_GetSeqBufLen` of blast_query_info.c:246-250: the
//! end of the last context's letters). The contexts lie one after the other with one NULLB between
//! two of them (api/blast_setup_cxx.cpp:65-88: `shift = prev_len ? prev_len + 1 : 0`, so a context
//! without letters takes no room; LOSAT makes no such context, `translation::generate_frames`, and
//! `QueryContext::frame_base` is NCBI's `query_offset`, lookup/backbone.rs). The NULLB's row scores
//! `BLAST_SCORE_MIN` (-32768) against every letter (blast_stat.c:1559-1576, kept for `-`).
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
//! (`InBatch::left_past_start`, `InBatch::right_past_end`). With a raw X-drop of 32,768 or less the
//! NULLB always breaks the walk (`maxscore - score >= 32768`), on the left with the same best and
//! length, on the right also with the same stop: `*s_last_off = s_off + i` with `i` the NULLB's
//! step, the count of the walk that ended at the context's end.
//!
//! A walk past the NULLB that scores no higher than its best leaves the extension's result as it
//! was, except that the right walk's stop (`s_last_off`, and through it the diagonal's `last_hit`
//! for the later hits of that diagonal, aa_ungapped.c:593-597) moves past the NULLB: that is NCBI's
//! output, which LOSAT reproduces (PD-LOSAT-NCBI-DEFECTS, rule 2). A walk that scores higher there
//! moves the HSP's start into the previous context or its end into the next: NCBI keeps the HSP in
//! the seed's context with offsets outside it (`s_AdjustInitialHSPOffsets`, whose
//! `ASSERT(q_start >= 0)` is compiled out of release builds, blast_gapalign.c:2384-2392), and its
//! later steps write query coordinates past the query or drop the alignment from every report.
//! LOSAT refuses such a search (PD-LOSAT-NCBI-DEFECTS, rule 3: `refusal`).

use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::scoring::TbxScore;
use crate::algorithm::tblastx::subjects::XReadFailure;

/// The byte between two contexts of NCBI's query batch (`NULLB`, ncbi_std.h), which is also each
/// context buffer's sentinel in LOSAT.
const NULLB: u8 = 0;

/// The query batch of one subject's search, as its two-hit and one-hit extensions read it: the
/// batch's contexts, and where the search reports an extension that LOSAT refuses.
pub(crate) struct QueryBatch<'a> {
    contexts: &'a [QueryContext],
    failure: &'a XReadFailure,
    s_idx: usize,
}

impl<'a> QueryBatch<'a> {
    /// The batch's `contexts` (in the order of their `frame_base`), for the search of subject
    /// `s_idx`, whose refusal goes to `failure`.
    pub(crate) fn new(
        contexts: &'a [QueryContext],
        failure: &'a XReadFailure,
        s_idx: usize,
    ) -> Self {
        Self {
            contexts,
            failure,
            s_idx,
        }
    }

    /// Context `ctx` of the batch, as an extension in it reads the batch.
    #[inline(always)]
    pub(crate) fn at(&'a self, ctx: usize) -> InBatch<'a> {
        InBatch {
            batch: Some(self),
            ctx,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:897-912
    // ```c
    //     n = MIN(s_off, q_off);
    //     best_i = n + 1;
    //
    //     s = subject->sequence + s_off - n;
    //     q = query->sequence + q_off - n;
    //
    //     for (i = n; i >= 0; i--) {
    //         score += matrix[q[i]][s[i]];
    //
    //         if (score > maxscore) {
    //             maxscore = score;
    //             best_i = i;
    //         }
    // ...
    //         if ((maxscore - score) >= dropoff)
    //             break;
    //     }
    // ```
    // The steps of NCBI's left walk after the context's first letter: `s_next` subject letters are
    // left (the next one is `s_seq[s_next - 1]`), and the batch's letters before the context are its
    // `frame_base`. The sums are `i64` (NCBI's `Int4` sum of many NULLBs would overflow).
    /// The left walk of context `ctx` past its start; refuses the search where it scores higher.
    #[cold]
    #[inline(never)]
    fn walk_left<S: TbxScore>(
        &self,
        ctx: usize,
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
                    self.refuse();
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
                    self.refuse();
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
    // ```c
    //     n = MIN(subject->length - s_off, query->length - q_off);
    // ...
    //     for (i = 0; i < n; i++) {
    //         score += matrix[q[i]][s[i]];
    //
    //         if (score > maxscore) {
    //             maxscore = score;
    //             best_i = i;
    //         }
    // ...
    //         if (score <= 0 || (maxscore - score) >= dropoff)
    //             break;
    //     }
    //
    //     *length = best_i + 1;
    //     *s_last_off = s_off + i;
    // ```
    // The steps of NCBI's right walk after the context's letters: the walk has taken `j` steps, the
    // next subject letter is `s_seq[s_at]` and the next query byte is at offset `q_at` of the
    // context (its letters' count, the NULLB after them; or one more, where the one-hit extension
    // starts its walk past its context's end, `extend_hit_one_hit`), and the batch ends at its last
    // context's letters (`query->length`). Returns NCBI's `i` at the end of the walk.
    /// The right walk of context `ctx` past its end, and its stop; refuses the search where it
    /// scores higher.
    #[cold]
    #[inline(never)]
    fn walk_right<S: TbxScore>(
        &self,
        ctx: usize,
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
                    self.refuse();
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
                    self.refuse();
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
    // ```c
    // static NCBI_INLINE void
    // s_AdjustInitialHSPOffsets(BlastInitHSP* init_hsp, Int4 query_start)
    // {
    //     init_hsp->offsets.qs_offsets.q_off -= query_start;
    //     if (init_hsp->ungapped_data) {
    //         init_hsp->ungapped_data->q_start -= query_start;
    //     }
    //     ASSERT(init_hsp->offsets.qs_offsets.q_off >= 0);
    //     ASSERT(init_hsp->ungapped_data ? init_hsp->ungapped_data->q_start >= 0 : TRUE);
    // }
    // ```
    // An HSP that a walk past a NULLB makes span two contexts keeps the seed's context with offsets
    // outside it (the debug-only assertions above), and NCBI's later steps report query
    // coordinates past the query or no alignment (PD-LOSAT-NCBI-DEFECTS, rule 3).
    /// Refuses the search of this subject (and so the query batch's search).
    fn refuse(&self) {
        self.failure.fail(self.s_idx, refusal());
    }
}

/// The error of a search that `QueryBatch` refuses: an ungapped extension that scores higher past
/// the end of its query frame, in the frame next to it in the query batch.
pub(crate) fn refusal() -> anyhow::Error {
    anyhow::anyhow!(
        "an ungapped extension that passes the end of a query frame into the frame next to it in \
         the query batch and scores higher there (a raw -xdrop_ungap above 32,768, with BLOSUM62 \
         about 15,000 bits or more; NCBI BLAST+ 2.17.0 then keeps an alignment of two frames with \
         query coordinates outside the query's frame and writes them past the query or drops the \
         alignment) is not supported by LOSAT's TBLASTX"
    )
}

/// Where the letters an ungapped extension walks lie: context `ctx` of a query batch
/// (`QueryBatch::at`), or, in the tests, a context alone (`ALONE`: the walks stop at its ends, as
/// NCBI's do when the NULLB breaks them).
#[derive(Clone, Copy)]
pub(crate) struct InBatch<'a> {
    batch: Option<&'a QueryBatch<'a>>,
    ctx: usize,
}

#[cfg(test)]
impl InBatch<'static> {
    /// A context alone (the tests; `x_two_hit_cells` counts the cells of the same walks).
    pub(crate) const ALONE: Self = Self {
        batch: None,
        ctx: 0,
    };
}

impl InBatch<'_> {
    /// The left walk read every letter of the context before its start offset without a break,
    /// and the subject has `s_next` letters more: NCBI's walk goes on into the batch
    /// (`QueryBatch::walk_left`). A context alone, or the batch's first, ends here.
    #[inline(always)]
    pub(super) fn left_past_start<S: TbxScore>(
        self,
        sc: S,
        s_seq: &[u8],
        s_next: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) {
        if let (Some(batch), true) = (self.batch, self.ctx > 0) {
            batch.walk_left(self.ctx, sc, s_seq, s_next, score, max_score, x_drop);
        }
    }

    /// The right walk read every letter of the context after its start offset without a break in
    /// `j` steps (it started at or past the context's end where `j` is 0), its next query byte is at
    /// offset `q_at` of the context and the subject's next letter is `s_seq[s_at]`: NCBI's walk
    /// goes on into the batch (`QueryBatch::walk_right`). Returns the walk's count at its end (`j`
    /// for a context alone or the batch's last).
    #[inline(always)]
    pub(super) fn right_past_end<S: TbxScore>(
        self,
        sc: S,
        s_seq: &[u8],
        s_at: usize,
        q_at: usize,
        j: usize,
        score: i32,
        max_score: i32,
        x_drop: i32,
    ) -> usize {
        match self.batch {
            Some(batch) => {
                batch.walk_right(self.ctx, sc, s_seq, s_at, q_at, j, score, max_score, x_drop)
            }
            None => j,
        }
    }
}

#[cfg(test)]
#[path = "query_batch_tests.rs"]
mod tests;
