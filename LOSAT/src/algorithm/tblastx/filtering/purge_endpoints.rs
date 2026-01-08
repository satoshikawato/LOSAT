//! HSP endpoint purging - removes HSPs with common start/end points
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
//! Function: Blast_HSPListPurgeHSPsWithCommonEndpoints

use super::super::chaining::UngappedHit;

/// NCBI Blast_HSPListPurgeHSPsWithCommonEndpoints equivalent
/// Reference: blast_hits.c:2455-2535
///
/// This function removes HSPs that share common endpoints to reduce redundancy.
/// It performs two passes:
/// 1. Remove HSPs with common start points (same context, query.offset, subject.offset, subject.frame)
/// 2. Remove HSPs with common end points (same context, query.end, subject.end, subject.frame)
///
/// In each pass, the highest-scoring HSP is kept among duplicates.
///
/// NCBI field mapping:
/// - `hsp->context` → `UngappedHit.ctx_idx`
/// - `hsp->query.offset` → `UngappedHit.q_aa_start`
/// - `hsp->query.end` → `UngappedHit.q_aa_end`
/// - `hsp->subject.offset` → `UngappedHit.s_aa_start`
/// - `hsp->subject.end` → `UngappedHit.s_aa_end`
/// - `hsp->subject.frame` → `UngappedHit.s_frame`
/// - `hsp->score` → `UngappedHit.raw_score`
///
/// **Note:** `subject.frame` is NOT in the sort comparators, but IS used
/// in duplicate detection (NCBI lines 2487, 2513).
#[allow(dead_code)]
pub fn purge_hsps_with_common_endpoints(hits: &mut Vec<UngappedHit>) {
    if hits.len() <= 1 {
        return;
    }

    // =========================================================================
    // Phase 1: Remove HSPs with common start points
    // =========================================================================
    // NCBI s_QueryOffsetCompareHSPs (blast_hits.c:2267-2321):
    // Sort order: context ASC → query.offset ASC → subject.offset ASC →
    //             score DESC → query.end ASC → subject.end ASC
    // ```c
    // if (h1->context < h2->context) return -1;
    // if (h1->context > h2->context) return 1;
    // if (h1->query.offset < h2->query.offset) return -1;
    // if (h1->query.offset > h2->query.offset) return 1;
    // if (h1->subject.offset < h2->subject.offset) return -1;
    // if (h1->subject.offset > h2->subject.offset) return 1;
    // if (h1->score < h2->score) return 1;   // DESC
    // if (h1->score > h2->score) return -1;  // DESC
    // if (h1->query.end < h2->query.end) return 1;   // shorter range first
    // if (h1->query.end > h2->query.end) return -1;
    // if (h1->subject.end < h2->subject.end) return 1;
    // if (h1->subject.end > h2->subject.end) return -1;
    // ```
    hits.sort_by(|a, b| {
        a.ctx_idx
            .cmp(&b.ctx_idx)
            .then_with(|| a.q_aa_start.cmp(&b.q_aa_start))
            .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC: higher score first
            .then_with(|| b.q_aa_end.cmp(&a.q_aa_end)) // DESC: larger end first (NCBI line 2310-2313)
            .then_with(|| b.s_aa_end.cmp(&a.s_aa_end)) // DESC: larger end first (NCBI line 2315-2318)
    });

    // NCBI duplicate detection (blast_hits.c:2482-2487):
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i] && hsp_array[i+j] &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
    //        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame)
    // ```
    let mut write_idx = 0;
    let mut i = 0;
    while i < hits.len() {
        // Keep the first HSP of each duplicate group (highest score due to sort)
        hits.swap(write_idx, i);
        let kept_ctx = hits[write_idx].ctx_idx;
        let kept_q_off = hits[write_idx].q_aa_start;
        let kept_s_off = hits[write_idx].s_aa_start;
        let kept_s_frame = hits[write_idx].s_frame;
        write_idx += 1;

        // Skip duplicates with same (context, query.offset, subject.offset, subject.frame)
        let mut j = 1;
        while i + j < hits.len() {
            let h = &hits[i + j];
            if h.ctx_idx == kept_ctx
                && h.q_aa_start == kept_q_off
                && h.s_aa_start == kept_s_off
                && h.s_frame == kept_s_frame
            {
                j += 1;
            } else {
                break;
            }
        }
        i += j;
    }
    hits.truncate(write_idx);

    // =========================================================================
    // Phase 2: Remove HSPs with common end points
    // =========================================================================
    // NCBI s_QueryEndCompareHSPs (blast_hits.c:2332-2387):
    // Sort order: context ASC → query.end ASC → subject.end ASC →
    //             score DESC → query.offset DESC → subject.offset DESC
    // ```c
    // if (h1->context < h2->context) return -1;
    // if (h1->context > h2->context) return 1;
    // if (h1->query.end < h2->query.end) return -1;
    // if (h1->query.end > h2->query.end) return 1;
    // if (h1->subject.end < h2->subject.end) return -1;
    // if (h1->subject.end > h2->subject.end) return 1;
    // if (h1->score < h2->score) return 1;   // DESC
    // if (h1->score > h2->score) return -1;  // DESC
    // if (h1->query.offset < h2->query.offset) return 1;   // shorter range first (larger offset)
    // if (h1->query.offset > h2->query.offset) return -1;
    // if (h1->subject.offset < h2->subject.offset) return 1;
    // if (h1->subject.offset > h2->subject.offset) return -1;
    // ```
    hits.sort_by(|a, b| {
        a.ctx_idx
            .cmp(&b.ctx_idx)
            .then_with(|| a.q_aa_end.cmp(&b.q_aa_end))
            .then_with(|| a.s_aa_end.cmp(&b.s_aa_end))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC: higher score first
            .then_with(|| b.q_aa_start.cmp(&a.q_aa_start)) // DESC: larger offset first (NCBI line 2376-2379)
            .then_with(|| b.s_aa_start.cmp(&a.s_aa_start)) // DESC: larger offset first (NCBI line 2381-2384)
    });

    // NCBI duplicate detection (blast_hits.c:2508-2513):
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i] && hsp_array[i+j] &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
    //        hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame)
    // ```
    write_idx = 0;
    i = 0;
    while i < hits.len() {
        hits.swap(write_idx, i);
        let kept_ctx = hits[write_idx].ctx_idx;
        let kept_q_end = hits[write_idx].q_aa_end;
        let kept_s_end = hits[write_idx].s_aa_end;
        let kept_s_frame = hits[write_idx].s_frame;
        write_idx += 1;

        // Skip duplicates with same (context, query.end, subject.end, subject.frame)
        let mut j = 1;
        while i + j < hits.len() {
            let h = &hits[i + j];
            if h.ctx_idx == kept_ctx
                && h.q_aa_end == kept_q_end
                && h.s_aa_end == kept_s_end
                && h.s_frame == kept_s_frame
            {
                j += 1;
            } else {
                break;
            }
        }
        i += j;
    }
    hits.truncate(write_idx);
}

// =============================================================================
// Tests for purge_hsps_with_common_endpoints (NCBI parity)
// =============================================================================
#[cfg(test)]
mod tests {
    use super::purge_hsps_with_common_endpoints;
    use crate::algorithm::tblastx::chaining::UngappedHit;

    /// Helper to create a minimal UngappedHit for testing purge logic.
    fn make_hit(
        ctx_idx: usize,
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        s_frame: i8,
        raw_score: i32,
    ) -> UngappedHit {
        UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx,
            s_f_idx: 0,
            q_frame: 1,
            s_frame,
            q_aa_start: q_start,
            q_aa_end: q_end,
            s_aa_start: s_start,
            s_aa_end: s_end,
            q_orig_len: 1000,
            s_orig_len: 1000,
            raw_score,
            e_value: 0.0,
            num_ident: 0,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        }
    }

    /// NCBI parity test: equivalent to `testCheckHSPCommonEndpoints` from
    /// `blasthits_unit_test.cpp` lines 1163-1221.
    #[test]
    fn test_purge_ncbi_parity() {
        let scores = [1044, 995, 965, 219, 160, 125, 110, 107, 103];
        let q_offsets = [2, 2, 2, 236, 88, 259, 278, 259, 278];
        let q_ends = [322, 336, 300, 322, 182, 322, 341, 341, 341];
        let s_offsets = [7, 7, 7, 194, 2, 194, 197, 194, 197];
        let s_ends = [292, 293, 301, 292, 96, 292, 260, 260, 266];

        let mut hits: Vec<UngappedHit> = (0..9)
            .map(|i| {
                make_hit(
                    0,
                    q_offsets[i],
                    q_ends[i],
                    s_offsets[i],
                    s_ends[i],
                    0,
                    scores[i] as i32,
                )
            })
            .collect();

        purge_hsps_with_common_endpoints(&mut hits);

        assert_eq!(hits.len(), 3, "Expected 3 surviving HSPs, got {}", hits.len());

        let surviving_indices = [4, 0, 6];
        for &orig_idx in &surviving_indices {
            let found = hits.iter().any(|h| {
                h.q_aa_start == q_offsets[orig_idx]
                    && h.q_aa_end == q_ends[orig_idx]
                    && h.s_aa_start == s_offsets[orig_idx]
                    && h.s_aa_end == s_ends[orig_idx]
                    && h.raw_score == scores[orig_idx] as i32
            });
            assert!(
                found,
                "Expected HSP {} (score={}, q={}-{}, s={}-{}) not found in result",
                orig_idx,
                scores[orig_idx],
                q_offsets[orig_idx],
                q_ends[orig_idx],
                s_offsets[orig_idx],
                s_ends[orig_idx]
            );
        }
    }

    /// Test that subject.frame is used in duplicate detection.
    #[test]
    fn test_purge_different_s_frame_not_duplicate() {
        let mut hits = vec![
            make_hit(0, 10, 50, 20, 60, 1, 100),
            make_hit(0, 10, 50, 20, 60, 2, 90),
        ];

        purge_hsps_with_common_endpoints(&mut hits);

        assert_eq!(hits.len(), 2, "HSPs with different s_frame should both survive");
    }

    /// Test that context (ctx_idx) is used in duplicate detection.
    #[test]
    fn test_purge_different_context_not_duplicate() {
        let mut hits = vec![
            make_hit(0, 10, 50, 20, 60, 1, 100),
            make_hit(1, 10, 50, 20, 60, 1, 90),
        ];

        purge_hsps_with_common_endpoints(&mut hits);

        assert_eq!(hits.len(), 2, "HSPs with different context should both survive");
    }

    /// Regression test documenting why mixed-subject purge is WRONG.
    #[test]
    fn test_mixed_subject_purge_is_wrong() {
        let hit1 = UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx: 0,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: 10,
            q_aa_end: 50,
            s_aa_start: 20,
            s_aa_end: 60,
            q_orig_len: 100,
            s_orig_len: 100,
            raw_score: 100,
            e_value: 0.0,
            num_ident: 0,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        };
        let hit2 = UngappedHit {
            q_idx: 0,
            s_idx: 1,
            ctx_idx: 0,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: 10,
            q_aa_end: 50,
            s_aa_start: 20,
            s_aa_end: 60,
            q_orig_len: 100,
            s_orig_len: 100,
            raw_score: 90,
            e_value: 0.0,
            num_ident: 0,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        };

        let mut subject0_hits = vec![hit1.clone()];
        let mut subject1_hits = vec![hit2.clone()];
        purge_hsps_with_common_endpoints(&mut subject0_hits);
        purge_hsps_with_common_endpoints(&mut subject1_hits);
        let correct_count = subject0_hits.len() + subject1_hits.len();
        assert_eq!(correct_count, 2, "Per-subject purge should preserve both HSPs");

        let mut mixed_hits = vec![hit1, hit2];
        purge_hsps_with_common_endpoints(&mut mixed_hits);
    }
}
