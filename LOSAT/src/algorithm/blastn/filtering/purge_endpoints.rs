//! HSP endpoint purging for blastn - removes HSPs with common start/end points
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
//! Function: Blast_HSPListPurgeHSPsWithCommonEndpoints

use crate::common::Hit;

/// Purge HSPs with common start or end positions.
///
/// NCBI reference: blast_hits.c:2455-2535 Blast_HSPListPurgeHSPsWithCommonEndpoints
///
/// This function removes redundant HSPs that share the same:
/// 1. Start position: (strand, q_start, s_start) - keeps highest scoring
/// 2. End position: (strand, q_end, s_end) - keeps highest scoring
///
/// For blastn, "frame" is represented by strand (plus vs minus).
/// Minus strand is indicated by s_start > s_end in the output.
pub fn purge_hsps_with_common_endpoints(mut hits: Vec<Hit>) -> Vec<Hit> {
    if hits.len() <= 1 {
        return hits;
    }

    // Helper to determine strand: minus if s_start > s_end
    let is_minus_strand = |h: &Hit| h.s_start > h.s_end;

    // Helper to get canonical subject start (smaller value for sorting)
    let subject_start = |h: &Hit| h.s_start.min(h.s_end);

    // Helper to get canonical subject end (larger value for sorting)
    let subject_end = |h: &Hit| h.s_start.max(h.s_end);

    // Pass 1: Remove HSPs with common START positions
    // NCBI reference: blast_hits.c:2478 - s_QueryOffsetCompareHSPs
    // Sort by: (strand, q_start, s_start) with score DESC as tiebreaker
    hits.sort_by(|a, b| {
        let a_minus = is_minus_strand(a);
        let b_minus = is_minus_strand(b);

        a_minus.cmp(&b_minus)
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| subject_start(a).cmp(&subject_start(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC - higher score first
    });

    // Remove consecutive HSPs with same (strand, q_start, s_start)
    // Keep first (highest score) for each unique combination
    let mut i = 0;
    while i < hits.len() {
        let mut j = i + 1;
        while j < hits.len() {
            let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
            let same_q_start = hits[i].q_start == hits[j].q_start;
            let same_s_start = subject_start(&hits[i]) == subject_start(&hits[j]);

            if same_strand && same_q_start && same_s_start {
                hits.remove(j);
            } else {
                break;
            }
        }
        i += 1;
    }

    // Pass 2: Remove HSPs with common END positions
    // NCBI reference: blast_hits.c:2504 - s_QueryEndCompareHSPs
    // Sort by: (strand, q_end, s_end) with score DESC as tiebreaker
    hits.sort_by(|a, b| {
        let a_minus = is_minus_strand(a);
        let b_minus = is_minus_strand(b);

        a_minus.cmp(&b_minus)
            .then_with(|| a.q_end.cmp(&b.q_end))
            .then_with(|| subject_end(a).cmp(&subject_end(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC - higher score first
    });

    // Remove consecutive HSPs with same (strand, q_end, s_end)
    let mut i = 0;
    while i < hits.len() {
        let mut j = i + 1;
        while j < hits.len() {
            let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
            let same_q_end = hits[i].q_end == hits[j].q_end;
            let same_s_end = subject_end(&hits[i]) == subject_end(&hits[j]);

            if same_strand && same_q_end && same_s_end {
                hits.remove(j);
            } else {
                break;
            }
        }
        i += 1;
    }

    hits
}
