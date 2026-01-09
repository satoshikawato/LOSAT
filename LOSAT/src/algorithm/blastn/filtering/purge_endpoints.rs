//! HSP endpoint purging for blastn - removes HSPs with common start/end points
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
//! Function: Blast_HSPListPurgeHSPsWithCommonEndpoints

use crate::common::Hit;

/// Purge HSPs with common start or end positions.
///
/// NCBI reference: blast_hits.c:2455-2535 Blast_HSPListPurgeHSPsWithCommonEndpoints
///
/// ```c
/// // NCBI s_QueryOffsetCompareHSPs (blast_hits.c:2268-2320)
/// // Sort by: context, query.offset, subject.offset, score DESC
/// if (h1->context < h2->context) return -1;
/// if (h1->context > h2->context) return 1;
/// if (h1->query.offset < h2->query.offset) return -1;
/// if (h1->query.offset > h2->query.offset) return 1;
/// if (h1->subject.offset < h2->subject.offset) return -1;
/// if (h1->subject.offset > h2->subject.offset) return 1;
/// if (h1->score < h2->score) return 1;  // DESC
/// if (h1->score > h2->score) return -1;
/// ```
///
/// ```c
/// // NCBI purge condition (blast_hits.c:2482-2487)
/// hsp_array[i]->context == hsp_array[i+j]->context &&
/// hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
/// hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
/// hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame
/// ```
///
/// For blastn:
/// - context = query strand (0=plus, 1=minus per query)
/// - subject.frame = subject strand (+1 or -1)
/// Since LOSAT outputs query always forward, we use subject coord ordering for strand.
/// IMPORTANT: Use ACTUAL s_start/s_end, NOT canonical (min/max).
pub fn purge_hsps_with_common_endpoints(mut hits: Vec<Hit>) -> Vec<Hit> {
    if hits.len() <= 1 {
        return hits;
    }

    // Helper to determine strand: minus if s_start > s_end
    // NCBI reference: For blastn, subject.frame encodes strand direction
    let is_minus_strand = |h: &Hit| h.s_start > h.s_end;

    // Pass 1: Remove HSPs with common START positions
    // NCBI reference: blast_hits.c:2478 - s_QueryOffsetCompareHSPs
    // Sort by: (strand, q_start, s_start, score DESC)
    // CRITICAL: Use ACTUAL s_start, not canonical min(s_start, s_end)
    hits.sort_by(|a, b| {
        let a_minus = is_minus_strand(a);
        let b_minus = is_minus_strand(b);

        // NCBI sorts by context (strand), then query.offset, then subject.offset
        a_minus.cmp(&b_minus)
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| a.s_start.cmp(&b.s_start))  // ACTUAL s_start, not canonical
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
            let same_s_start = hits[i].s_start == hits[j].s_start;  // ACTUAL, not canonical

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
    // Sort by: (strand, q_end, s_end, score DESC)
    // CRITICAL: Use ACTUAL s_end, not canonical max(s_start, s_end)
    hits.sort_by(|a, b| {
        let a_minus = is_minus_strand(a);
        let b_minus = is_minus_strand(b);

        a_minus.cmp(&b_minus)
            .then_with(|| a.q_end.cmp(&b.q_end))
            .then_with(|| a.s_end.cmp(&b.s_end))  // ACTUAL s_end, not canonical
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC - higher score first
    });

    // Remove consecutive HSPs with same (strand, q_end, s_end)
    let mut i = 0;
    while i < hits.len() {
        let mut j = i + 1;
        while j < hits.len() {
            let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
            let same_q_end = hits[i].q_end == hits[j].q_end;
            let same_s_end = hits[i].s_end == hits[j].s_end;  // ACTUAL, not canonical

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

