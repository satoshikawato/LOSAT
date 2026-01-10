//! HSP endpoint purging for blastn - removes HSPs with common start/end points
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
//! Function: Blast_HSPListPurgeHSPsWithCommonEndpoints

use crate::common::Hit;
use rustc_hash::FxHashMap;

/// Purge HSPs with common start or end positions.
///
/// NCBI reference: blast_hits.c:2455-2535 Blast_HSPListPurgeHSPsWithCommonEndpoints
///
/// CRITICAL: NCBI processes HSPs per-subject (BlastHSPList is per-subject).
/// This function groups hits by subject_id before purging to match NCBI behavior.
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
/// - context = query_id + query strand (we use query_id as first group key)
/// - subject.frame = subject strand (+1 or -1), encoded by s_start > s_end
/// IMPORTANT: Use ACTUAL s_start/s_end, NOT canonical (min/max).
pub fn purge_hsps_with_common_endpoints(hits: Vec<Hit>) -> Vec<Hit> {
    if hits.len() <= 1 {
        return hits;
    }

    // NCBI reference: blast_hits.c:2455-2535
    // Blast_HSPListPurgeHSPsWithCommonEndpoints operates on BlastHSPList which is per-subject.
    // We must group hits by subject_id and process each group separately.

    // Group hits by subject_id
    let mut subject_groups: FxHashMap<String, Vec<Hit>> = FxHashMap::default();
    for hit in hits {
        subject_groups.entry(hit.subject_id.clone()).or_default().push(hit);
    }

    // Process each subject group independently and collect results
    let mut result: Vec<Hit> = Vec::new();
    for (_subject_id, group_hits) in subject_groups {
        let purged = purge_hsps_for_subject(group_hits);
        result.extend(purged);
    }

    result
}

/// Purge HSPs with common endpoints for a single subject.
///
/// NCBI reference: blast_hits.c:2455-2535
/// This matches NCBI's per-BlastHSPList processing.
fn purge_hsps_for_subject(mut hits: Vec<Hit>) -> Vec<Hit> {
    if hits.len() <= 1 {
        return hits;
    }

    let initial_count = hits.len();
    let mut start_purged = 0usize;
    let mut end_purged = 0usize;

    // Helper to determine strand: minus if s_start > s_end
    // NCBI reference: For blastn, subject.frame encodes strand direction
    let is_minus_strand = |h: &Hit| h.s_start > h.s_end;

    // NCBI uses CANONICAL coordinates: subject.offset < subject.end always
    // ASSERT(hsp->subject.offset < hsp->subject.end) at blast_engine.c:1312
    // We need to use min/max for comparison to match NCBI behavior
    let s_offset = |h: &Hit| h.s_start.min(h.s_end);  // NCBI subject.offset
    let s_end_canon = |h: &Hit| h.s_start.max(h.s_end);  // NCBI subject.end

    // NCBI context for blastn:
    // - context = query_index * NUM_STRANDS + strand_index
    // - NUM_STRANDS = 2 for blastn
    // - strand_index = 0 for plus, 1 for minus
    //
    // Since LOSAT doesn't store context explicitly, we use (query_id, is_minus_strand)
    // to emulate the NCBI context behavior.

    // Pass 1: Remove HSPs with common START positions
    // NCBI reference: blast_hits.c:2478 - s_QueryOffsetCompareHSPs
    // Sort by: (context, query.offset, subject.offset, score DESC)
    // CRITICAL: NCBI does NOT include subject.frame in sort key!
    // This means hits with same (context, q.offset, s.offset) but different frames
    // are interleaved by score. The purge loop exits when frame differs,
    // preventing purge of same-frame hits that are separated by different-frame hits.
    hits.sort_by(|a, b| {
        // NCBI sorts by context (query_id for us), then offsets, then score
        // NO STRAND/FRAME in sort key - this is critical for NCBI parity!
        a.query_id.cmp(&b.query_id)
            .then_with(|| a.q_start.cmp(&b.q_start))
            // NCBI uses subject.offset which is CANONICAL (min)
            .then_with(|| s_offset(a).cmp(&s_offset(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // score DESC
    });

    // Remove consecutive HSPs with same (context, query.offset, subject.offset, subject.frame)
    // NCBI purge condition at blast_hits.c:2482-2487
    let mut i = 0;
    while i < hits.len() {
        let mut j = i + 1;
        while j < hits.len() {
            let same_query = hits[i].query_id == hits[j].query_id;
            let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
            let same_q_start = hits[i].q_start == hits[j].q_start;
            // Use CANONICAL subject.offset for comparison
            let same_s_offset = s_offset(&hits[i]) == s_offset(&hits[j]);

            if same_query && same_strand && same_q_start && same_s_offset {
                hits.remove(j);
                start_purged += 1;
            } else {
                break;
            }
        }
        i += 1;
    }

    // Pass 2: Remove HSPs with common END positions
    // NCBI reference: blast_hits.c:2504 - s_QueryEndCompareHSPs
    // Sort by: (context, query.end, subject.end, score DESC)
    // CRITICAL: NCBI does NOT include subject.frame in sort key (same as pass 1)
    hits.sort_by(|a, b| {
        // NO STRAND/FRAME in sort key - matches NCBI behavior
        a.query_id.cmp(&b.query_id)
            .then_with(|| a.q_end.cmp(&b.q_end))
            // NCBI uses subject.end which is CANONICAL (max)
            .then_with(|| s_end_canon(a).cmp(&s_end_canon(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // score DESC
    });

    // Remove consecutive HSPs with same (context, query.end, subject.end, subject.frame)
    let mut i = 0;
    while i < hits.len() {
        let mut j = i + 1;
        while j < hits.len() {
            let same_query = hits[i].query_id == hits[j].query_id;
            let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
            let same_q_end = hits[i].q_end == hits[j].q_end;
            // Use CANONICAL subject.end for comparison
            let same_s_end = s_end_canon(&hits[i]) == s_end_canon(&hits[j]);

            if same_query && same_strand && same_q_end && same_s_end {
                hits.remove(j);
                end_purged += 1;
            } else {
                break;
            }
        }
        i += 1;
    }

    // Debug output (uncomment if needed for investigation)
    // if initial_count > hits.len() {
    //     eprintln!(
    //         "[DEBUG] purge_hsps_for_subject: {} -> {} (start_purged={}, end_purged={})",
    //         initial_count, hits.len(), start_purged, end_purged
    //     );
    // }
    let _ = (initial_count, start_purged, end_purged); // Silence unused warnings

    hits
}

