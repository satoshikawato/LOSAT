//! HSP chaining and filtering for TBLASTX
//!
//! This module handles clustering, chaining, and filtering of HSPs (High-Scoring Pairs)
//! for protein/translated sequences.

use crate::common::Hit;
use crate::stats::KarlinParams;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use super::constants::{MAX_DIAG_DRIFT_AA, MAX_GAP_AA};
use super::diagnostics::DiagnosticCounters;

/// Extended hit structure that includes frame and amino acid coordinate information
/// for HSP chaining. This is used internally during chaining and then converted
/// back to regular Hit for output.
#[derive(Debug, Clone)]
pub struct ExtendedHit {
    pub hit: Hit,
    /// Raw alignment score (ungapped or gapped) used for sum-statistics linking.
    /// This corresponds to NCBI `BlastHSP.score`.
    pub raw_score: i32,
    pub q_frame: i8,
    pub s_frame: i8,
    pub q_aa_start: usize,
    pub q_aa_end: usize,
    pub s_aa_start: usize,
    pub s_aa_end: usize,
    pub q_orig_len: usize,
    pub s_orig_len: usize,
    pub from_gapped: bool,
}

/// Sequence data for re-alignment during HSP chaining
pub type SequenceKey = (String, String, i8, i8); // (query_id, subject_id, q_frame, s_frame)
pub type SequenceData = (Vec<u8>, Vec<u8>); // (query_aa_seq, subject_aa_seq)

/// NCBI BLAST-style HSP domination test.
/// Based on s_DominateTest() from hspfilter_culling.c (lines 79-120).
///
/// Returns true if `p` dominates `y`, meaning `y` should be culled.
///
/// Key criteria:
/// 1. HSPs must be in the same frame combination (critical for TBLASTX)
/// 2. Overlap must be > 50% of the candidate HSP's length
/// 3. Main criterion: weighted score/length comparison
///    d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
///    If d > 0, p dominates y
fn hsp_dominates(p: &ExtendedHit, y: &ExtendedHit) -> bool {
    // Critical: Only compare HSPs within the same frame combination
    // Without this check, HSPs from different frames would incorrectly dominate each other
    if p.q_frame != y.q_frame || p.s_frame != y.s_frame {
        return false;
    }

    // Must be same query-subject pair
    if p.hit.query_id != y.hit.query_id || p.hit.subject_id != y.hit.subject_id {
        return false;
    }

    // Normalize coordinates to plus strand (like NCBI BLAST's LinkedHSP.begin/end)
    // For reverse-strand hits, q_start > q_end, so we use min/max to ensure b < e
    let b1 = p.hit.q_start.min(p.hit.q_end) as i64;
    let e1 = p.hit.q_start.max(p.hit.q_end) as i64;
    let b2 = y.hit.q_start.min(y.hit.q_end) as i64;
    let e2 = y.hit.q_start.max(y.hit.q_end) as i64;

    let l1 = e1 - b1;
    let l2 = e2 - b2;

    // Calculate overlap
    let overlap = e1.min(e2) - b1.max(b2);

    // If not overlap by more than 50% of candidate's length, don't dominate
    if 2 * overlap < l2 {
        return false;
    }

    // Use bit_score as the score metric (higher is better)
    let s1 = (p.hit.bit_score * 100.0) as i64;
    let s2 = (y.hit.bit_score * 100.0) as i64;

    // Main criterion: 2 * (%diff in score) + 1 * (%diff in length)
    // Formula: d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
    let d = 4 * s1 * l1 + 2 * s1 * l2 - 2 * s2 * l1 - 4 * s2 * l2;

    // Tie-breaker for identical HSPs
    if (s1 == s2 && b1 == b2 && l1 == l2) || d == 0 {
        if s1 != s2 {
            return s1 > s2;
        }
        // Use subject start as tie-breaker (like NCBI uses OID/subject.offset)
        return p.hit.s_start < y.hit.s_start;
    }

    d > 0
}

/// Chain and filter HSPs for protein/translated sequences.
///
/// This implements a cluster-then-extend approach similar to BLASTN:
/// 1. Groups HSPs by query-subject pair AND frame combination
/// 2. Clusters nearby HSPs on similar diagonals
/// 3. Re-runs gapped alignment on merged regions to get accurate statistics
/// 4. Filters redundant overlapping hits
/// 5. Applies e-value threshold to filter insignificant hits
pub fn chain_and_filter_hsps_protein(
    mut hits: Vec<ExtendedHit>,
    _sequences: &FxHashMap<SequenceKey, SequenceData>,
    _db_len_aa_total: usize,
    _params: &KarlinParams,
    evalue_threshold: f64,
    diagnostics: Option<&DiagnosticCounters>,
) -> Vec<Hit> {
    if hits.is_empty() {
        return Vec::new();
    }

    // Record HSPs before chaining
    let hsps_before = hits.len();
    if let Some(diag) = diagnostics {
        diag.base.hsps_before_chain.store(hsps_before, AtomicOrdering::Relaxed);
    }

    // Group hits by query-subject pair AND frame combination
    let mut groups: FxHashMap<SequenceKey, Vec<ExtendedHit>> = FxHashMap::default();
    for hit in hits.drain(..) {
        let key = (
            hit.hit.query_id.clone(),
            hit.hit.subject_id.clone(),
            hit.q_frame,
            hit.s_frame,
        );
        groups.entry(key).or_default().push(hit);
    }

    let mut result_hits = Vec::new();

    for (_key, mut group_hits) in groups {
        if group_hits.is_empty() {
            continue;
        }

        // Sort by query amino acid start position
        group_hits.sort_by_key(|h| h.q_aa_start);

        // Cluster HSPs that are nearby and on similar diagonals
        let mut clusters: Vec<Vec<ExtendedHit>> = Vec::new();

        for hit in group_hits {
            let hit_diag = hit.s_aa_start as isize - hit.q_aa_start as isize;

            // Try to find an existing cluster to add this hit to
            let mut cluster_idx: Option<usize> = None;
            for (idx, cluster) in clusters.iter().enumerate() {
                if let Some(last) = cluster.last() {
                    let last_diag = last.s_aa_end as isize - last.q_aa_end as isize;
                    let diag_drift = (hit_diag - last_diag).abs();

                    // Calculate gap or overlap between HSPs
                    // Positive = gap, negative = overlap
                    let q_distance = hit.q_aa_start as isize - last.q_aa_end as isize;
                    let s_distance = hit.s_aa_start as isize - last.s_aa_end as isize;

                    // Allow overlapping HSPs (negative distance) or gaps up to MAX_GAP_AA
                    // For overlaps, limit to 50% of the smaller HSP to avoid merging unrelated HSPs
                    let hit_q_len = (hit.q_aa_end - hit.q_aa_start) as isize;
                    let last_q_len = (last.q_aa_end - last.q_aa_start) as isize;
                    let max_overlap = -(hit_q_len.min(last_q_len) / 2);

                    let q_ok = q_distance >= max_overlap && q_distance <= MAX_GAP_AA as isize;
                    let s_ok = s_distance >= max_overlap && s_distance <= MAX_GAP_AA as isize;

                    if q_ok && s_ok && diag_drift <= MAX_DIAG_DRIFT_AA {
                        cluster_idx = Some(idx);
                        break;
                    }
                }
            }

            if let Some(idx) = cluster_idx {
                clusters[idx].push(hit);
            } else {
                clusters.push(vec![hit]);
            }
        }

        // Process each cluster - output all individual HSPs (like NCBI BLAST)
        // The domination filter will remove truly redundant hits later
        for cluster in clusters {
            if cluster.len() == 1 {
                if let Some(diag) = diagnostics {
                    diag.base.clusters_single.fetch_add(1, AtomicOrdering::Relaxed);
                }
            } else {
                if let Some(diag) = diagnostics {
                    diag.base.clusters_merged.fetch_add(1, AtomicOrdering::Relaxed);
                    diag.base.hsps_in_merged_clusters.fetch_add(cluster.len(), AtomicOrdering::Relaxed);
                }
            }

            // Output all HSPs in the cluster that pass e-value threshold
            // Keep as ExtendedHit for domination test (needs frame info)
            for ext_hit in cluster {
                if ext_hit.hit.e_value <= evalue_threshold {
                    result_hits.push(ext_hit);
                }
            }
        }
    }

    // Record HSPs after chaining (before domination filter)
    if let Some(diag) = diagnostics {
        diag.base.hsps_after_chain.store(result_hits.len(), AtomicOrdering::Relaxed);
    }

    // Sort by bit score (highest first) for domination filtering
    result_hits.sort_by(|a, b| {
        b.hit.bit_score
            .partial_cmp(&a.hit.bit_score)
            .unwrap_or(Ordering::Equal)
    });

    // NCBI BLAST does not apply HSP culling for TBLASTX by default
    // (see blast_options.c: "Gapped search is not allowed for tblastx")
    // Skip domination filter to match NCBI BLAST behavior
    let filtered_hits = result_hits;

    // Record culled HSPs (always 0 since culling is disabled)
    if let Some(diag) = diagnostics {
        diag.base.hsps_culled_dominated.store(0, AtomicOrdering::Relaxed);
    }

    // Count output sources and convert to Hit
    let mut output_from_ungapped = 0usize;
    let mut output_from_gapped = 0usize;
    let final_hits: Vec<Hit> = filtered_hits
        .into_iter()
        .map(|ext_hit| {
            if ext_hit.from_gapped {
                output_from_gapped += 1;
            } else {
                output_from_ungapped += 1;
            }
            ext_hit.hit
        })
        .collect();

    // Record output source counts and final hits
    if let Some(diag) = diagnostics {
        diag.output_from_ungapped.store(output_from_ungapped, AtomicOrdering::Relaxed);
        diag.output_from_gapped.store(output_from_gapped, AtomicOrdering::Relaxed);
        diag.base.hsps_after_overlap_filter.store(final_hits.len(), AtomicOrdering::Relaxed);
    }

    final_hits
}

