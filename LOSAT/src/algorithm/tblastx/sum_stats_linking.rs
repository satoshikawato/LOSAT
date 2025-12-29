//! NCBI-style sum-statistics (even-gap) HSP linking for TBLASTX.
//!
//! This module implements a simplified version of the even-gap linking algorithm 
//! from NCBI BLAST's link_hsps.c. For performance, we use a diagonal-based approach
//! that links nearby HSPs within the same context (query_id, subject_id, q_frame, s_frame).
//!
//! The key insight is that multiple HSPs in a chain have a combined E-value
//! that is better than any individual HSP's E-value.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c

use rustc_hash::FxHashMap;
use crate::stats::sum_statistics::{
    gap_decay_divisor, large_gap_sum_e, normalize_score,
    defaults::{GAP_DECAY_RATE_UNGAPPED, GAP_SIZE, OVERLAP_SIZE},
};
use crate::stats::KarlinParams;
use crate::stats::search_space::SearchSpace;

use super::chaining::ExtendedHit;

/// Key for grouping HSPs: (query_id, subject_id, q_frame, s_frame)
type ContextKey = (String, String, i8, i8);

/// Apply NCBI BLAST-style sum-statistics even-gap linking.
///
/// This implements a simplified greedy version of the algorithm from link_hsps.c:
/// 1. Group HSPs by context (query_id, subject_id, q_frame, s_frame)
/// 2. For each context, sort HSPs by query position
/// 3. Use diagonal-based indexing to efficiently find chainable HSPs
/// 4. Calculate set-level E-values using large-gap sum statistics
/// 5. Assign set-level E-values to all HSPs in each chain
pub fn apply_sum_stats_even_gap_linking(
    hits: Vec<ExtendedHit>,
    params: &KarlinParams,
) -> Vec<ExtendedHit> {
    if hits.is_empty() {
        return hits;
    }

    // Group HSPs by context
    let mut groups: FxHashMap<ContextKey, Vec<usize>> = FxHashMap::default();
    for (idx, hit) in hits.iter().enumerate() {
        let key = (
            hit.hit.query_id.clone(),
            hit.hit.subject_id.clone(),
            hit.q_frame,
            hit.s_frame,
        );
        groups.entry(key).or_default().push(idx);
    }

    // Process each context group
    let mut result_hits = hits;
    let gap_decay_rate = GAP_DECAY_RATE_UNGAPPED;
    // Use a larger window size for more aggressive linking
    // NCBI default: GAP_SIZE=40, OVERLAP_SIZE=9 -> window_size=50
    // We use a larger value to capture more chains
    let window_size = 100i64;
    
    for (_key, indices) in groups {
        if indices.is_empty() {
            continue;
        }

        // Get context-specific search space from the first hit
        let first_hit = &result_hits[indices[0]];
        let q_aa_len = first_hit.q_orig_len / 3;
        let s_aa_len = first_hit.s_orig_len / 3;
        
        // Calculate effective search space
        let search_space = SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, params);
        let eff_searchsp = search_space.effective_space as i64;

        // Sort indices by query amino acid position
        let mut sorted_indices = indices.clone();
        sorted_indices.sort_by_key(|&idx| result_hits[idx].q_aa_start);

        // Diagonal-based linking: use diagonal (q_pos - s_pos) as key
        // HSPs on similar diagonals are candidates for linking
        // Key: diagonal band (diagonal / window_size), Value: (chain_id, last_q_end, last_s_end)
        let mut diag_chains: FxHashMap<i64, Vec<(usize, usize, usize)>> = FxHashMap::default();
        let mut chains: Vec<Vec<usize>> = Vec::new();
        
        for &idx in &sorted_indices {
            let hit = &result_hits[idx];
            let diag = hit.q_aa_start as i64 - hit.s_aa_start as i64;
            let diag_band = diag / window_size;
            
            // Check nearby diagonal bands for chainable HSPs
            let mut best_chain: Option<usize> = None;
            let mut best_chain_last_q = 0usize;
            let mut best_chain_last_s = 0usize;
            
            for d in (diag_band - 1)..=(diag_band + 1) {
                if let Some(chain_entries) = diag_chains.get(&d) {
                    for &(chain_id, last_q_end, last_s_end) in chain_entries.iter().rev() {
                        // Check if this HSP can be linked to the chain
                        let q_gap = hit.q_aa_start as i64 - last_q_end as i64;
                        let s_gap = hit.s_aa_start as i64 - last_s_end as i64;
                        
                        // Both gaps must be positive (no overlap) and within window
                        if q_gap >= 0 && q_gap <= window_size 
                           && s_gap >= 0 && s_gap <= window_size {
                            // Found a valid chain to extend
                            if best_chain.is_none() || last_q_end > best_chain_last_q {
                                best_chain = Some(chain_id);
                                best_chain_last_q = last_q_end;
                                best_chain_last_s = last_s_end;
                            }
                        }
                    }
                }
            }
            
            let chain_id = if let Some(cid) = best_chain {
                // Add to existing chain
                chains[cid].push(idx);
                cid
            } else {
                // Start a new chain
                let cid = chains.len();
                chains.push(vec![idx]);
                cid
            };
            
            // Update diagonal index
            let entry = diag_chains.entry(diag_band).or_default();
            // Remove old entries for this chain and add new one
            entry.retain(|(cid, _, _)| *cid != chain_id);
            entry.push((chain_id, hit.q_aa_end, hit.s_aa_end));
        }

        // Calculate set-level E-values for each chain
        for chain in chains {
            let num_hsps = chain.len();
            
            // Calculate cumulative normalized score (xsum)
            let mut xsum = 0.0;
            for &idx in &chain {
                let hit = &result_hits[idx];
                xsum += normalize_score(hit.raw_score, params.lambda, params.k.ln());
            }
            
            // Calculate set-level E-value
            let set_evalue = if num_hsps == 1 {
                // Single HSP: use simple E-value formula
                search_space.effective_space * (-xsum).exp()
            } else {
                // Multiple HSPs: use large-gap sum statistics
                let divisor = gap_decay_divisor(gap_decay_rate, num_hsps);
                large_gap_sum_e(
                    num_hsps as i16,
                    xsum,
                    q_aa_len as i32,
                    s_aa_len as i32,
                    eff_searchsp,
                    divisor,
                )
            };
            
            // Assign set-level E-value to all HSPs in the chain
            for &idx in &chain {
                result_hits[idx].hit.e_value = set_evalue;
            }
        }
    }

    result_hits
}
