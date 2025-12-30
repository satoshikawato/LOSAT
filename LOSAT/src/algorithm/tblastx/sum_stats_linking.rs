//! NCBI-style sum-statistics (even-gap) HSP linking for TBLASTX.
//!
//! This module implements the even-gap linking algorithm from NCBI BLAST's link_hsps.c.
//! HSPs within the same context (query_id, subject_id, q_frame, s_frame) are linked.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c

use rustc_hash::FxHashMap;
use crate::stats::sum_statistics::{
    gap_decay_divisor, small_gap_sum_e, large_gap_sum_e, normalize_score,
    defaults::{GAP_DECAY_RATE_UNGAPPED, GAP_SIZE, OVERLAP_SIZE},
};
use crate::stats::KarlinParams;
use crate::stats::search_space::SearchSpace;

use super::chaining::ExtendedHit;

/// Key for grouping HSPs: (query_id, subject_id, q_frame, s_frame)
type ContextKey = (String, String, i8, i8);

/// NCBI BLAST default window size = gap_size + overlap_size + 1 = 40 + 9 + 1 = 50
const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;

/// Apply NCBI BLAST-style sum-statistics even-gap linking (per-frame).
pub fn apply_sum_stats_even_gap_linking(
    hits: Vec<ExtendedHit>,
    params: &KarlinParams,
) -> Vec<ExtendedHit> {
    if hits.is_empty() {
        return hits;
    }

    // Group HSPs by context (query_id, subject_id, q_frame, s_frame)
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

    let mut result_hits = hits;
    let gap_decay_rate = GAP_DECAY_RATE_UNGAPPED;
    let window_size = WINDOW_SIZE as i64;
    
    for (_key, indices) in groups {
        if indices.len() < 2 {
            continue; // No linking needed for single HSP
        }

        let first_hit = &result_hits[indices[0]];
        let q_aa_len = first_hit.q_orig_len / 3;
        let s_aa_len = first_hit.s_orig_len / 3;
        
        let search_space = SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, params);
        let eff_searchsp = search_space.effective_space as i64;

        // Sort by query position
        let mut sorted_indices = indices.clone();
        sorted_indices.sort_by_key(|&idx| {
            let hit = &result_hits[idx];
            (hit.q_aa_start, hit.s_aa_start)
        });

        // Greedy chaining
        let num_hsps = sorted_indices.len();
        let mut chain_id: Vec<usize> = (0..num_hsps).collect();
        let mut chain_score: Vec<i32> = sorted_indices.iter()
            .map(|&idx| result_hits[idx].raw_score)
            .collect();
        let mut chain_xsum: Vec<f64> = sorted_indices.iter()
            .map(|&idx| normalize_score(result_hits[idx].raw_score, params.lambda, params.k.ln()))
            .collect();
        let mut chain_len: Vec<usize> = vec![1; num_hsps];

        for i in 1..num_hsps {
            let hit_i = &result_hits[sorted_indices[i]];
            let i_q_start = hit_i.q_aa_start as i64;
            let i_s_start = hit_i.s_aa_start as i64;
            
            let mut best_j: Option<usize> = None;
            let mut best_score = 0i32;

            let start_j = if i > 50 { i - 50 } else { 0 };
            for j in (start_j..i).rev() {
                let hit_j = &result_hits[sorted_indices[j]];
                let j_q_end = hit_j.q_aa_end as i64;
                let j_s_end = hit_j.s_aa_end as i64;

                let q_gap = i_q_start - j_q_end;
                let s_gap = i_s_start - j_s_end;

                if q_gap >= 0 && q_gap <= window_size 
                   && s_gap >= 0 && s_gap <= window_size {
                    if chain_score[j] > best_score {
                        best_j = Some(j);
                        best_score = chain_score[j];
                    }
                }
            }

            if let Some(j) = best_j {
                let score_i = result_hits[sorted_indices[i]].raw_score;
                let xsum_i = normalize_score(score_i, params.lambda, params.k.ln());
                
                chain_id[i] = chain_id[j];
                chain_score[i] = chain_score[j] + score_i;
                chain_xsum[i] = chain_xsum[j] + xsum_i;
                chain_len[i] = chain_len[j] + 1;
            }
        }

        // Find final HSP for each chain
        let mut chain_final: FxHashMap<usize, usize> = FxHashMap::default();
        for i in 0..num_hsps {
            chain_final.insert(chain_id[i], i);
        }

        // Assign E-values
        for (&cid, &final_idx) in &chain_final {
            let num_in_chain = chain_len[final_idx];
            let xsum = chain_xsum[final_idx];
            let divisor = gap_decay_divisor(gap_decay_rate, num_in_chain);

            let e_small = if num_in_chain == 1 {
                search_space.effective_space * (-xsum).exp()
            } else {
                small_gap_sum_e(
                    WINDOW_SIZE,
                    num_in_chain as i16,
                    xsum,
                    q_aa_len as i32,
                    s_aa_len as i32,
                    eff_searchsp,
                    divisor,
                )
            };

            let e_large = if num_in_chain == 1 {
                search_space.effective_space * (-xsum).exp()
            } else {
                large_gap_sum_e(
                    num_in_chain as i16,
                    xsum,
                    q_aa_len as i32,
                    s_aa_len as i32,
                    eff_searchsp,
                    divisor,
                )
            };

            let set_evalue = e_small.min(e_large);

            for i in 0..num_hsps {
                if chain_id[i] == cid {
                    result_hits[sorted_indices[i]].hit.e_value = set_evalue;
                }
            }
        }
    }

    result_hits
}
