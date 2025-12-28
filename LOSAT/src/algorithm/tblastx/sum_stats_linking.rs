//! NCBI-style sum-statistics (even-gap) HSP linking for TBLASTX.
//!
//! This module applies a BLAST+ compatible *effect*:
//! - HSPs are linked into chains (small-gap / large-gap models)
//! - A *set-level* E-value is computed for each chain using BLAST sum statistics
//! - That set-level E-value is assigned to all HSPs in the chain
//!
//! Primary references (NCBI BLAST+):
//! - `core/link_hsps.c`  (even-gap linking / DP)
//! - `core/blast_stat.c` (BLAST_SmallGapSumE / BLAST_LargeGapSumE / SumP)
//! - `core/blast_parameters.c:CalculateLinkHSPCutoffs` (cutoff_small_gap/big_gap)

use crate::stats::karlin::raw_score_from_bit_score;
use crate::stats::search_space::SearchSpace;
use crate::stats::sum_statistics::{
    gap_decay_divisor, large_gap_sum_e, small_gap_sum_e,
};
use crate::stats::KarlinParams;
use rustc_hash::FxHashMap;

use super::chaining::{ExtendedHit, SequenceKey};
use super::constants::GAP_TRIGGER_BIT_SCORE;

/// BLAST default parameters for linking HSPs (ungapped).
const GAP_PROB: f64 = 0.5; // BLAST_GAP_PROB
const GAP_DECAY_RATE: f64 = 0.5; // BLAST_GAP_DECAY_RATE
const GAP_SIZE: i32 = 40; // BLAST_GAP_SIZE
const OVERLAP_SIZE: i32 = 9; // BLAST_OVERLAP_SIZE

#[derive(Debug, Clone, Copy)]
struct LinkCutoffs {
    cutoff_small_gap: i32,
    cutoff_big_gap: i32,
    gap_prob: f64,
}

fn frame_aa_len(dna_len: usize, frame: i8) -> usize {
    let shift = (frame.unsigned_abs() as usize).saturating_sub(1);
    dna_len.saturating_sub(shift) / 3
}

fn blast_nint(x: f64) -> i64 {
    // NCBI BLAST_Nint: x += (x >= 0 ? 0.5 : -0.5); return (long)x;
    if x >= 0.0 {
        (x + 0.5) as i64
    } else {
        (x - 0.5) as i64
    }
}

/// Port of NCBI `CalculateLinkHSPCutoffs` for the ungapped case (tblastx).
fn calculate_link_hsp_cutoffs(
    params: &KarlinParams,
    query_len_aa: i32,
    subject_len_aa: i32,
    db_len_aa: i64, // 0 for one-subject searches
    word_cutoff_score_min: i32,
) -> LinkCutoffs {
    let window_size = GAP_SIZE + OVERLAP_SIZE + 1;
    let mut gap_prob = GAP_PROB;
    let gap_decay_rate = GAP_DECAY_RATE;
    let eps = 1.0e-9;

    // Subtract off the expected score.
    let mut q_len = query_len_aa.max(1) as i64;
    let mut s_len = subject_len_aa.max(1) as i64;

    let k = params.k;
    let lambda = params.lambda;
    let h = params.h;

    // expected_length = Nint(log(K*m*n)/H)
    let expected_length = {
        let mn = (q_len as f64) * (s_len as f64);
        if mn <= 0.0 || k <= 0.0 || h <= 0.0 {
            0_i64
        } else {
            blast_nint(((k * mn).ln()) / h)
        }
    };

    q_len = (q_len - expected_length).max(1);
    s_len = (s_len - expected_length).max(1);

    // y_variable depends on whether this is a DB search (db_len_aa > subject_len_aa)
    let y_variable = if db_len_aa > s_len {
        ((db_len_aa as f64) / (s_len as f64)).ln() * k / gap_decay_rate
    } else {
        // subject_length + expected_length is the original subject length (before subtract)
        (((s_len + expected_length) as f64) / (s_len as f64)).ln() * k / gap_decay_rate
    };

    let search_sp: i64 = q_len.saturating_mul(s_len);
    let mut x_variable = 0.25 * y_variable * (search_sp as f64);

    let mut cutoff_big_gap: i32;
    let mut cutoff_small_gap: i32;

    // If sequences are large compared to window_size, allow both models and adjust cutoffs.
    if (search_sp as f64) > 8.0 * (window_size as f64) * (window_size as f64) {
        // big gaps
        x_variable /= 1.0 - gap_prob + eps;
        cutoff_big_gap = (x_variable.ln() / lambda).floor() as i32 + 1;

        // small gaps
        x_variable = y_variable * (window_size as f64) * (window_size as f64);
        x_variable /= gap_prob + eps;
        cutoff_small_gap = (x_variable.ln() / lambda).floor() as i32 + 1;
        cutoff_small_gap = cutoff_small_gap.max(word_cutoff_score_min);
    } else {
        cutoff_big_gap = (x_variable.ln() / lambda).floor() as i32 + 1;
        // Force small-gap rule to be ignored
        gap_prob = 0.0;
        cutoff_small_gap = 0;
    }

    // scale_factor is 1.0 for non-RPS searches
    LinkCutoffs {
        cutoff_small_gap,
        cutoff_big_gap,
        gap_prob,
    }
}

fn set_level_evalue_small_gap(
    window_size: i32,
    num: i16,
    xsum: f64,
    query_eff_len: i32,
    subject_eff_len: i32,
    searchsp_eff: i64,
    gap_prob: f64,
) -> f64 {
    let weight_div = gap_decay_divisor(GAP_DECAY_RATE, num as usize);
    let mut e = small_gap_sum_e(
        window_size,
        num,
        xsum,
        query_eff_len,
        subject_eff_len,
        searchsp_eff,
        weight_div,
    );
    if num > 1 {
        if gap_prob == 0.0 {
            e = i32::MAX as f64;
        } else {
            e /= gap_prob;
            if e > i32::MAX as f64 {
                e = i32::MAX as f64;
            }
        }
    }
    e
}

fn set_level_evalue_large_gap(
    num: i16,
    xsum: f64,
    query_eff_len: i32,
    subject_eff_len: i32,
    searchsp_eff: i64,
    gap_prob: f64,
) -> f64 {
    let weight_div = gap_decay_divisor(GAP_DECAY_RATE, num as usize);
    let mut e = large_gap_sum_e(
        num,
        xsum,
        query_eff_len,
        subject_eff_len,
        searchsp_eff,
        weight_div,
    );
    if num > 1 {
        let denom = 1.0 - gap_prob;
        if denom == 0.0 {
            e = i32::MAX as f64;
        } else {
            e /= denom;
            if e > i32::MAX as f64 {
                e = i32::MAX as f64;
            }
        }
    }
    e
}

/// Apply sum-statistics (even-gap) linking and assign set-level E-values to all hits.
///
/// This consumes `hits` and returns them with `hit.e_value` overwritten by set-level E-values.
pub fn apply_sum_stats_even_gap_linking(
    hits: Vec<ExtendedHit>,
    params: &KarlinParams,
) -> Vec<ExtendedHit> {
    if hits.is_empty() {
        return hits;
    }

    // Group by (query_id, subject_id, q_frame, s_frame) (frame-separated, LOSAT-style).
    let mut groups: FxHashMap<SequenceKey, Vec<ExtendedHit>> = FxHashMap::default();
    for h in hits {
        let key: SequenceKey = (
            h.hit.query_id.clone(),
            h.hit.subject_id.clone(),
            h.q_frame,
            h.s_frame,
        );
        groups.entry(key).or_default().push(h);
    }

    let mut out: Vec<ExtendedHit> = Vec::new();
    for (_key, mut hsps) in groups {
        if hsps.is_empty() {
            continue;
        }

        // Effective lengths / search space (NCBI-style length adjustment).
        let q_len_aa = frame_aa_len(hsps[0].q_orig_len, hsps[0].q_frame) as i32;
        let s_len_aa = frame_aa_len(hsps[0].s_orig_len, hsps[0].s_frame) as i32;
        let ss = SearchSpace::with_length_adjustment(q_len_aa as usize, s_len_aa as usize, params);
        let length_adj = ss.length_adjustment;
        let query_eff_len = ((q_len_aa as i64) - length_adj).max(1) as i32;
        let subject_eff_len = ((s_len_aa as i64) - length_adj).max(1) as i32;
        let searchsp_eff = (query_eff_len as i64).saturating_mul(subject_eff_len as i64);

        // Compute link cutoffs (NCBI CalculateLinkHSPCutoffs analogue).
        // word_cutoff_score_min is bounded by the gap-trigger raw score in NCBI.
        let word_cutoff_score_min = raw_score_from_bit_score(GAP_TRIGGER_BIT_SCORE, params);
        let cutoffs =
            calculate_link_hsp_cutoffs(params, q_len_aa, s_len_aa, 0, word_cutoff_score_min);
        let ignore_small_gaps = cutoffs.cutoff_small_gap == 0;

        // Sort by (reverse) position like NCBI s_RevCompareHSPsTbx within a frame.
        hsps.sort_by(|a, b| {
            b.q_aa_start
                .cmp(&a.q_aa_start)
                .then_with(|| b.q_aa_end.cmp(&a.q_aa_end))
                .then_with(|| b.s_aa_start.cmp(&a.s_aa_start))
                .then_with(|| b.s_aa_end.cmp(&a.s_aa_end))
        });

        let n = hsps.len();
        let window_size = GAP_SIZE + OVERLAP_SIZE + 1;
        let trim_size: i32 = (OVERLAP_SIZE + 1) / 2;

        let lambda = params.lambda;
        let log_k = params.k.ln();

        // Precompute trimmed coords
        let mut q_off_trim = vec![0_i32; n];
        let mut q_end_trim = vec![0_i32; n];
        let mut s_off_trim = vec![0_i32; n];
        let mut s_end_trim = vec![0_i32; n];
        for (i, h) in hsps.iter().enumerate() {
            let q_len = (h.q_aa_end as i32) - (h.q_aa_start as i32);
            let s_len = (h.s_aa_end as i32) - (h.s_aa_start as i32);
            let q_trim = q_len.min(trim_size);
            let s_trim = s_len.min(trim_size);
            q_off_trim[i] = (h.q_aa_start as i32) + q_trim;
            q_end_trim[i] = (h.q_aa_end as i32) - q_trim;
            s_off_trim[i] = (h.s_aa_start as i32) + s_trim;
            s_end_trim[i] = (h.s_aa_end as i32) - s_trim;
        }

        // DP arrays (reused; recomputed after each chain removal).
        let mut link: Vec<[Option<usize>; 2]> = vec![[None, None]; n];
        let mut num: Vec<[i16; 2]> = vec![[0, 0]; n];
        let mut sum: Vec<[i32; 2]> = vec![[0, 0]; n];
        let mut xsum: Vec<[f64; 2]> = vec![[0.0, 0.0]; n];

        let mut active = vec![true; n];
        let mut remaining = n;

        while remaining > 0 {
            // Reset DP state
            for i in 0..n {
                link[i] = [None, None];
                num[i] = [0, 0];
                sum[i] = [0, 0];
                xsum[i] = [0.0, 0.0];
            }

            // Method 0: small gaps
            if !ignore_small_gaps {
                let cutoff0 = cutoffs.cutoff_small_gap;
                for i in 0..n {
                    if !active[i] {
                        continue;
                    }
                    let score_i = hsps[i].raw_score;

                    let mut best_j: Option<usize> = None;
                    let mut best_sum: i32 = 0;
                    let mut best_num: i16 = 0;
                    let mut best_xsum: f64 = 0.0;

                    if score_i > cutoff0 {
                        let h_q_etrim = q_end_trim[i];
                        let h_s_etrim = s_end_trim[i];
                        let h_q_et_gap = h_q_etrim + window_size;
                        let h_s_et_gap = h_s_etrim + window_size;
                        for j in (0..i).rev() {
                            let q_off_t = q_off_trim[j];
                            if q_off_t > h_q_et_gap + trim_size {
                                break;
                            }
                            if !active[j] {
                                continue;
                            }
                            let s_off_t = s_off_trim[j];
                            if q_off_t <= h_q_etrim
                                || s_off_t <= h_s_etrim
                                || q_off_t > h_q_et_gap
                                || s_off_t > h_s_et_gap
                            {
                                continue;
                            }
                            if sum[j][0] > best_sum {
                                best_sum = sum[j][0];
                                best_num = num[j][0];
                                best_xsum = xsum[j][0];
                                best_j = Some(j);
                            }
                        }
                    }

                    link[i][0] = best_j;
                    num[i][0] = best_num + 1;
                    sum[i][0] = best_sum + (score_i - cutoff0);
                    xsum[i][0] = best_xsum + (score_i as f64) * lambda - log_k;
                }
            }

            // Method 1: large gaps
            let cutoff1 = cutoffs.cutoff_big_gap;
            for i in 0..n {
                if !active[i] {
                    continue;
                }
                let score_i = hsps[i].raw_score;

                let mut best_j: Option<usize> = None;
                let mut best_sum: i32 = 0;
                let mut best_num: i16 = 0;
                let mut best_xsum: f64 = 0.0;

                if score_i > cutoff1 {
                    let h_q_etrim = q_end_trim[i];
                    let h_s_etrim = s_end_trim[i];
                    for j in (0..i).rev() {
                        if !active[j] {
                            continue;
                        }
                        let q_off_t = q_off_trim[j];
                        let s_off_t = s_off_trim[j];
                        if q_off_t <= h_q_etrim || s_off_t <= h_s_etrim {
                            continue;
                        }
                        if sum[j][1] > best_sum {
                            best_sum = sum[j][1];
                            best_num = num[j][1];
                            best_xsum = xsum[j][1];
                            best_j = Some(j);
                        }
                    }
                }

                link[i][1] = best_j;
                num[i][1] = best_num + 1;
                sum[i][1] = best_sum + (score_i - cutoff1);
                xsum[i][1] = best_xsum + (score_i as f64) * lambda - log_k;
            }

            // Pick best chain endpoints by sum score (NCBI-style), then decide small vs large by E-value.
            let mut best0: Option<usize> = None;
            let mut best0_sum: i32 = i32::MIN;
            let mut best1: Option<usize> = None;
            let mut best1_sum: i32 = i32::MIN;

            for i in 0..n {
                if !active[i] {
                    continue;
                }
                if !ignore_small_gaps && sum[i][0] >= best0_sum {
                    best0_sum = sum[i][0];
                    best0 = Some(i);
                }
                if sum[i][1] >= best1_sum {
                    best1_sum = sum[i][1];
                    best1 = Some(i);
                }
            }

            let (method, end_idx, prob) = if ignore_small_gaps {
                let end = best1.expect("remaining>0 implies at least one active node");
                let p = set_level_evalue_large_gap(
                    num[end][1],
                    xsum[end][1],
                    query_eff_len,
                    subject_eff_len,
                    searchsp_eff,
                    cutoffs.gap_prob,
                );
                (1usize, end, p)
            } else {
                let end0 = best0.expect("remaining>0 implies at least one active node");
                let end1 = best1.expect("remaining>0 implies at least one active node");
                let p0 = set_level_evalue_small_gap(
                    window_size,
                    num[end0][0],
                    xsum[end0][0],
                    query_eff_len,
                    subject_eff_len,
                    searchsp_eff,
                    cutoffs.gap_prob,
                );
                let p1 = set_level_evalue_large_gap(
                    num[end1][1],
                    xsum[end1][1],
                    query_eff_len,
                    subject_eff_len,
                    searchsp_eff,
                    cutoffs.gap_prob,
                );
                if p0 <= p1 {
                    (0usize, end0, p0)
                } else {
                    (1usize, end1, p1)
                }
            };

            // Remove the chosen chain from the active set, assigning the chain E-value to all members.
            let mut cur = Some(end_idx);
            let mut removed_any = false;
            while let Some(i) = cur {
                if !active[i] {
                    break;
                }
                hsps[i].hit.e_value = prob;
                active[i] = false;
                remaining -= 1;
                removed_any = true;
                cur = link[i][method];
            }

            debug_assert!(removed_any, "must remove at least one node per iteration");
        }

        out.extend(hsps);
    }

    out
}


