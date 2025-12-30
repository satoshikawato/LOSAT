//! NCBI-style sum-statistics HSP linking with lh_helper optimization
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c
//! 
//! Key optimizations from NCBI:
//! 1. `changed` flag: skip recomputation when previous link is still valid
//! 2. `linked_to` counter: track how many HSPs link to this one
//! 3. `path_changed` flag: skip recomputation when no chains were affected
//! 4. `use_current_max`: reuse max chains when they weren't affected by removal
//! 5. `next_larger`: skip HSPs with too-small sum
//! 6. Dual-index approach: index=0 (small gaps), index=1 (large gaps)

use rustc_hash::FxHashMap;
use crate::stats::sum_statistics::{
    gap_decay_divisor, small_gap_sum_e, large_gap_sum_e, normalize_score,
    defaults::{GAP_DECAY_RATE_UNGAPPED, GAP_SIZE, OVERLAP_SIZE, GAP_PROB_UNGAPPED},
};
use crate::stats::KarlinParams;
use crate::stats::search_space::SearchSpace;

use super::chaining::UngappedHit;

type ContextKey = (u32, u32, i8, i8);

const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
const TRIM_SIZE: i32 = (OVERLAP_SIZE + 1) / 2;

/// Ordering method indices (NCBI ELinkOrderingMethod)
const LINK_SMALL_GAPS: usize = 0;
const LINK_LARGE_GAPS: usize = 1;

/// Calculate cutoff_big_gap, cutoff_small_gap and ignore_small_gaps for HSP linking
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
/// CalculateLinkHSPCutoffs()
/// Returns (cutoff[0], cutoff[1], ignore_small_gaps)
fn calculate_link_hsp_cutoffs(
    q_aa_len: usize,
    s_aa_len: usize,
    params: &KarlinParams,
    gap_decay_rate: f64,
) -> ([i32; 2], bool) {
    const EPSILON: f64 = 1.0e-9;
    
    // Calculate expected_length (length adjustment)
    let query_length = q_aa_len as f64;
    let subject_length = s_aa_len as f64;
    let expected_length = (params.k * query_length * subject_length).ln() / params.h;
    
    // Adjusted lengths (NCBI lines 1036-1042)
    let adj_query_len = (query_length - expected_length).max(1.0);
    let adj_subject_len = (subject_length - expected_length).max(1.0);
    
    // y_variable calculation (NCBI lines 1046-1052)
    let y_variable = ((subject_length + expected_length) / subject_length).ln() 
        * params.k / gap_decay_rate;
    
    // search_sp = adjusted query_length * adjusted subject_length
    let search_sp = adj_query_len * adj_subject_len;
    
    // x_variable = 0.25 * y_variable * search_sp (NCBI line 1055)
    let x_variable_base = 0.25 * y_variable * search_sp;
    
    // Check if search space is large enough for small gaps (NCBI line 1062)
    let window_sq = (WINDOW_SIZE * WINDOW_SIZE) as f64;
    let gap_prob = GAP_PROB_UNGAPPED;
    
    let ignore_small_gaps = search_sp <= 8.0 * window_sq;
    
    let cutoffs = if !ignore_small_gaps {
        // Large search space: calculate both cutoffs (NCBI lines 1063-1070)
        let x_big = x_variable_base / (1.0 - gap_prob + EPSILON);
        let cutoff_big = (x_big.ln() / params.lambda).floor() as i32 + 1;
        
        // cutoff_small_gap uses window_size^2 (NCBI lines 1066-1068)
        let x_small = y_variable * window_sq / (gap_prob + EPSILON);
        let cutoff_small = (x_small.ln() / params.lambda).floor() as i32 + 1;
        
        [cutoff_small.max(0), cutoff_big]
    } else {
        // Small search space: only big gap cutoff, small gap disabled (NCBI lines 1072-1078)
        let cutoff_big = (x_variable_base.ln() / params.lambda).floor() as i32 + 1;
        [0, cutoff_big]
    };
    
    (cutoffs, ignore_small_gaps)
}

/// lh_helper structure (NCBI link_hsps.c:100-109)
/// This is the "hot" data used in the inner loop
#[derive(Clone)]
struct LhHelper {
    q_off_trim: i32,
    s_off_trim: i32,
    sum: [i32; 2],      // sum for both indices (small gap, large gap)
    next_larger: usize, // index of next HSP with larger sum[1]
}

/// HSP link information (NCBI BlastHSPLink structure)
/// Reference: link_hsps.c:63-71
#[derive(Clone)]
struct HspLink {
    score: i32,
    q_off_trim: i32,
    s_off_trim: i32,
    q_end_trim: i32,
    s_end_trim: i32,
    sum: [i32; 2],           // sum for both indices
    xsum: [f64; 2],          // normalized sum for both indices
    num: [i16; 2],           // number of HSPs in chain for both indices
    link: [Option<usize>; 2], // Best previous HSP to link with for both indices
    changed: bool,           // Has link been changed? (NCBI: hsp_link.changed)
    linked_to: i32,          // How many HSPs link TO this one (NCBI: linked_to)
}

pub fn apply_sum_stats_even_gap_linking(
    hits: Vec<UngappedHit>,
    params: &KarlinParams,
) -> Vec<UngappedHit> {
    if hits.is_empty() {
        return hits;
    }

    // Strand-based grouping (NCBI style)
    let mut groups: FxHashMap<ContextKey, Vec<usize>> = FxHashMap::default();
    for (idx, hit) in hits.iter().enumerate() {
        let q_strand: i8 = if hit.q_frame > 0 { 1 } else { -1 };
        let s_strand: i8 = if hit.s_frame > 0 { 1 } else { -1 };
        let key = (hit.q_idx, hit.s_idx, q_strand, s_strand);
        groups.entry(key).or_default().push(idx);
    }

    let mut result_hits = hits;
    let gap_decay_rate = GAP_DECAY_RATE_UNGAPPED;
    
    for (_key, indices) in groups {
        if indices.is_empty() { continue; }

        let first_hit = &result_hits[indices[0]];
        let q_aa_len = first_hit.q_orig_len / 3;
        let s_aa_len = first_hit.s_orig_len / 3;
        let search_space = SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, params);
        
        // Calculate cutoffs for both indices (NCBI link_hsps.c:492-495)
        let (cutoffs, ignore_small_gaps) = calculate_link_hsp_cutoffs(q_aa_len, s_aa_len, params, gap_decay_rate);

        // Sort by query position DESCENDING (NCBI s_RevCompareHSPsTbx: link_hsps.c:359-374)
        let mut sorted_indices = indices.clone();
        sorted_indices.sort_by(|&a, &b| {
            let ha = &result_hits[a];
            let hb = &result_hits[b];
            hb.q_aa_start.cmp(&ha.q_aa_start)
                .then(hb.q_aa_end.cmp(&ha.q_aa_end))
                .then(hb.s_aa_start.cmp(&ha.s_aa_start))
                .then(hb.s_aa_end.cmp(&ha.s_aa_end))
        });

        let n = sorted_indices.len();
        
        // Initialize HspLink array with both indices
        let mut hsp_links: Vec<HspLink> = sorted_indices.iter().map(|&idx| {
            let hit = &result_hits[idx];
            let ql = (hit.q_aa_end as i32 - hit.q_aa_start as i32).abs();
            let sl = (hit.s_aa_end as i32 - hit.s_aa_start as i32).abs();
            let qt = TRIM_SIZE.min(ql / 4);
            let st = TRIM_SIZE.min(sl / 4);
            let score = hit.raw_score;
            let xsum_val = normalize_score(score, params.lambda, params.k.ln());
            
            // NCBI initializes sum with (score - cutoff[index]) for each index
            HspLink {
                score,
                q_off_trim: hit.q_aa_start as i32 + qt,
                s_off_trim: hit.s_aa_start as i32 + st,
                q_end_trim: hit.q_aa_end as i32 - qt,
                s_end_trim: hit.s_aa_end as i32 - st,
                sum: [score - cutoffs[0], score - cutoffs[1]],
                xsum: [xsum_val, xsum_val],
                num: [1, 1],
                link: [None, None],
                changed: true,
                linked_to: 0,
            }
        }).collect();

        // lh_helper array
        let mut lh_helpers: Vec<LhHelper> = Vec::with_capacity(n);

        let mut remaining = n;
        let mut first_pass = true;
        let mut path_changed = true;

        while remaining > 0 {
            let mut best: [Option<usize>; 2] = [None, None];
            let mut best_sum: [i32; 2] = [i32::MIN, i32::MIN];
            let mut use_current_max = false;
            
            if !first_pass {
                // Find current max sum among active HSPs for both indices
                if !ignore_small_gaps {
                    for i in 0..n {
                        if hsp_links[i].linked_to != -1000 {
                            if hsp_links[i].sum[0] > best_sum[0] {
                                best_sum[0] = hsp_links[i].sum[0];
                                best[0] = Some(i);
                            }
                        }
                    }
                }
                for i in 0..n {
                    if hsp_links[i].linked_to != -1000 {
                        if hsp_links[i].sum[1] > best_sum[1] {
                            best_sum[1] = hsp_links[i].sum[1];
                            best[1] = Some(i);
                        }
                    }
                }
                
                if !path_changed {
                    use_current_max = true;
                } else {
                    // Check if best chains were affected by removal
                    use_current_max = true;
                    if !ignore_small_gaps {
                        if let Some(bi) = best[0] {
                            let mut cur = bi;
                            loop {
                                if hsp_links[cur].linked_to == -1000 {
                                    use_current_max = false;
                                    break;
                                }
                                match hsp_links[cur].link[0] {
                                    Some(p) => cur = p,
                                    None => break,
                                }
                            }
                        }
                    }
                    if use_current_max {
                        if let Some(bi) = best[1] {
                            let mut cur = bi;
                            loop {
                                if hsp_links[cur].linked_to == -1000 {
                                    use_current_max = false;
                                    break;
                                }
                                match hsp_links[cur].link[1] {
                                    Some(p) => cur = p,
                                    None => break,
                                }
                            }
                        }
                    }
                }
            }

            if !use_current_max {
                // Rebuild lh_helper array
                lh_helpers.clear();
                
                let mut max_sum1: [i32; 3] = [i32::MIN; 3]; // For frame tracking (not used in tblastx)
                
                for i in 0..n {
                    if hsp_links[i].linked_to == -1000 { 
                        lh_helpers.push(LhHelper {
                            q_off_trim: i32::MAX,
                            s_off_trim: i32::MAX,
                            sum: [i32::MIN, i32::MIN],
                            next_larger: 0,
                        });
                        continue; 
                    }
                    
                    // Compute next_larger based on sum[1]
                    let sum1 = hsp_links[i].sum[1];
                    let mut prev = if i > 0 { i - 1 } else { 0 };
                    if i > 0 {
                        while prev > 0 {
                            if hsp_links[prev].linked_to != -1000 {
                                if hsp_links[prev].sum[1] > sum1 {
                                    break;
                                }
                                let nl = if prev < lh_helpers.len() { lh_helpers[prev].next_larger } else { 0 };
                                if nl == 0 || nl >= prev { break; }
                                prev = nl;
                            } else {
                                prev -= 1;
                            }
                        }
                    }
                    let next_larger = if prev > 0 && hsp_links[prev].linked_to != -1000 && hsp_links[prev].sum[1] > sum1 {
                        prev
                    } else {
                        0
                    };
                    
                    // Track max sum for frame (simplified for tblastx)
                    max_sum1[1] = max_sum1[1].max(sum1);
                    
                    lh_helpers.push(LhHelper {
                        q_off_trim: hsp_links[i].q_off_trim,
                        s_off_trim: hsp_links[i].s_off_trim,
                        sum: [hsp_links[i].sum[0], hsp_links[i].sum[1]],
                        next_larger,
                    });
                    
                    hsp_links[i].linked_to = 0;
                }

                best = [None, None];
                best_sum = [i32::MIN, i32::MIN];
                
                // ============ INDEX 0 LOOP (SMALL GAPS) - NCBI lines 691-768 ============
                if !ignore_small_gaps {
                    for i in 0..n {
                        if hsp_links[i].linked_to == -1000 { continue; }
                        
                        let mut h_hsp_num: i16 = 0;
                        let mut h_hsp_sum: i32 = 0;
                        let mut h_hsp_xsum: f64 = 0.0;
                        let mut h_hsp_link: Option<usize> = None;
                        
                        if hsp_links[i].score > cutoffs[0] {
                            let h_qe = hsp_links[i].q_end_trim;
                            let h_se = hsp_links[i].s_end_trim;
                            let h_qe_gap = h_qe + WINDOW_SIZE;
                            let h_se_gap = h_se + WINDOW_SIZE;
                            
                            // Inner loop for small gaps
                            if i > 0 {
                                let mut j = i - 1;
                                loop {
                                    if hsp_links[j].linked_to == -1000 {
                                        if j == 0 { break; }
                                        j -= 1;
                                        continue;
                                    }
                                    
                                    let qo = lh_helpers[j].q_off_trim;
                                    let so = lh_helpers[j].s_off_trim;
                                    let sum = lh_helpers[j].sum[0];
                                    
                                    // Position constraints for small gaps (NCBI lines 720-734)
                                    let b1 = qo <= h_qe;     // overlap in query
                                    let b2 = so <= h_se;     // overlap in subject
                                    let b4 = qo > h_qe_gap;  // beyond window in query
                                    let b5 = so > h_se_gap;  // beyond window in subject
                                    
                                    // Early termination (NCBI line 733)
                                    if qo > h_qe_gap + TRIM_SIZE {
                                        break;
                                    }
                                    
                                    if b1 || b2 || b4 || b5 {
                                        if j == 0 { break; }
                                        j -= 1;
                                        continue;
                                    }
                                    
                                    if sum > h_hsp_sum {
                                        h_hsp_num = hsp_links[j].num[0];
                                        h_hsp_sum = sum;
                                        h_hsp_xsum = hsp_links[j].xsum[0];
                                        h_hsp_link = Some(j);
                                    }
                                    
                                    if j == 0 { break; }
                                    j -= 1;
                                }
                            }
                        }
                        
                        // Update link info for index 0
                        let score = hsp_links[i].score;
                        let new_sum = h_hsp_sum + (score - cutoffs[0]);
                        let new_xsum = h_hsp_xsum + normalize_score(score, params.lambda, params.k.ln());
                        
                        hsp_links[i].sum[0] = new_sum;
                        hsp_links[i].num[0] = h_hsp_num + 1;
                        hsp_links[i].link[0] = h_hsp_link;
                        hsp_links[i].xsum[0] = new_xsum;
                        
                        if i < lh_helpers.len() {
                            lh_helpers[i].sum[0] = new_sum;
                        }
                        
                        if new_sum > best_sum[0] {
                            best_sum[0] = new_sum;
                            best[0] = Some(i);
                        }
                        
                        if let Some(link) = h_hsp_link {
                            hsp_links[link].linked_to += 1;
                        }
                    }
                }
                
                // ============ INDEX 1 LOOP (LARGE GAPS) - NCBI lines 771-896 ============
                for i in 0..n {
                    if hsp_links[i].linked_to == -1000 { continue; }
                    
                    let mut h_hsp_num: i16 = 0;
                    let mut h_hsp_sum: i32 = 0;
                    let mut h_hsp_xsum: f64 = 0.0;
                    let mut h_hsp_link: Option<usize> = None;
                    
                    // Check if previous link is still valid
                    let prev_link = hsp_links[i].link[1];
                    let can_skip = !first_pass && 
                        (prev_link.is_none() || 
                         (prev_link.is_some() && !hsp_links[prev_link.unwrap()].changed));
                    
                    if can_skip {
                        if let Some(pl) = prev_link {
                            h_hsp_num = hsp_links[pl].num[1];
                            h_hsp_sum = hsp_links[pl].sum[1];
                            h_hsp_xsum = hsp_links[pl].xsum[1];
                        }
                        h_hsp_link = prev_link;
                        hsp_links[i].changed = false;
                    } else if hsp_links[i].score > cutoffs[1] {
                        hsp_links[i].changed = true;
                        
                        let h_qe = hsp_links[i].q_end_trim;
                        let h_se = hsp_links[i].s_end_trim;
                        
                        // Initialize with previous best if still active
                        if !first_pass {
                            if let Some(pl) = prev_link {
                                if hsp_links[pl].linked_to >= 0 {
                                    h_hsp_sum = hsp_links[pl].sum[1] - 1;
                                }
                            }
                        }
                        
                        // Inner loop for large gaps (no window constraint)
                        if i > 0 {
                            let mut j = i - 1;
                            loop {
                                if hsp_links[j].linked_to == -1000 {
                                    if j == 0 { break; }
                                    j -= 1;
                                    continue;
                                }
                                
                                let sum = lh_helpers[j].sum[1];
                                let next_larger = lh_helpers[j].next_larger;
                                let qo = lh_helpers[j].q_off_trim;
                                let so = lh_helpers[j].s_off_trim;
                                
                                let b0 = sum <= h_hsp_sum;
                                
                                let next_j = if b0 && next_larger > 0 && next_larger < j {
                                    next_larger
                                } else if j > 0 {
                                    j - 1
                                } else {
                                    0
                                };
                                
                                // Position constraints (overlap only, no window)
                                let b1 = qo <= h_qe;
                                let b2 = so <= h_se;
                                
                                if !(b0 || b1 || b2) {
                                    h_hsp_num = hsp_links[j].num[1];
                                    h_hsp_sum = hsp_links[j].sum[1];
                                    h_hsp_xsum = hsp_links[j].xsum[1];
                                    h_hsp_link = Some(j);
                                }
                                
                                if next_j == 0 && j == 0 { break; }
                                if next_j >= j && j > 0 { j -= 1; } else { j = next_j; }
                                if j == 0 && hsp_links[0].linked_to == -1000 { break; }
                            }
                        }
                    } else {
                        hsp_links[i].changed = true;
                    }
                    
                    // Update link info for index 1
                    let score = hsp_links[i].score;
                    let new_sum = h_hsp_sum + (score - cutoffs[1]);
                    let new_xsum = h_hsp_xsum + normalize_score(score, params.lambda, params.k.ln());
                    
                    hsp_links[i].sum[1] = new_sum;
                    hsp_links[i].num[1] = h_hsp_num + 1;
                    hsp_links[i].link[1] = h_hsp_link;
                    hsp_links[i].xsum[1] = new_xsum;
                    
                    if i < lh_helpers.len() {
                        lh_helpers[i].sum[1] = new_sum;
                        
                        // Update next_larger
                        let mut prev = if i > 0 { i - 1 } else { 0 };
                        while prev > 0 {
                            if hsp_links[prev].linked_to != -1000 {
                                if lh_helpers[prev].sum[1] > new_sum {
                                    break;
                                }
                                let nl = lh_helpers[prev].next_larger;
                                if nl == 0 || nl >= prev { break; }
                                prev = nl;
                            } else {
                                prev -= 1;
                            }
                        }
                        lh_helpers[i].next_larger = if prev > 0 && hsp_links[prev].linked_to != -1000 && lh_helpers[prev].sum[1] > new_sum {
                            prev
                        } else {
                            0
                        };
                    }
                    
                    if new_sum > best_sum[1] {
                        best_sum[1] = new_sum;
                        best[1] = Some(i);
                    }
                    
                    if let Some(link) = h_hsp_link {
                        hsp_links[link].linked_to += 1;
                    }
                }
                
                path_changed = false;
                first_pass = false;
            }
            
            // ============ SELECT BEST ORDERING METHOD (NCBI lines 901-952) ============
            // Calculate E-values for both methods and select the better one
            
            let gap_prob = if ignore_small_gaps { 0.0 } else { GAP_PROB_UNGAPPED };
            let mut prob: [f64; 2] = [f64::MAX, f64::MAX];
            
            // Calculate E-value for small gaps (index 0)
            if !ignore_small_gaps {
                if let Some(bi) = best[0] {
                    // Restore sum by adding back cutoff * num (NCBI line 907-908)
                    let restored_sum = hsp_links[bi].sum[0] + (hsp_links[bi].num[0] as i32) * cutoffs[0];
                    let _ = restored_sum; // Used for debugging/verification
                    
                    let num = hsp_links[bi].num[0] as usize;
                    let xsum = hsp_links[bi].xsum[0];
                    let divisor = gap_decay_divisor(gap_decay_rate, num);
                    
                    prob[0] = small_gap_sum_e(WINDOW_SIZE, num as i16, xsum,
                        search_space.effective_query_len as i32,
                        search_space.effective_db_len as i32,
                        search_space.effective_space as i64, divisor);
                    
                    // Adjust for gap probability (NCBI lines 918-922)
                    if num > 1 && gap_prob > 0.0 {
                        prob[0] /= gap_prob;
                    }
                }
            }
            
            // Calculate E-value for large gaps (index 1)
            if let Some(bi) = best[1] {
                // Restore sum (NCBI line 942-943)
                let restored_sum = hsp_links[bi].sum[1] + (hsp_links[bi].num[1] as i32) * cutoffs[1];
                let _ = restored_sum;
                
                let num = hsp_links[bi].num[1] as usize;
                let xsum = hsp_links[bi].xsum[1];
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                
                prob[1] = large_gap_sum_e(num as i16, xsum,
                    search_space.effective_query_len as i32,
                    search_space.effective_db_len as i32,
                    search_space.effective_space as i64, divisor);
                
                // Adjust for gap probability (NCBI lines 931-934)
                if num > 1 && gap_prob > 0.0 && gap_prob < 1.0 {
                    prob[1] /= 1.0 - gap_prob;
                }
            }
            
            // Select the ordering method with better (smaller) E-value (NCBI line 936-937)
            let ordering_method = if !ignore_small_gaps && prob[0] <= prob[1] {
                LINK_SMALL_GAPS
            } else {
                LINK_LARGE_GAPS
            };
            
            let best_i = match best[ordering_method] {
                Some(i) => i,
                None => {
                    // Fallback to other method if available
                    if let Some(i) = best[1 - ordering_method] {
                        i
                    } else {
                        break;
                    }
                }
            };
            
            let evalue = prob[ordering_method];

            // Mark chain as removed and assign E-value
            if hsp_links[best_i].linked_to > 0 {
                path_changed = true;
            }
            
            let mut cur = best_i;
            loop {
                if hsp_links[cur].linked_to > 1 {
                    path_changed = true;
                }
                
                hsp_links[cur].linked_to = -1000;
                hsp_links[cur].changed = true;
                
                let orig_idx = sorted_indices[cur];
                result_hits[orig_idx].e_value = evalue;
                remaining -= 1;
                
                match hsp_links[cur].link[ordering_method] {
                    Some(p) => cur = p,
                    None => break,
                }
            }
        }
    }

    result_hits
}
