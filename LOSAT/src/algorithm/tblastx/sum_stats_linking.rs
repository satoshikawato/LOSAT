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
use rayon::prelude::*;
use crate::stats::sum_statistics::{
    gap_decay_divisor, small_gap_sum_e, large_gap_sum_e, normalize_score,
    defaults::{GAP_DECAY_RATE_UNGAPPED, GAP_SIZE, OVERLAP_SIZE, GAP_PROB_UNGAPPED},
};
use crate::stats::KarlinParams;
use crate::stats::search_space::SearchSpace;

use super::chaining::UngappedHit;

// NCBI groups by (query_context, subject_frame_sign) for linking
// For tblastx, query_context is the query frame (1-6), subject_frame_sign is +/-
// We extend this to include q_idx and s_idx for multi-sequence support
type ContextKey = (u32, u32, i8, i8); // (q_idx, s_idx, q_frame, s_frame)

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

    // Group by (q_idx, s_idx, q_strand, s_strand) to match NCBI link_hsps.c
    // NCBI groups by context/3 (query strand) and SIGN(subject.frame) (subject strand)
    // Reference: link_hsps.c lines 524-525
    let mut groups: FxHashMap<ContextKey, Vec<UngappedHit>> = FxHashMap::default();
    for hit in hits {
        let q_strand: i8 = if hit.q_frame > 0 { 1 } else { -1 };
        let s_strand: i8 = if hit.s_frame > 0 { 1 } else { -1 };
        let key = (hit.q_idx, hit.s_idx, q_strand, s_strand);
        groups.entry(key).or_default().push(hit);
    }
    
    // Print group count for debugging
    static PRINTED_GROUPS: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if !PRINTED_GROUPS.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[DEBUG] Number of (q_strand, s_strand) groups: {}", groups.len());
    }

    // Process each group in parallel - all groups use NCBI-compliant linking
    let results: Vec<Vec<UngappedHit>> = groups
        .into_par_iter()
        .map(|(_, group_hits)| {
            link_hsp_group_ncbi(group_hits, params)
        })
        .collect();

    results.into_iter().flatten().collect()
}

/// NCBI-style HSP linking for moderate-sized groups
fn link_hsp_group_ncbi(
    mut group_hits: Vec<UngappedHit>,
    params: &KarlinParams,
) -> Vec<UngappedHit> {
    if group_hits.is_empty() {
        return group_hits;
    }
    
    let n = group_hits.len();
    let gap_decay_rate = GAP_DECAY_RATE_UNGAPPED;
    
    let first_hit = &group_hits[0];
    let q_aa_len = first_hit.q_orig_len / 3;
    let s_aa_len = first_hit.s_orig_len / 3;
    let search_space = SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, params);
    
    // Debug output for first group only (to avoid spam)
    static DEBUG_ONCE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if !DEBUG_ONCE.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[DEBUG] link_hsp_group_ncbi first group:");
        eprintln!("  q_aa_len={}, s_aa_len={}", q_aa_len, s_aa_len);
        eprintln!("  effective_query_len={:.2}, effective_db_len={:.2}", 
            search_space.effective_query_len, search_space.effective_db_len);
        eprintln!("  effective_space={:.2e}, length_adjustment={}", 
            search_space.effective_space, search_space.length_adjustment);
        eprintln!("  lambda={}, k={}, logK={:.4}", params.lambda, params.k, params.k.ln());
    }
    
    let (cutoffs, ignore_small_gaps) = calculate_link_hsp_cutoffs(q_aa_len, s_aa_len, params, gap_decay_rate);
    
    // Debug: print cutoffs for first group
    static DEBUG_CUTOFFS: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if !DEBUG_CUTOFFS.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[DEBUG] Linking cutoffs: small_gap={}, big_gap={}, ignore_small_gaps={}", 
            cutoffs[0], cutoffs[1], ignore_small_gaps);
    }

    // Sort by query position DESCENDING (NCBI s_RevCompareHSPsTbx)
    // Note: NCBI sorts by context and subject.frame first, but we've already grouped by those
    group_hits.sort_by(|a, b| {
        b.q_aa_start.cmp(&a.q_aa_start)
            .then(b.q_aa_end.cmp(&a.q_aa_end))
            .then(b.s_aa_start.cmp(&a.s_aa_start))
            .then(b.s_aa_end.cmp(&a.s_aa_end))
    });

    // Initialize HspLink array (indices 0..n correspond to HSPs)
    // NCBI uses 1-based indexing with lh_helper[0] as sentinel
    // We use 0-based but handle the logic appropriately
    let mut hsp_links: Vec<HspLink> = group_hits.iter().map(|hit| {
        let ql = (hit.q_aa_end as i32 - hit.q_aa_start as i32).abs();
        let sl = (hit.s_aa_end as i32 - hit.s_aa_start as i32).abs();
        let qt = TRIM_SIZE.min(ql / 4);
        let st = TRIM_SIZE.min(sl / 4);
        let score = hit.raw_score;
        let xsum_val = normalize_score(score, params.lambda, params.k.ln());
        
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

    // lh_helper array with sentinel at index 0
    // Size is n+1: [sentinel, hsp0, hsp1, ..., hsp(n-1)]
    let mut lh_helpers: Vec<LhHelper> = Vec::with_capacity(n + 1);
    // Sentinel (NCBI lh_helper[0])
    lh_helpers.push(LhHelper {
        q_off_trim: 0,
        s_off_trim: 0,
        sum: [i32::MIN / 2, i32::MIN / 2],
        next_larger: 0,
    });

    let mut remaining = n;
    let mut first_pass = true;
    // NCBI: path_changed tracks whether any removed HSP had linked_to > 0
    // If path_changed == 0, we can skip chain validation entirely
    let mut path_changed = true;
    let gap_prob = if ignore_small_gaps { 0.0 } else { GAP_PROB_UNGAPPED };

    while remaining > 0 {
        let mut best: [Option<usize>; 2] = [None, None];
        let mut best_sum: [i32; 2] = [-cutoffs[0], -cutoffs[1]];
        let mut use_current_max = false;
        
        // NCBI lines 597-652: Try to reuse previous max
        if !first_pass {
            // Find the current max sums
            if !ignore_small_gaps {
                for i in 0..n {
                    if hsp_links[i].linked_to >= 0 && hsp_links[i].sum[0] >= best_sum[0] {
                        best_sum[0] = hsp_links[i].sum[0];
                        best[0] = Some(i);
                    }
                }
            }
            for i in 0..n {
                if hsp_links[i].linked_to >= 0 && hsp_links[i].sum[1] >= best_sum[1] {
                    best_sum[1] = hsp_links[i].sum[1];
                    best[1] = Some(i);
                }
            }
            
            // NCBI line 635: if(path_changed==0) use_current_max=1
            if !path_changed {
                // No path was changed, use these max sums directly
                use_current_max = true;
            } else {
                // NCBI lines 639-650: If max path hasn't changed, we can use it
                // Walk down best, give up if we find a removed item in path
                use_current_max = true;
                if !ignore_small_gaps {
                    if let Some(bi) = best[0] {
                        let mut cur = bi;
                        while let Some(p) = hsp_links[cur].link[0] {
                            if hsp_links[p].linked_to < 0 { use_current_max = false; break; }
                            cur = p;
                        }
                    }
                }
                if use_current_max {
                    if let Some(bi) = best[1] {
                        let mut cur = bi;
                        while let Some(p) = hsp_links[cur].link[1] {
                            if hsp_links[p].linked_to < 0 { use_current_max = false; break; }
                            cur = p;
                        }
                    }
                }
            }
        }

        if !use_current_max {
            // Rebuild lh_helper array (NCBI lines 626-656)
            lh_helpers.truncate(1); // Keep sentinel
            let mut max_sum1 = i32::MIN / 2;
            
            for i in 0..n {
                if hsp_links[i].linked_to < 0 {
                    lh_helpers.push(LhHelper {
                        q_off_trim: i32::MAX,
                        s_off_trim: i32::MAX,
                        sum: [i32::MIN / 2, i32::MIN / 2],
                        next_larger: 0,
                    });
                } else {
                    let sum1 = hsp_links[i].sum[1];
                    max_sum1 = max_sum1.max(sum1);
                    
                    // Compute next_larger (NCBI lines 647-655)
                    let cur_idx = lh_helpers.len();
                    let mut prev = cur_idx.saturating_sub(1);
                    while prev > 0 && lh_helpers[prev].sum[1] <= sum1 {
                        prev = lh_helpers[prev].next_larger;
                    }
                    
                    lh_helpers.push(LhHelper {
                        q_off_trim: hsp_links[i].q_off_trim,
                        s_off_trim: hsp_links[i].s_off_trim,
                        sum: [hsp_links[i].sum[0], hsp_links[i].sum[1]],
                        next_larger: prev,
                    });
                    hsp_links[i].linked_to = 0;
                }
            }

            best = [None, None];
            best_sum = [-cutoffs[0], -cutoffs[1]];

            // INDEX 0 LOOP (small gaps) - NCBI lines 691-768
            if !ignore_small_gaps {
                for i in 0..n {
                    if hsp_links[i].linked_to < 0 { continue; }
                    
                    let mut h_sum = 0i32;
                    let mut h_xsum = 0.0f64;
                    let mut h_num = 0i16;
                    let mut h_link: Option<usize> = None;
                    
                    if hsp_links[i].score > cutoffs[0] {
                        let h_qe = hsp_links[i].q_end_trim;
                        let h_se = hsp_links[i].s_end_trim;
                        let h_qe_gap = h_qe + WINDOW_SIZE;
                        let h_se_gap = h_se + WINDOW_SIZE;
                        
                        // Inner loop: check previous HSPs (j < i)
                        let h_idx = i + 1; // lh_helpers index (1-based)
                        for j_idx in (1..h_idx).rev() {
                            let qo = lh_helpers[j_idx].q_off_trim;
                            let so = lh_helpers[j_idx].s_off_trim;
                            let sum = lh_helpers[j_idx].sum[0];
                            
                            if qo > h_qe_gap + TRIM_SIZE { break; }
                            if qo <= h_qe || so <= h_se || qo > h_qe_gap || so > h_se_gap { continue; }
                            
                            if sum > h_sum {
                                let j = j_idx - 1; // Convert back to hsp_links index
                                h_num = hsp_links[j].num[0];
                                h_sum = hsp_links[j].sum[0];
                                h_xsum = hsp_links[j].xsum[0];
                                h_link = Some(j);
                            }
                        }
                    }
                    
                    let score = hsp_links[i].score;
                    let new_sum = h_sum + (score - cutoffs[0]);
                    let new_xsum = h_xsum + normalize_score(score, params.lambda, params.k.ln());
                    
                    hsp_links[i].sum[0] = new_sum;
                    hsp_links[i].num[0] = h_num + 1;
                    hsp_links[i].link[0] = h_link;
                    hsp_links[i].xsum[0] = new_xsum;
                    lh_helpers[i + 1].sum[0] = new_sum;
                    
                    if new_sum >= best_sum[0] {
                        best_sum[0] = new_sum;
                        best[0] = Some(i);
                    }
                    if let Some(l) = h_link { hsp_links[l].linked_to += 1; }
                }
            }

            // INDEX 1 LOOP (large gaps) with next_larger optimization - NCBI lines 771-896
            for i in 0..n {
                if hsp_links[i].linked_to < 0 { continue; }
                
                let mut h_sum = 0i32;
                let mut h_xsum = 0.0f64;
                let mut h_num = 0i16;
                let mut h_link: Option<usize> = None;
                
                hsp_links[i].changed = true;
                let prev_link = hsp_links[i].link[1];
                
                // Skip recomputation if previous link unchanged (NCBI lines 786-798)
                let can_skip = !first_pass && 
                    (prev_link.is_none() || 
                     (prev_link.is_some() && !hsp_links[prev_link.unwrap()].changed && hsp_links[prev_link.unwrap()].linked_to >= 0));
                
                if can_skip {
                    if let Some(pl) = prev_link {
                        h_num = hsp_links[pl].num[1];
                        h_sum = hsp_links[pl].sum[1];
                        h_xsum = hsp_links[pl].xsum[1];
                    }
                    h_link = prev_link;
                    hsp_links[i].changed = false;
                } else if hsp_links[i].score > cutoffs[1] {
                    let h_qe = hsp_links[i].q_end_trim;
                    let h_se = hsp_links[i].s_end_trim;
                    
                    // Initialize with previous best if still valid (NCBI lines 812-824)
                    if !first_pass {
                        if let Some(pl) = prev_link {
                            if hsp_links[pl].linked_to >= 0 {
                                h_sum = hsp_links[pl].sum[1] - 1;
                            }
                        }
                    }
                    
                    // Inner loop with next_larger jump optimization (NCBI lines 829-859)
                    let h_idx = i + 1;
                    let mut j_idx = h_idx.saturating_sub(1);
                    while j_idx > 0 {
                        let sum = lh_helpers[j_idx].sum[1];
                        let next_larger = lh_helpers[j_idx].next_larger;
                        let qo = lh_helpers[j_idx].q_off_trim;
                        let so = lh_helpers[j_idx].s_off_trim;
                        
                        let skip = sum <= h_sum;
                        if skip {
                            j_idx = next_larger; // Jump to larger sum
                            continue;
                        }
                        
                        j_idx -= 1;
                        
                        if qo <= h_qe || so <= h_se { continue; }
                        
                        let j = j_idx; // Already decremented, so this is correct hsp index
                        if j < n && hsp_links[j].linked_to >= 0 {
                            h_num = hsp_links[j].num[1];
                            h_sum = hsp_links[j].sum[1];
                            h_xsum = hsp_links[j].xsum[1];
                            h_link = Some(j);
                        }
                    }
                }
                
                let score = hsp_links[i].score;
                let new_sum = h_sum + (score - cutoffs[1]);
                let new_xsum = h_xsum + normalize_score(score, params.lambda, params.k.ln());
                
                hsp_links[i].sum[1] = new_sum;
                hsp_links[i].num[1] = h_num + 1;
                hsp_links[i].link[1] = h_link;
                hsp_links[i].xsum[1] = new_xsum;
                lh_helpers[i + 1].sum[1] = new_sum;
                
                // Update next_larger for this entry (NCBI lines 872-881)
                let cur_idx = i + 1;
                let mut prev = cur_idx.saturating_sub(1);
                while prev > 0 && lh_helpers[prev].sum[1] <= new_sum {
                    prev = lh_helpers[prev].next_larger;
                }
                lh_helpers[cur_idx].next_larger = prev;
                
                if new_sum >= best_sum[1] {
                    best_sum[1] = new_sum;
                    best[1] = Some(i);
                }
                if let Some(l) = h_link { hsp_links[l].linked_to += 1; }
            }
            
            // NCBI line 897: path_changed=0 after first pass computation
            path_changed = false;
            first_pass = false;
        }

        // Select best ordering method (NCBI lines 901-952)
        let mut prob = [f64::MAX, f64::MAX];
        
        if !ignore_small_gaps {
            if let Some(bi) = best[0] {
                let num = hsp_links[bi].num[0] as usize;
                let xsum = hsp_links[bi].xsum[0];
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[0] = small_gap_sum_e(WINDOW_SIZE, num as i16, xsum,
                    search_space.effective_query_len as i32,
                    search_space.effective_db_len as i32,
                    search_space.effective_space as i64, divisor);
                if num > 1 && gap_prob > 0.0 { prob[0] /= gap_prob; }
            }
        }
        
        if let Some(bi) = best[1] {
            let num = hsp_links[bi].num[1] as usize;
            let xsum = hsp_links[bi].xsum[1];
            let divisor = gap_decay_divisor(gap_decay_rate, num);
            prob[1] = large_gap_sum_e(num as i16, xsum,
                search_space.effective_query_len as i32,
                search_space.effective_db_len as i32,
                search_space.effective_space as i64, divisor);
            if num > 1 && gap_prob > 0.0 && gap_prob < 1.0 { prob[1] /= 1.0 - gap_prob; }
            
            // Debug: log the first few E-value calculations
            static DEBUG_EVALUE_COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let count = DEBUG_EVALUE_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if count < 3 {
                eprintln!("[DEBUG] E-value calc #{}: score={}, num={}, xsum={:.4}, divisor={:.4}, e={:.4e}",
                    count, hsp_links[bi].score, num, xsum, divisor, prob[1]);
            }
        }

        let ordering = if !ignore_small_gaps && prob[0] <= prob[1] { 0 } else { 1 };
        let best_i = match best[ordering] {
            Some(i) => i,
            None => match best[1 - ordering] {
                Some(i) => i,
                None => break,
            }
        };
        
        let evalue = prob[ordering];
        
        // Debug: count chain statistics
        static CHAIN_STATS: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if !CHAIN_STATS.swap(true, std::sync::atomic::Ordering::Relaxed) {
            // Count chain size for first group only
            let chain_len = hsp_links[best_i].num[ordering];
            eprintln!("[DEBUG] First chain: len={}, e-value={:.2e}, ordering={}", 
                chain_len, evalue, ordering);
        }

        // NCBI lines 964-968: If best has linked_to > 0, set path_changed
        // This means some other HSP was linking to this chain
        if hsp_links[best_i].linked_to > 0 {
            path_changed = true;
        }

        // Mark chain as removed (NCBI lines 961-1000)
        let mut cur = best_i;
        loop {
            // NCBI line 968: if (H->linked_to>1) path_changed=1
            if hsp_links[cur].linked_to > 1 {
                path_changed = true;
            }
            hsp_links[cur].linked_to = -1000;
            hsp_links[cur].changed = true;
            group_hits[cur].e_value = evalue;
            remaining -= 1;
            
            match hsp_links[cur].link[ordering] {
                Some(p) => cur = p,
                None => break,
            }
        }
    }

    group_hits
}

