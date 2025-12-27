use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::path::PathBuf;
use std::sync::mpsc::channel;
use crate::common::{write_output, Hit};
use crate::config::NuclScoringSpec;
use crate::stats::karlin::{bit_score as calc_bit_score, evalue as calc_evalue};
use crate::stats::search_space::SearchSpace;
use crate::stats::{lookup_nucl_params, KarlinParams};
use crate::utils::dust::{DustMasker, MaskedInterval};
use super::args::BlastnArgs;
use super::alignment::extend_gapped_heuristic;
use super::extension::extend_hit_ungapped;
use super::constants::{X_DROP_GAPPED_FINAL, X_DROP_UNGAPPED, TWO_HIT_WINDOW, MIN_UNGAPPED_SCORE_MEGABLAST, MIN_UNGAPPED_SCORE_BLASTN, MAX_DIRECT_LOOKUP_WORD_SIZE};
use super::lookup::{build_lookup, build_pv_direct_lookup, build_two_stage_lookup, build_direct_lookup, reverse_complement, TwoStageLookup, PvDirectLookup, DirectKmerLookup, KmerLookup, pack_diag_key, is_kmer_masked};
use super::sequence_compare::{compare_28mer_simd, compare_sequences_simd, find_first_mismatch};

pub fn calculate_evalue(
    score: i32,
    q_len: usize,
    db_len: usize,
    db_num_seqs: usize,
    params: &KarlinParams,
) -> (f64, f64) {
    // Use NCBI-compatible length adjustment for database search
    let search_space = SearchSpace::for_database_search(q_len, db_len, db_num_seqs, params, true);
    let bs = calc_bit_score(score, params);
    let ev = calc_evalue(bs, &search_space);
    (bs, ev)
}
pub fn chain_and_filter_hsps(
    mut hits: Vec<Hit>,
    sequences: &FxHashMap<(String, String), (Vec<u8>, Vec<u8>)>,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    db_len_total: usize,
    db_num_seqs: usize,
    params: &KarlinParams,
    _use_dp: bool,
    verbose: bool,
    chain_enabled: bool,
) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    let total_start = std::time::Instant::now();

    // BLAST-compatible mode: when chaining is disabled, skip clustering/merging
    // but still apply overlap filtering to remove redundant hits
    // This matches NCBI BLAST's behavior of outputting individual HSPs without merging
    let mut result_hits: Vec<Hit> = if !chain_enabled {
        if verbose {
            eprintln!(
                "[INFO] Chaining disabled (BLAST-compatible mode): skipping clustering, {} raw HSPs",
                hits.len()
            );
        }
        hits // Use raw hits directly, skip to filtering step below
    } else {
        // === BEGIN CHAINING LOGIC ===

        // Parameters for chaining - adaptive based on HSP quality
    // Key insight: use bit_score/length as a proxy for identity
    // High-quality HSPs (high identity) can be merged across larger gaps
    // Low-quality HSPs (low identity) should use strict gap limits
    const MAX_GAP_STRICT: usize = 50; // For low-quality HSPs
    const MAX_GAP_PERMISSIVE: usize = 5000; // For high-quality HSPs
    const MAX_DIAG_DRIFT_STRICT: isize = 30; // For low-quality HSPs
    const MAX_DIAG_DRIFT_PERMISSIVE: isize = 500; // For high-quality HSPs
    const DIAG_BIN_SIZE: isize = 25; // Bin size for diagonal bucketing optimization

    // Threshold for "high quality" HSP: bit_score / length
    // For blastn with reward=2, penalty=-3:
    //   - 95% identity: expected score ≈ 0.95*2 + 0.05*(-3) = 1.75 per base
    //   - 74% identity: expected score ≈ 0.74*2 + 0.26*(-3) = 0.70 per base
    // Bit score = (raw_score * lambda - ln(K)) / ln(2)
    // For high identity (~95%), bit_score/length ≈ 1.5-1.8
    // For low identity (~74%), bit_score/length ≈ 0.5-0.8
    const HIGH_QUALITY_THRESHOLD: f64 = 1.2; // bit_score/length threshold

    // Group hits by query-subject pair
    let grouping_start = std::time::Instant::now();
    let mut groups: FxHashMap<(String, String), Vec<Hit>> = FxHashMap::default();
    for hit in hits.drain(..) {
        let key = (hit.query_id.clone(), hit.subject_id.clone());
        groups.entry(key).or_default().push(hit);
    }
    let grouping_time = grouping_start.elapsed();

    let mut result_hits = Vec::new();
    let mut total_clustering_time = std::time::Duration::ZERO;
    let mut total_align_time = std::time::Duration::ZERO;
    let mut align_calls = 0usize;
    let mut max_region_len = 0usize;
    let mut total_clusters = 0usize;
    let mut max_cluster_size = 0usize;

    if verbose {
        eprintln!(
            "[INFO] chain_and_filter_hsps: {} groups, grouping took {:.3}s",
            groups.len(),
            grouping_time.as_secs_f64()
        );
    }

    for ((query_id, subject_id), mut group_hits) in groups {
        if group_hits.is_empty() {
            continue;
        }

        let group_size = group_hits.len();
        let clustering_start = std::time::Instant::now();

        // Sort by query start position
        group_hits.sort_by_key(|h| h.q_start);

        // Optimized clustering using diagonal bins and active cluster tracking
        // Key insight: hits are sorted by q_start, so clusters whose last hit's q_end
        // is more than MAX_GAP behind the current hit's q_start can never accept new hits

        // Structure: diag_bin -> list of (cluster_idx, max_q_end)
        // We track cluster quality separately using sum_score/span
        let mut diag_bins: FxHashMap<isize, Vec<(usize, usize)>> = FxHashMap::default();
        let mut clusters: Vec<Vec<Hit>> = Vec::new();
        // Track cluster statistics for quality assessment: (sum_bit_score, q_min, q_max, current_bin)
        let mut cluster_stats: Vec<(f64, usize, usize, isize)> = Vec::new();

        for hit in group_hits {
            let hit_diag = hit.s_start as isize - hit.q_start as isize;
            let hit_diag_bin = hit_diag / DIAG_BIN_SIZE;

            // Determine if this hit is high-quality based on bit_score / length
            let hit_quality = hit.bit_score / (hit.length as f64);
            let hit_is_high_quality = hit_quality >= HIGH_QUALITY_THRESHOLD;

            // Use permissive search range for high-quality HSPs
            let search_max_diag_drift = if hit_is_high_quality {
                MAX_DIAG_DRIFT_PERMISSIVE
            } else {
                MAX_DIAG_DRIFT_STRICT
            };

            // Only search nearby diagonal bins (within search_max_diag_drift)
            let bin_range = (search_max_diag_drift / DIAG_BIN_SIZE) + 1;
            let mut cluster_idx: Option<usize> = None;

            'bin_search: for bin_offset in -bin_range..=bin_range {
                let search_bin = hit_diag_bin + bin_offset;
                if let Some(bin_clusters) = diag_bins.get(&search_bin) {
                    for &(idx, max_q_end) in bin_clusters.iter().rev() {
                        // Calculate cluster quality: average score density = sum_score / span
                        let (sum_score, q_min, q_max, _) = cluster_stats[idx];
                        let cluster_span = (q_max - q_min + 1) as f64;
                        let cluster_density = sum_score / cluster_span;
                        let cluster_is_high_quality = cluster_density >= HIGH_QUALITY_THRESHOLD;

                        // Both the hit and the cluster must be high-quality to use permissive parameters
                        let both_high_quality = hit_is_high_quality && cluster_is_high_quality;
                        let effective_max_gap = if both_high_quality {
                            MAX_GAP_PERMISSIVE
                        } else {
                            MAX_GAP_STRICT
                        };
                        let effective_max_diag = if both_high_quality {
                            MAX_DIAG_DRIFT_PERMISSIVE
                        } else {
                            MAX_DIAG_DRIFT_STRICT
                        };

                        // Skip clusters that are too far behind (can never match)
                        if hit.q_start > max_q_end + effective_max_gap {
                            continue;
                        }

                        let cluster = &clusters[idx];
                        if let Some(last) = cluster.last() {
                            let last_diag = last.s_end as isize - last.q_end as isize;
                            let diag_drift = (hit_diag - last_diag).abs();

                            let q_distance = hit.q_start as isize - last.q_end as isize;
                            let s_distance = hit.s_start as isize - last.s_end as isize;

                            let hit_q_len = (hit.q_end - hit.q_start) as isize;
                            let last_q_len = (last.q_end - last.q_start) as isize;
                            let max_overlap = -(hit_q_len.min(last_q_len) / 2);

                            let q_ok = q_distance >= max_overlap
                                && q_distance <= effective_max_gap as isize;
                            let s_ok = s_distance >= max_overlap
                                && s_distance <= effective_max_gap as isize;

                            if q_ok && s_ok && diag_drift <= effective_max_diag as isize {
                                cluster_idx = Some(idx);
                                break 'bin_search;
                            }
                        }
                    }
                }
            }

            if let Some(idx) = cluster_idx {
                clusters[idx].push(hit.clone());

                // Update cluster stats
                let (sum_score, q_min, q_max, old_bin) = cluster_stats[idx];
                let new_sum_score = sum_score + hit.bit_score;
                let new_q_min = q_min.min(hit.q_start);
                let new_q_max = q_max.max(hit.q_end);
                let new_diag = hit.s_end as isize - hit.q_end as isize;
                let new_bin = new_diag / DIAG_BIN_SIZE;
                cluster_stats[idx] = (new_sum_score, new_q_min, new_q_max, new_bin);

                // Fix: properly re-index cluster when diagonal bin changes
                if new_bin != old_bin {
                    // Remove from old bin
                    if let Some(old_bin_clusters) = diag_bins.get_mut(&old_bin) {
                        old_bin_clusters.retain(|(i, _)| *i != idx);
                    }
                    // Add to new bin
                    diag_bins.entry(new_bin).or_default().push((idx, new_q_max));
                } else {
                    // Just update max_q_end in the same bin
                    if let Some(bin_clusters) = diag_bins.get_mut(&new_bin) {
                        if let Some(entry) = bin_clusters.iter_mut().find(|(i, _)| *i == idx) {
                            entry.1 = new_q_max;
                        }
                    }
                }
            } else {
                let new_idx = clusters.len();
                clusters.push(vec![hit.clone()]);
                cluster_stats.push((hit.bit_score, hit.q_start, hit.q_end, hit_diag_bin));
                diag_bins
                    .entry(hit_diag_bin)
                    .or_default()
                    .push((new_idx, hit.q_end));
            }
        }

        let clustering_time = clustering_start.elapsed();
        total_clustering_time += clustering_time;
        total_clusters += clusters.len();
        max_cluster_size =
            max_cluster_size.max(clusters.iter().map(|c| c.len()).max().unwrap_or(0));

        if verbose && group_size > 1000 {
            eprintln!(
                "[INFO] Group {}/{}: {} hits -> {} clusters in {:.3}s",
                query_id,
                subject_id,
                group_size,
                clusters.len(),
                clustering_time.as_secs_f64()
            );
        }

        // PERFORMANCE OPTIMIZATION: Only run align_region on the best clusters
        // NCBI BLAST is fast because it only does gapped extension on a small subset of candidates
        // We score clusters by their total bit_score and only align the top ones
        const MAX_ALIGN_REGION_CALLS: usize = 20; // Limit expensive align_region calls
        const MIN_CLUSTER_SCORE_FOR_ALIGN: f64 = 500.0; // Minimum score to consider for align_region
        
        // Maximum region size for greedy re-alignment to prevent over-extension
        // NCBI BLAST typically produces alignments up to ~50kbp for blastn task
        // Larger regions should use the best HSP instead of re-alignment to match NCBI behavior
        const MAX_REGION_SIZE_FOR_ALIGN: usize = 50000;

        // Collect multi-HSP clusters with their scores for prioritization
        let mut multi_hsp_clusters_with_scores: Vec<(usize, f64, Vec<Hit>)> = Vec::new();
        let mut single_hsp_hits: Vec<Hit> = Vec::new();

        for (idx, cluster) in clusters.into_iter().enumerate() {
            if cluster.is_empty() {
                continue;
            }

            if cluster.len() == 1 {
                // Single HSP, no re-alignment needed
                single_hsp_hits.push(cluster.into_iter().next().unwrap());
            } else {
                // Multi-HSP cluster - calculate total score
                let (sum_score, _, _, _) = cluster_stats[idx];
                multi_hsp_clusters_with_scores.push((idx, sum_score, cluster));
            }
        }

        // Sort multi-HSP clusters by score (descending) and take top N
        multi_hsp_clusters_with_scores
            .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        // Process each cluster with progress logging
        let num_multi_clusters = multi_hsp_clusters_with_scores.len();
        let mut processed_clusters = 0usize;
        let mut aligned_clusters = 0usize;
        let cluster_process_start = std::time::Instant::now();

        for (_idx, sum_score, cluster) in multi_hsp_clusters_with_scores {
            processed_clusters += 1;

            // Progress logging every 100 clusters for multi-HSP
            if verbose && processed_clusters % 100 == 0 {
                eprintln!("[INFO] Processed {}/{} multi-HSP clusters, {} align_region calls, {:.3}s elapsed",
                          processed_clusters, num_multi_clusters, align_calls,
                          cluster_process_start.elapsed().as_secs_f64());
            }

            // Calculate cluster region size to check if it's too large for re-alignment
            let cluster_q_min = cluster.iter().map(|h| h.q_start).min().unwrap_or(0);
            let cluster_q_max = cluster.iter().map(|h| h.q_end).max().unwrap_or(0);
            let cluster_s_min = cluster.iter().map(|h| h.s_start).min().unwrap_or(0);
            let cluster_s_max = cluster.iter().map(|h| h.s_end).max().unwrap_or(0);
            let cluster_region_size = (cluster_q_max.saturating_sub(cluster_q_min))
                .max(cluster_s_max.saturating_sub(cluster_s_min));

            // Only run align_region for top clusters with high scores AND reasonable region size
            // For lower-scoring clusters or very large regions, just use the best single HSP
            // This prevents over-extension that produces alignments much longer than NCBI BLAST
            let should_align = aligned_clusters < MAX_ALIGN_REGION_CALLS
                && sum_score >= MIN_CLUSTER_SCORE_FOR_ALIGN
                && cluster_region_size <= MAX_REGION_SIZE_FOR_ALIGN;

            if should_align {
                aligned_clusters += 1;

                // NCBI BLAST optimization: Instead of aligning the entire merged region
                // (which can be 657kbp and takes O(n * band_size) time), we use greedy
                // extension from the best HSP's center point. This is O(alignment_length)
                // and much faster for high-identity sequences.
                //
                // Reference: NCBI BLAST's BLAST_GappedAlignmentWithTraceback in blast_gapalign.c
                // extends LEFT and RIGHT from the HSP center using ALIGN_EX or greedy alignment.

                // Get sequences for this query-subject pair
                let key = (query_id.clone(), subject_id.clone());
                if let Some((q_seq, s_seq)) = sequences.get(&key) {
                    // Find the best HSP in the cluster to use as seed for extension
                    let best_hsp = cluster
                        .iter()
                        .max_by(|a, b| {
                            a.bit_score
                                .partial_cmp(&b.bit_score)
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap();

                    // Use the center of the best HSP as the seed point for extension
                    // Convert 1-based coordinates to 0-based
                    let seed_q_start = best_hsp.q_start.saturating_sub(1);
                    let seed_q_end = best_hsp.q_end;
                    let seed_s_start = best_hsp.s_start.saturating_sub(1);
                    let seed_s_end = best_hsp.s_end;

                    // Calculate seed center and length
                    let seed_q_center = (seed_q_start + seed_q_end) / 2;
                    let seed_s_center = (seed_s_start + seed_s_end) / 2;
                    let seed_len = 1; // Start from a single point and extend

                    let align_start = std::time::Instant::now();

                    // ALWAYS use greedy extension in chaining phase (NCBI BLAST approach)
                    // Greedy is better for chaining because:
                    // 1. It doesn't have the MAX_WINDOW_SIZE limit that caps DP at ~10kbp
                    // 2. It's faster and handles high-identity sequences well
                    // 3. NCBI BLAST uses greedy for HSP extension even in blastn task
                    // The use_dp parameter only affects the initial seed extension, not chaining
                    let (
                        final_q_start_0,
                        final_q_end_0,
                        final_s_start_0,
                        final_s_end_0,
                        score,
                        matches,
                        mismatches,
                        gap_opens,
                        gap_letters,
                        _dp_cells,
                    ) = extend_gapped_heuristic(
                        q_seq,
                        s_seq,
                        seed_q_center,
                        seed_s_center,
                        seed_len,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        X_DROP_GAPPED_FINAL,
                        false, // Always use greedy for chaining re-alignment
                    );

                    total_align_time += align_start.elapsed();
                    align_calls += 1;

                    // Track region length for logging
                    let region_len =
                        (final_q_end_0 - final_q_start_0).max(final_s_end_0 - final_s_start_0);
                    max_region_len = max_region_len.max(region_len);

                    // Skip if no valid alignment found
                    if score <= 0 {
                        // Fall back to best HSP
                        result_hits.push(best_hsp.clone());
                        continue;
                    }

                    // Calculate statistics using BLAST's definition:
                    // alignment_length = matches + mismatches + gap_letters (total aligned columns)
                    let aln_len = matches + mismatches + gap_letters;

                    let identity = if aln_len > 0 {
                        ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                    } else {
                        0.0
                    };

                    let (bit_score, e_value) =
                        calculate_evalue(score, q_seq.len(), db_len_total, db_num_seqs, params);

                    // Convert 0-based coordinates to 1-based for output
                    let final_q_start = final_q_start_0 + 1;
                    let final_q_end = final_q_end_0;
                    let final_s_start = final_s_start_0 + 1;
                    let final_s_end = final_s_end_0;

                    result_hits.push(Hit {
                        query_id: query_id.clone(),
                        subject_id: subject_id.clone(),
                        identity,
                        length: aln_len,
                        mismatch: mismatches,
                        gapopen: gap_opens,
                        q_start: final_q_start,
                        q_end: final_q_end,
                        s_start: final_s_start,
                        s_end: final_s_end,
                        e_value,
                        bit_score,
                    });
                } else {
                    // Fallback: just use the best HSP (shouldn't happen)
                    let best_hsp = cluster
                        .into_iter()
                        .max_by(|a, b| {
                            a.bit_score
                                .partial_cmp(&b.bit_score)
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap();
                    result_hits.push(best_hsp);
                }
            } else {
                // PERFORMANCE: For lower-scoring clusters, just use the best single HSP
                // This avoids expensive align_region calls for clusters that won't produce top hits
                let best_hsp = cluster
                    .into_iter()
                    .max_by(|a, b| {
                        a.bit_score
                            .partial_cmp(&b.bit_score)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap();
                result_hits.push(best_hsp);
            }
        }

        // Add single-HSP hits
        result_hits.extend(single_hsp_hits);
    }

        // Log after cluster processing
        if verbose {
            eprintln!("[INFO] Cluster processing done: {} result_hits, {} align_region calls in {:.3}s, max_region={}bp",
                      result_hits.len(), align_calls, total_align_time.as_secs_f64(), max_region_len);
        }

        result_hits // Return chained hits from else block
        // === END CHAINING LOGIC ===
    };

    // Final pass: remove redundant overlapping hits using spatial binning for O(n) instead of O(n²)
    //
    // BLAST-compatible mode (chain_enabled = false):
    //   Uses NCBI BLAST's s_DominateTest algorithm which:
    //   - Compares HSPs based on query coordinates only (not subject)
    //   - Uses a score/length tradeoff formula to determine dominance
    //   - Does NOT use diagonal gating (all HSPs are compared)
    //   Reference: hspfilter_culling.c s_DominateTest()
    //
    // Chaining mode (chain_enabled = true):
    //   Uses diagonal gating to preserve hits on different diagonals
    //   (representing different biological alignments like repeats)
    let filter_start = std::time::Instant::now();
    let result_hits_count = result_hits.len();

    result_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(Ordering::Equal)
    });

    // Use spatial binning to avoid O(n²) comparisons
    // Bin size of 1000bp - hits can only overlap if they share a bin
    const FILTER_BIN_SIZE: usize = 1000;
    const MAX_DIAG_DIFF_FOR_SAME_ALIGNMENT: isize = 100; // Only used when chaining is enabled
    let mut q_bins: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    let mut final_hits: Vec<Hit> = Vec::new();

    for hit in result_hits {
        // Calculate which q bins this hit spans
        let q_bin_start = hit.q_start / FILTER_BIN_SIZE;
        let q_bin_end = hit.q_end / FILTER_BIN_SIZE;

        // Calculate diagonal for this hit (s_start - q_start) - only used when chaining enabled
        let hit_diag = hit.s_start as isize - hit.q_start as isize;

        // Only check hits in overlapping bins
        let mut dominated = false;
        for bin in q_bin_start..=q_bin_end {
            if let Some(kept_indices) = q_bins.get(&bin) {
                for &kept_idx in kept_indices {
                    let kept = &final_hits[kept_idx];

                    if chain_enabled {
                        // Chaining mode: use diagonal gating to preserve different alignments
                        // Hits on different diagonals represent different biological alignments
                        let kept_diag = kept.s_start as isize - kept.q_start as isize;
                        let diag_diff = (hit_diag - kept_diag).abs();
                        if diag_diff > MAX_DIAG_DIFF_FOR_SAME_ALIGNMENT {
                            continue; // Different diagonal = different alignment, don't filter
                        }

                        // Check q overlap (fixed off-by-one: use >= for inclusive ranges)
                        let q_overlap_start = hit.q_start.max(kept.q_start);
                        let q_overlap_end = hit.q_end.min(kept.q_end);
                        let q_overlap = if q_overlap_end >= q_overlap_start {
                            q_overlap_end - q_overlap_start + 1
                        } else {
                            0
                        };
                        let q_hit_len = hit.q_end - hit.q_start + 1;
                        let q_overlap_frac = if q_hit_len > 0 {
                            q_overlap as f64 / q_hit_len as f64
                        } else {
                            0.0
                        };

                        if q_overlap_frac < 0.5 {
                            continue; // No significant q overlap, skip s check
                        }

                        // Check s overlap (fixed off-by-one)
                        let s_overlap_start = hit.s_start.max(kept.s_start);
                        let s_overlap_end = hit.s_end.min(kept.s_end);
                        let s_overlap = if s_overlap_end >= s_overlap_start {
                            s_overlap_end - s_overlap_start + 1
                        } else {
                            0
                        };
                        let s_hit_len = hit.s_end - hit.s_start + 1;
                        let s_overlap_frac = if s_hit_len > 0 {
                            s_overlap as f64 / s_hit_len as f64
                        } else {
                            0.0
                        };

                        if s_overlap_frac >= 0.5 {
                            dominated = true;
                            break;
                        }
                    } else {
                        // BLAST-compatible mode: diagonal-aware overlap filtering
                        // 
                        // NCBI BLAST's culling considers both query AND subject coordinates.
                        // HSPs on different diagonals represent different biological alignments
                        // (e.g., repeats, duplications) and should NOT be filtered even if
                        // they overlap in query space.
                        //
                        // Key insight: For self-comparison, the full sequence match (diagonal 0)
                        // should NOT dominate HSPs on other diagonals representing repeats.
                        //
                        // We use diagonal difference to determine if HSPs are in the same
                        // alignment region. HSPs with large diagonal differences are preserved.

                        // Calculate diagonals (s_start - q_start)
                        let kept_diag = kept.s_start as isize - kept.q_start as isize;
                        let diag_diff = (hit_diag - kept_diag).abs();

                        // For BLAST-compatible mode, use a larger diagonal tolerance
                        // This allows filtering of truly redundant HSPs while preserving
                        // HSPs on different diagonals (repeats, etc.)
                        // 
                        // The key difference from chaining mode:
                        // - Chaining mode uses MAX_DIAG_DIFF_FOR_SAME_ALIGNMENT = 100
                        // - BLAST-compatible mode uses a proportional threshold based on
                        //   the HSP length to handle both short and long alignments
                        let hit_len = (hit.q_end - hit.q_start + 1).max(hit.s_end - hit.s_start + 1);
                        let diag_tolerance = (hit_len / 10).max(50) as isize; // At least 50bp or 10% of length

                        if diag_diff > diag_tolerance {
                            continue; // Different diagonal = different alignment, don't filter
                        }

                        // Check query overlap (fixed off-by-one: use >= for inclusive ranges)
                        let q_overlap_start = hit.q_start.max(kept.q_start);
                        let q_overlap_end = hit.q_end.min(kept.q_end);
                        let q_overlap = if q_overlap_end >= q_overlap_start {
                            q_overlap_end - q_overlap_start + 1
                        } else {
                            0
                        };
                        let q_hit_len = hit.q_end - hit.q_start + 1;
                        let q_overlap_frac = if q_hit_len > 0 {
                            q_overlap as f64 / q_hit_len as f64
                        } else {
                            0.0
                        };

                        // Need >50% query overlap to be considered redundant
                        if q_overlap_frac < 0.5 {
                            continue;
                        }

                        // Check subject overlap as well
                        let s_overlap_start = hit.s_start.max(kept.s_start);
                        let s_overlap_end = hit.s_end.min(kept.s_end);
                        let s_overlap = if s_overlap_end >= s_overlap_start {
                            s_overlap_end - s_overlap_start + 1
                        } else {
                            0
                        };
                        let s_hit_len = hit.s_end - hit.s_start + 1;
                        let s_overlap_frac = if s_hit_len > 0 {
                            s_overlap as f64 / s_hit_len as f64
                        } else {
                            0.0
                        };

                        // Need >50% subject overlap as well
                        if s_overlap_frac < 0.5 {
                            continue;
                        }

                        // Both query and subject have >50% overlap AND similar diagonal
                        // This is a redundant HSP - filter it out
                        dominated = true;
                        break;
                    }
                }
                if dominated {
                    break;
                }
            }
        }

        if !dominated {
            let new_idx = final_hits.len();
            // Register this hit in all bins it spans
            for bin in q_bin_start..=q_bin_end {
                q_bins.entry(bin).or_default().push(new_idx);
            }
            final_hits.push(hit);
        }
    }
    let filter_time = filter_start.elapsed();

    if verbose {
        eprintln!("[INFO] chain_and_filter_hsps summary:");
        eprintln!(
            "[INFO]   Total time: {:.3}s",
            total_start.elapsed().as_secs_f64()
        );
        if chain_enabled {
            eprintln!("[INFO]   Chaining: enabled");
        } else {
            eprintln!("[INFO]   Chaining: disabled (BLAST-compatible mode)");
        }
        eprintln!(
            "[INFO]   Filtering: {:.3}s ({} -> {} hits)",
            filter_time.as_secs_f64(),
            result_hits_count,
            final_hits.len()
        );
    }

    final_hits
}

pub fn run(args: BlastnArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    // Determine effective word size based on task
    // If user specified default word_size (28) and task is not megablast, use task-appropriate defaults
    let effective_word_size = match args.task.as_str() {
        "megablast" => args.word_size,
        "blastn" => {
            if args.word_size == 28 {
                11
            } else {
                args.word_size
            }
        }
        "dc-megablast" => {
            if args.word_size == 28 {
                11
            } else {
                args.word_size
            }
        }
        _ => args.word_size,
    };

    // Determine effective scoring parameters based on task
    // NCBI BLAST uses different defaults for different tasks:
    // - megablast: reward=1, penalty=-2, gapopen=0, gapextend=0 (linear gap)
    // - blastn: reward=2, penalty=-3, gapopen=5, gapextend=2
    // - dc-megablast: reward=2, penalty=-3, gapopen=5, gapextend=2
    // - blastn-short: reward=1, penalty=-3, gapopen=5, gapextend=2
    // Note: gap costs are stored as negative values internally for DP scoring
    let (reward, penalty, gap_open, gap_extend) = match args.task.as_str() {
        "megablast" => {
            // Use user-specified values or megablast defaults
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -2 } else { args.penalty };
            let go = if args.gap_open == 0 { 0 } else { args.gap_open };
            let ge = if args.gap_extend == 0 {
                0
            } else {
                args.gap_extend
            }; // Linear gap (0 triggers non-affine greedy)
            (r, p, go, ge)
        }
        "blastn" | "dc-megablast" => {
            // Use user-specified values or blastn defaults (reward=2, penalty=-3)
            let r = if args.reward == 1 { 2 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            let go = if args.gap_open == 0 {
                -5
            } else {
                args.gap_open
            };
            let ge = if args.gap_extend == 0 {
                -2
            } else {
                args.gap_extend
            };
            (r, p, go, ge)
        }
        "blastn-short" => {
            // blastn-short: reward=1, penalty=-3, gapopen=5, gapextend=2
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            let go = if args.gap_open == 0 {
                -5
            } else {
                args.gap_open
            };
            let ge = if args.gap_extend == 0 {
                -2
            } else {
                args.gap_extend
            };
            (r, p, go, ge)
        }
        _ => (args.reward, args.penalty, args.gap_open, args.gap_extend),
    };

    // Task-specific ungapped score threshold for triggering gapped extension
    // Higher threshold for blastn reduces the number of gapped extensions significantly
    // This is critical for performance on self-comparison with many off-diagonal matches
    let min_ungapped_score = match args.task.as_str() {
        "megablast" => MIN_UNGAPPED_SCORE_MEGABLAST,
        _ => MIN_UNGAPPED_SCORE_BLASTN,
    };

    // NCBI BLAST uses different extension algorithms based on task:
    // - megablast: greedy alignment (eGreedyScoreOnly) - fast, good for high-identity sequences
    // - blastn: DP-based alignment (eDynProgScoreOnly) - slower, but handles divergent sequences better
    // This follows NCBI BLAST's approach for better alignment of divergent sequences
    let use_dp = match args.task.as_str() {
        "megablast" => false, // Use greedy for megablast (high-identity sequences)
        _ => true,            // Use DP for blastn, dc-megablast, blastn-short (divergent sequences)
    };

    // SCAN STRIDE OPTIMIZATION (NCBI BLAST algorithm)
    // NCBI BLAST uses: scan_step = word_length - lut_word_length + 1
    // This will be recalculated after two-stage lookup is built if use_two_stage is true
    let scan_step = if args.scan_step > 0 {
        args.scan_step // User-specified value
    } else {
        // Auto-calculate based on word_size (will be overridden for two-stage lookup)
        if effective_word_size >= 16 {
            4 // For megablast-like word sizes, use stride of 4
        } else if effective_word_size >= 11 {
            // For blastn (word_size=11), initial value (will be recalculated for two-stage lookup)
            2
        } else {
            1 // For very small word sizes, check every position
        }
    };

    if args.verbose {
        eprintln!("[INFO] Scan stride optimization: scan_step={} (word_size={})", scan_step, effective_word_size);
    }

    eprintln!("Reading query & subject...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries
        .iter()
        .map(|r| {
            r.id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string()
        })
        .collect();

    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();

    if queries.is_empty() || subjects.is_empty() {
        return Ok(());
    }

    let db_len_total: usize = subjects.iter().map(|r| r.seq().len()).sum();
    let db_num_seqs: usize = subjects.len();

    let scoring_spec = NuclScoringSpec {
        reward,
        penalty,
        gap_open,
        gap_extend,
    };
    let params = lookup_nucl_params(&scoring_spec);

    // Choose between direct address table (O(1) lookup) and HashMap based on word size
    // Direct address table is much faster but requires 4^word_size memory
    let use_direct_lookup = effective_word_size <= MAX_DIRECT_LOOKUP_WORD_SIZE;

    // TWO-STAGE LOOKUP TABLE (NCBI BLAST algorithm)
    // NCBI BLAST uses two-stage lookup for word_size >= 11:
    // - For word_size >= 16 (megablast): lut_word_length = 8
    // - For word_size = 11 (blastn): lut_word_length = 8 (when approx_table_entries < 12000)
    // - word_length = effective_word_size for extension triggering
    // This provides O(1) lookup performance even for large word sizes
    // Following NCBI BLAST's approach: use two-stage lookup when word_size >= 11
    let use_two_stage = effective_word_size >= 11;
    let lut_word_length = if use_two_stage {
        if effective_word_size >= 16 {
            8 // NCBI BLAST default for megablast (word_size >= 16)
        } else if effective_word_size == 11 {
            // For blastn (word_size=11), use lut_word_length=8 following NCBI BLAST
            // This matches NCBI BLAST's choice when approx_table_entries < 12000
            8
        } else {
            // For other word sizes (12-15), use 8 as default
            8
        }
    } else {
        effective_word_size.min(MAX_DIRECT_LOOKUP_WORD_SIZE)
    };

    // Apply DUST filter to mask low-complexity regions in query sequences
    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        eprintln!(
            "Applying DUST filter (level={}, window={}, linker={})...",
            args.dust_level, args.dust_window, args.dust_linker
        );
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        let masks: Vec<Vec<MaskedInterval>> = queries
            .iter()
            .map(|record| masker.mask_sequence(record.seq()))
            .collect();
        
        // Report DUST masking statistics
        let total_masked: usize = masks.iter().map(|m| m.iter().map(|i| i.end - i.start).sum::<usize>()).sum();
        let total_bases: usize = queries.iter().map(|r| r.seq().len()).sum();
        if total_bases > 0 {
            eprintln!(
                "DUST masked {} bases ({:.2}%) across {} sequences",
                total_masked,
                100.0 * total_masked as f64 / total_bases as f64,
                queries.len()
            );
        }
        masks
    } else {
        vec![Vec::new(); queries.len()]
    };

    eprintln!(
        "Building lookup (Task: {}, Word: {}, TwoStage: {}, LUTWord: {}, Direct: {}, DUST: {})...",
        args.task, effective_word_size, use_two_stage, lut_word_length, use_direct_lookup, args.dust
    );

    // Build the appropriate lookup table
    let two_stage_lookup: Option<TwoStageLookup> = if use_two_stage {
        Some(build_two_stage_lookup(&queries, effective_word_size, lut_word_length, &query_masks))
    } else {
        None
    };
    let pv_direct_lookup: Option<PvDirectLookup> = if !use_two_stage && use_direct_lookup {
        Some(build_pv_direct_lookup(&queries, effective_word_size, &query_masks))
    } else {
        None
    };
    let hash_lookup: Option<KmerLookup> = if !use_two_stage && !use_direct_lookup {
        Some(build_lookup(&queries, effective_word_size, &query_masks))
    } else {
        None
    };
    
    // Update scan_step for two-stage lookup
    let scan_step = if use_two_stage {
        // NCBI BLAST formula: scan_step = word_length - lut_word_length + 1
        (effective_word_size as isize - lut_word_length as isize + 1).max(1) as usize
    } else {
        scan_step
    };
    
    if args.verbose && use_two_stage {
        eprintln!("[INFO] Using two-stage lookup: lut_word_length={}, word_length={}, scan_step={}", 
                  lut_word_length, effective_word_size, scan_step);
    }

    if args.verbose {
        eprintln!("Searching...");
    }

    let bar = ProgressBar::new(subjects.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    // Channel for sending hits and sequence data
    // Use Option to signal completion: None means "all subjects processed"
    let (tx, rx) = channel::<Option<(Vec<Hit>, Vec<(String, String, Vec<u8>, Vec<u8>)>)>>();
    let out_path = args.out.clone();
    // Note: reward, penalty, gap_open, gap_extend are already defined above with task-adjusted values
    let params_clone = params.clone();

    // Keep a sender for the main thread to send the completion signal
    let tx_main = tx.clone();
    let verbose = args.verbose;
    let chain_enabled = args.chain;

    let writer_handle = std::thread::spawn(move || -> Result<()> {
        if verbose {
            eprintln!("[INFO] Writer thread started, waiting for hits...");
        }
        let mut all_hits = Vec::new();
        let mut all_sequences: FxHashMap<(String, String), (Vec<u8>, Vec<u8>)> =
            FxHashMap::default();
        let mut messages_received = 0usize;

        while let Ok(msg) = rx.recv() {
            match msg {
                Some((hits, seq_data)) => {
                    messages_received += 1;
                    if verbose && (messages_received == 1 || messages_received % 100 == 0) {
                        eprintln!(
                            "[INFO] Received message #{}, {} hits so far",
                            messages_received,
                            all_hits.len() + hits.len()
                        );
                    }
                    all_hits.extend(hits);
                    for (q_id, s_id, q_seq, s_seq) in seq_data {
                        all_sequences.entry((q_id, s_id)).or_insert((q_seq, s_seq));
                    }
                }
                None => {
                    // Completion signal received
                    if verbose {
                        eprintln!(
                            "[INFO] Completion signal received after {} messages",
                            messages_received
                        );
                    }
                    break;
                }
            }
        }

        if verbose {
            eprintln!(
                "[INFO] Received {} raw hits total, starting post-processing...",
                all_hits.len()
            );
        }
        let chain_start = std::time::Instant::now();

        // Chain nearby HSPs into longer alignments using cluster-then-extend
        // When chain_enabled is false (default), this just returns individual HSPs for BLAST compatibility
        let filtered_hits = chain_and_filter_hsps(
            all_hits,
            &all_sequences,
            reward,
            penalty,
            gap_open,
            gap_extend,
            db_len_total,
            db_num_seqs,
            &params_clone,
            use_dp,
            verbose,
            chain_enabled,
        );

        if verbose {
            eprintln!(
                "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
                chain_start.elapsed().as_secs_f64(),
                filtered_hits.len()
            );
        }
        let write_start = std::time::Instant::now();

        write_output(&filtered_hits, out_path.as_ref())?;

        if verbose {
            eprintln!(
                "[INFO] Output written in {:.2}s",
                write_start.elapsed().as_secs_f64()
            );
        }
        Ok(())
    });

    // 修正: s_idx -> _s_idx (未使用抑制)
    // Debug mode: set BLEMIR_DEBUG=1 to enable, BLEMIR_DEBUG_WINDOW="q_start-q_end,s_start-s_end" to focus on a region
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    let debug_window: Option<(usize, usize, usize, usize)> =
        std::env::var("BLEMIR_DEBUG_WINDOW").ok().and_then(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                let q_parts: Vec<usize> =
                    parts[0].split('-').filter_map(|x| x.parse().ok()).collect();
                let s_parts: Vec<usize> =
                    parts[1].split('-').filter_map(|x| x.parse().ok()).collect();
                if q_parts.len() == 2 && s_parts.len() == 2 {
                    return Some((q_parts[0], q_parts[1], s_parts[0], s_parts[1]));
                }
            }
            None
        });

    // Check if mask should be disabled for debugging
    let disable_mask = std::env::var("BLEMIR_DEBUG_NO_MASK").is_ok();

    if debug_mode {
        // Build marker to verify correct code is running
        eprintln!(
            "[DEBUG] BLEMIR build: 2024-12-24-v7 (adaptive banding with MAX_WINDOW_SIZE=50000)"
        );
        eprintln!(
            "[DEBUG] Task: {}, Scoring: reward={}, penalty={}, gap_open={}, gap_extend={}",
            args.task, reward, penalty, gap_open, gap_extend
        );
        if disable_mask {
            eprintln!("[DEBUG] Mask DISABLED for debugging");
        }
        if let Some((q_start, q_end, s_start, s_end)) = debug_window {
            eprintln!(
                "[DEBUG] Focusing on window: query {}-{}, subject {}-{}",
                q_start, q_end, s_start, s_end
            );
        }
    }

    // 修正: tx はここで所有権ごと渡す。for_each_with 終了時に自動でDropされるため、明示的な drop(tx) は不要。
    // Also pass two_stage_lookup reference and queries for use in the closure
    let two_stage_lookup_ref = two_stage_lookup.as_ref();
    let queries_ref = &queries;
    let query_ids_ref = &query_ids;
    subjects
        .par_iter()
        .enumerate()
        .for_each_with(tx, |tx, (_s_idx, s_record)| {
            let queries = queries_ref;
            let query_ids = query_ids_ref;
            let mut local_hits: Vec<Hit> = Vec::new();
            let mut local_sequences: Vec<(String, String, Vec<u8>, Vec<u8>)> = Vec::new();
            let mut seen_pairs: std::collections::HashSet<(String, String)> =
                std::collections::HashSet::new();

            let s_seq = s_record.seq();
            let s_len = s_seq.len();
            let s_id = s_record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();

            if s_len < effective_word_size {
                return;
            }

            // Debug counters for this subject
            let mut dbg_total_s_positions = 0usize;
            let dbg_ambiguous_skipped = 0usize;
            let dbg_no_lookup_match = 0usize;
            let mut dbg_seeds_found = 0usize;
            let mut dbg_mask_skipped = 0usize;
            let mut dbg_ungapped_low = 0usize;
            let mut dbg_two_hit_failed = 0usize;
            let mut dbg_gapped_attempted = 0usize;
            let mut dbg_window_seeds = 0usize;

            let safe_k = effective_word_size.min(31);

            // PERFORMANCE OPTIMIZATION: Rolling k-mer scanner with O(1) sliding window
            // Instead of packing the entire sequence first (which adds O(n) overhead),
            // we compute k-mers on-the-fly using a rolling window approach.
            // This achieves O(1) per-position k-mer extraction without allocation overhead.
            //
            // Encoding: A=0, C=1, G=2, T=3 (same as NCBI BLAST's ncbi2na)
            // The mask ensures we only keep the rightmost 2*k bits.
            let kmer_mask: u64 = (1u64 << (2 * safe_k)) - 1;

            // Lookup table for ASCII to 2-bit encoding (0xFF = invalid/ambiguous)
            const ENCODE_LUT: [u8; 256] = {
                let mut lut = [0xFFu8; 256];
                lut[b'A' as usize] = 0;
                lut[b'a' as usize] = 0;
                lut[b'C' as usize] = 1;
                lut[b'c' as usize] = 1;
                lut[b'G' as usize] = 2;
                lut[b'g' as usize] = 2;
                lut[b'T' as usize] = 3;
                lut[b't' as usize] = 3;
                lut[b'U' as usize] = 3;
                lut[b'u' as usize] = 3;
                lut
            };

            // Generate reverse complement for minus strand search
            let s_seq_rc = reverse_complement(s_seq);

            // Search both strands: plus (forward) and minus (reverse complement)
            // is_minus_strand: false = plus strand, true = minus strand
            for is_minus_strand in [false, true] {
                // Select the sequence to search based on strand
                let search_seq: &[u8] = if is_minus_strand {
                    &s_seq_rc
                } else {
                    s_seq
                };

                // PERFORMANCE OPTIMIZATION: For single query, use direct array indexing instead of HashMap
                // This eliminates millions of HashMap lookups and provides O(1) array access
                // Diagonal range: from -s_len (query before subject) to +s_len (query after subject)
                // We use an offset to map negative diagonals to positive array indices
                // However, for very large subject sequences, array initialization overhead can be significant,
                // so we only use array indexing for smaller sequences (< 1M bases)
                const MAX_ARRAY_DIAG_SIZE: usize = 2_000_000; // ~1M diagonals max (for ~500KB subject)
                let use_array_indexing = queries.len() == 1 && (s_len * 2 + 1) <= MAX_ARRAY_DIAG_SIZE;
                let diag_offset = if use_array_indexing {
                    s_len as isize // Offset to make all diagonals non-negative
                } else {
                    0
                };
                let diag_array_size = if use_array_indexing {
                    (s_len * 2) + 1 // Range: [-s_len, +s_len] -> [0, 2*s_len]
                } else {
                    0
                };

                // Mask to track already-extended regions on each diagonal
                // For single query: use Vec for O(1) access, otherwise use HashMap
                let mut mask_array: Vec<Option<usize>> = if use_array_indexing {
                    vec![None; diag_array_size]
                } else {
                    Vec::new()
                };
                let mut mask_hash: FxHashMap<u64, usize> = if !use_array_indexing {
                    FxHashMap::default()
                } else {
                    FxHashMap::default()
                };

                // Two-hit tracking: stores the last seed position on each diagonal for each query
                // For single query: use Vec for O(1) access, otherwise use HashMap
                let mut last_seed_array: Vec<Option<usize>> = if use_array_indexing {
                    vec![None; diag_array_size]
                } else {
                    Vec::new()
                };
                let mut last_seed_hash: FxHashMap<u64, usize> = if !use_array_indexing {
                    FxHashMap::default()
                } else {
                    FxHashMap::default()
                };

                // TWO-STAGE LOOKUP: Use separate rolling k-mer for lut_word_length
                if let Some(two_stage) = two_stage_lookup_ref {
                    // For two-stage lookup, use lut_word_length (8) for scanning
                    let lut_word_length = two_stage.lut_word_length();
                    let word_length = two_stage.word_length();
                    let lut_kmer_mask: u64 = (1u64 << (2 * lut_word_length)) - 1;
                    
                    // Rolling k-mer state for lut_word_length
                    let mut current_lut_kmer: u64 = 0;
                    let mut valid_bases: usize = 0;
                    
                    // Scan through the subject sequence with rolling lut_word_length k-mer
                    for s_pos in 0..s_len {
                        let base = search_seq[s_pos];
                        let code = ENCODE_LUT[base as usize];
                        
                        if code == 0xFF {
                            // Ambiguous base - reset the rolling window
                            valid_bases = 0;
                            current_lut_kmer = 0;
                            continue;
                        }
                        
                        // Shift in the new base
                        current_lut_kmer = ((current_lut_kmer << 2) | (code as u64)) & lut_kmer_mask;
                        valid_bases += 1;
                        
                        // Only process if we have a complete lut_word_length k-mer
                        if valid_bases < lut_word_length {
                            continue;
                        }
                        
                        // Calculate the starting position of this k-mer
                        let kmer_start = s_pos + 1 - lut_word_length;
                        
                        // SCAN STRIDE OPTIMIZATION: Skip positions based on scan_step
                        if kmer_start % scan_step != 0 {
                            continue;
                        }
                        
                        dbg_total_s_positions += 1;
                        
                        // Lookup using lut_word_length k-mer
                        let matches_slice = two_stage.get_hits(current_lut_kmer);
                        
                        // For each match, check if word_length match exists
                        for &(q_idx, q_pos) in matches_slice {
                            dbg_seeds_found += 1;
                            
                            // Check if word_length match exists starting at these positions
                            // Need to verify that subject[kmer_start..kmer_start+word_length] 
                            // matches query[q_pos..q_pos+word_length]
                            let q_seq = queries[q_idx as usize].seq();
                            let q_pos_usize = q_pos as usize;
                            
                            // Skip if sequences are too short
                            if q_pos_usize + word_length > q_seq.len() || kmer_start + word_length > s_len {
                                continue;
                            }
                            
                            // Use optimized SIMD comparison for 28-mer match check
                            // For word_length == 28, use specialized 28-mer SIMD function
                            let matches = if word_length == 28 {
                                compare_28mer_simd(
                                    &q_seq[q_pos_usize..q_pos_usize + word_length],
                                    &search_seq[kmer_start..kmer_start + word_length],
                                )
                            } else {
                                compare_sequences_simd(
                                    &q_seq[q_pos_usize..q_pos_usize + word_length],
                                    &search_seq[kmer_start..kmer_start + word_length],
                                    word_length,
                                )
                            };
                            
                            if !matches {
                                continue; // Skip if full word_length match doesn't exist
                            }
                            
                            let diag = kmer_start as isize - q_pos_usize as isize;

                            // Check if this seed is in the debug window
                            let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                                q_pos_usize >= q_start && q_pos_usize <= q_end &&
                                kmer_start >= s_start && kmer_start <= s_end
                            } else {
                                false
                            };

                            if in_window {
                                dbg_window_seeds += 1;
                            }

                            // Check if this region was already extended (skip if mask is disabled for debugging)
                            // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                            if !disable_mask {
                                let should_skip = if use_array_indexing {
                                    let diag_idx = (diag + diag_offset) as usize;
                                    if diag_idx < diag_array_size {
                                        if let Some(last_s_end) = mask_array[diag_idx] {
                                            kmer_start < last_s_end
                                        } else {
                                            false
                                        }
                                    } else {
                                        false
                                    }
                                } else {
                                    let diag_key = pack_diag_key(q_idx, diag);
                                    if let Some(&last_s_end) = mask_hash.get(&diag_key) {
                                        kmer_start < last_s_end
                                    } else {
                                        false
                                    }
                                };
                                if should_skip {
                                    dbg_mask_skipped += 1;
                                    if in_window && debug_mode {
                                        eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED by mask", q_pos_usize, kmer_start);
                                    }
                                    continue;
                                }
                            }

                            // PERFORMANCE OPTIMIZATION: Apply two-hit filter BEFORE ungapped extension
                            // This is the key optimization that makes NCBI BLAST fast
                            // Most seeds don't have a nearby seed on the same diagonal, so we skip them early
                            // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                            let trigger_extension = if use_array_indexing {
                                let diag_idx = (diag + diag_offset) as usize;
                                if diag_idx < diag_array_size {
                                    if let Some(prev_s_pos) = last_seed_array[diag_idx] {
                                        kmer_start.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW
                                    } else {
                                        false
                                    }
                                } else {
                                    false
                                }
                            } else {
                                let diag_key = pack_diag_key(q_idx, diag);
                                if let Some(&prev_s_pos) = last_seed_hash.get(&diag_key) {
                                    kmer_start.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW
                                } else {
                                    false
                                }
                            };

                            // Update the last seed position for this diagonal
                            // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                            if use_array_indexing {
                                let diag_idx = (diag + diag_offset) as usize;
                                if diag_idx < diag_array_size {
                                    last_seed_array[diag_idx] = Some(kmer_start);
                                }
                            } else {
                                let diag_key = pack_diag_key(q_idx, diag);
                                last_seed_hash.insert(diag_key, kmer_start);
                            }

                            // Skip ungapped extension if two-hit requirement is not met
                            // This is the key optimization - we skip most seeds without doing any extension
                            if !trigger_extension {
                                dbg_two_hit_failed += 1;
                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: two-hit not met", q_pos_usize, kmer_start);
                                }
                                continue;
                            }

                            let q_record = &queries[q_idx as usize];

                            // Now do ungapped extension (only for seeds that passed two-hit filter)
                            let (qs, qe, ss, _se, ungapped_score) = extend_hit_ungapped(
                                q_record.seq(),
                                search_seq,
                                q_pos_usize,
                                kmer_start,
                                reward,
                                penalty,
                            );

                            // Skip if ungapped score is too low
                            // Use task-specific threshold: higher for blastn to reduce extension count
                            if ungapped_score < min_ungapped_score {
                                dbg_ungapped_low += 1;
                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos_usize, kmer_start, ungapped_score, min_ungapped_score);
                                }
                                continue;
                            }

                            dbg_gapped_attempted += 1;

                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> GAPPED EXTENSION (ungapped_score={}, seed_len={})", q_pos_usize, kmer_start, ungapped_score, qe - qs);
                            }

                            // Gapped extension with NCBI-style high X-drop for longer alignments
                            // Using X_DROP_GAPPED_FINAL directly to allow alignments to push through
                            // low-similarity regions (NCBI BLAST uses 100 for nucleotide)
                            let (
                                final_qs,
                                final_qe,
                                final_ss,
                                final_se,
                                score,
                                matches,
                                mismatches,
                                gaps,
                                gap_letters,
                                dp_cells,
                            ) = extend_gapped_heuristic(
                                q_record.seq(),
                                search_seq,
                                qs,
                                ss,
                                qe - qs,
                                reward,
                                penalty,
                                gap_open,
                                gap_extend,
                                X_DROP_GAPPED_FINAL,
                                use_dp, // Use DP for blastn task, greedy for megablast
                            );

                            // Calculate alignment length
                            let aln_len = matches + mismatches + gap_letters;

                            // Debug: show gapped extension results for window seeds
                            if in_window && debug_mode {
                                let identity = if aln_len > 0 { 100.0 * matches as f64 / aln_len as f64 } else { 0.0 };
                                eprintln!(
                                    "[DEBUG WINDOW] Gapped result: q={}-{}, s={}-{}, score={}, len={}, identity={:.1}%, gaps={}, dp_cells={}",
                                    final_qs, final_qe, final_ss, final_se, score, aln_len, identity, gap_letters, dp_cells
                                );
                            }

                            // Suppress unused variable warning when not in debug mode
                            let _ = dp_cells;

                            // Update mask (unless disabled for debugging via BLEMIR_DEBUG_NO_MASK=1)
                            // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                            if !disable_mask {
                                if use_array_indexing {
                                    let diag_idx = (diag + diag_offset) as usize;
                                    if diag_idx < diag_array_size {
                                        mask_array[diag_idx] = Some(final_se);
                                    }
                                } else {
                                    let diag_key = pack_diag_key(q_idx, diag);
                                    mask_hash.insert(diag_key, final_se);
                                }
                            }

                            let (bit_score, eval) = calculate_evalue(
                                score,
                                q_record.seq().len(),
                                db_len_total,
                                db_num_seqs,
                                &params,
                            );

                            if eval <= args.evalue {
                                // Identity is matches / alignment_length, capped at 100%
                                let identity = if aln_len > 0 {
                                    ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                                } else {
                                    0.0
                                };

                                let q_id = query_ids[q_idx as usize].clone();

                                // Store sequence data for cluster-then-extend chaining
                                let pair_key = (q_id.clone(), s_id.clone());
                                if !seen_pairs.contains(&pair_key) {
                                    seen_pairs.insert(pair_key);
                                    local_sequences.push((
                                        q_id.clone(),
                                        s_id.clone(),
                                        q_record.seq().to_vec(),
                                        s_seq.to_vec(),
                                    ));
                                }

                                // Convert coordinates for minus strand hits
                                // For minus strand: positions in reverse complement need to be converted
                                // back to original subject coordinates, with s_start > s_end to indicate minus strand
                                let (hit_s_start, hit_s_end) = if is_minus_strand {
                                    // Convert from reverse complement coordinates to original coordinates
                                    // In reverse complement: position 0 corresponds to original position s_len-1
                                    // final_ss and final_se are 0-based positions in the reverse complement
                                    // We need to convert them to 1-based positions in the original sequence
                                    // with s_start > s_end to indicate minus strand
                                    let orig_s_start = s_len - final_ss; // 1-based, larger value
                                    let orig_s_end = s_len - final_se + 1; // 1-based, smaller value
                                    (orig_s_start, orig_s_end)
                                } else {
                                    (final_ss + 1, final_se)
                                };

                                local_hits.push(Hit {
                                    query_id: q_id.clone(),
                                    subject_id: s_id.clone(),
                                    identity,
                                    length: aln_len,
                                    mismatch: mismatches,
                                    gapopen: gaps,
                                    q_start: final_qs + 1,
                                    q_end: final_qe,
                                    s_start: hit_s_start,
                                    s_end: hit_s_end,
                                    e_value: eval,
                                    bit_score,
                                });
                            } // end of if eval <= args.evalue
                        } // end of for matches_slice in two-stage lookup
                    } // end of s_pos loop for two-stage lookup
                }
                if two_stage_lookup_ref.is_none() {
                    // Original lookup method (for non-two-stage lookup)
                    // Rolling k-mer state
                    let mut current_kmer: u64 = 0;
                    let mut valid_bases: usize = 0; // Count of consecutive valid bases in current window

                    // Scan through the subject sequence with rolling k-mer
                    for s_pos in 0..s_len {
                        let base = search_seq[s_pos];
                        let code = ENCODE_LUT[base as usize];

                    if code == 0xFF {
                        // Ambiguous base - reset the rolling window
                        valid_bases = 0;
                        current_kmer = 0;
                        continue;
                    }

                    // Shift in the new base
                    current_kmer = ((current_kmer << 2) | (code as u64)) & kmer_mask;
                    valid_bases += 1;

                    // Only process if we have a complete k-mer
                    if valid_bases < safe_k {
                        continue;
                    }

                    // Calculate the starting position of this k-mer
                    let kmer_start = s_pos + 1 - safe_k;

                    // SCAN STRIDE OPTIMIZATION: Skip positions based on scan_step
                    // We continue updating the rolling k-mer at every position (for correctness),
                    // but only perform the expensive lookup/extension work every scan_step positions.
                    // This reduces k-mer lookups by ~scan_step times with minimal sensitivity loss
                    // for large word sizes (megablast).
                    if kmer_start % scan_step != 0 {
                        continue;
                    }

                        dbg_total_s_positions += 1;

                        // Phase 2: Use PV-based direct lookup (O(1) with fast PV filtering) for word_size <= 13
                        // For word_size > 13, use hash-based lookup
                        let matches_slice: &[(u32, u32)] = if use_direct_lookup {
                            // Use PV for fast filtering before accessing the lookup table
                            pv_direct_lookup.as_ref().map(|pv_dl| pv_dl.get_hits_checked(current_kmer)).unwrap_or(&[])
                        } else {
                            // Use hash-based lookup for larger word sizes
                            hash_lookup.as_ref().and_then(|hl| hl.get(&current_kmer).map(|v| v.as_slice())).unwrap_or(&[])
                        };

                        for &(q_idx, q_pos) in matches_slice {
                            dbg_seeds_found += 1;
                        let diag = kmer_start as isize - q_pos as isize;

                        // Check if this seed is in the debug window
                        let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                            (q_pos as usize) >= q_start && (q_pos as usize) <= q_end &&
                            kmer_start >= s_start && kmer_start <= s_end
                        } else {
                            false
                        };

                        if in_window {
                            dbg_window_seeds += 1;
                        }

                        // Check if this region was already extended (skip if mask is disabled for debugging)
                        // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                        if !disable_mask {
                            let should_skip = if use_array_indexing {
                                let diag_idx = (diag + diag_offset) as usize;
                                if diag_idx < diag_array_size {
                                    if let Some(last_s_end) = mask_array[diag_idx] {
                                        kmer_start < last_s_end
                                    } else {
                                        false
                                    }
                                } else {
                                    false
                                }
                            } else {
                                let diag_key = pack_diag_key(q_idx, diag);
                                if let Some(&last_s_end) = mask_hash.get(&diag_key) {
                                    kmer_start < last_s_end
                                } else {
                                    false
                                }
                            };
                            if should_skip {
                                dbg_mask_skipped += 1;
                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED by mask", q_pos, kmer_start);
                                }
                                continue;
                            }
                        }

                        // PERFORMANCE OPTIMIZATION: Apply two-hit filter BEFORE ungapped extension
                        // This is the key optimization that makes NCBI BLAST fast
                        // Most seeds don't have a nearby seed on the same diagonal, so we skip them early
                        // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                        let trigger_extension = if use_array_indexing {
                            let diag_idx = (diag + diag_offset) as usize;
                            if diag_idx < diag_array_size {
                                if let Some(prev_s_pos) = last_seed_array[diag_idx] {
                                    kmer_start.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW
                                } else {
                                    false
                                }
                            } else {
                                false
                            }
                        } else {
                            let diag_key = pack_diag_key(q_idx, diag);
                            if let Some(&prev_s_pos) = last_seed_hash.get(&diag_key) {
                                kmer_start.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW
                            } else {
                                false
                            }
                        };

                        // Update the last seed position for this diagonal
                        // PERFORMANCE OPTIMIZATION: Use direct array indexing for single query
                        if use_array_indexing {
                            let diag_idx = (diag + diag_offset) as usize;
                            if diag_idx < diag_array_size {
                                last_seed_array[diag_idx] = Some(kmer_start);
                            }
                        } else {
                            let diag_key = pack_diag_key(q_idx, diag);
                            last_seed_hash.insert(diag_key, kmer_start);
                        }

                        // Skip ungapped extension if two-hit requirement is not met
                        // This is the key optimization - we skip most seeds without doing any extension
                        if !trigger_extension {
                            dbg_two_hit_failed += 1;
                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: two-hit not met", q_pos, kmer_start);
                            }
                            continue;
                        }

                        let q_record = &queries[q_idx as usize];

                        // Now do ungapped extension (only for seeds that passed two-hit filter)
                        let (qs, qe, ss, _se, ungapped_score) = extend_hit_ungapped(
                            q_record.seq(),
                            search_seq,
                            q_pos as usize,
                            kmer_start,
                            reward,
                            penalty,
                        );

                        // Skip if ungapped score is too low
                        // Use task-specific threshold: higher for blastn to reduce extension count
                        if ungapped_score < min_ungapped_score {
                            dbg_ungapped_low += 1;
                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos, kmer_start, ungapped_score, min_ungapped_score);
                            }
                            continue;
                        }

                        dbg_gapped_attempted += 1;

                        if in_window && debug_mode {
                            eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> GAPPED EXTENSION (ungapped_score={}, seed_len={})", q_pos, kmer_start, ungapped_score, qe - qs);
                        }

                        // Gapped extension with NCBI-style high X-drop for longer alignments
                        // Using X_DROP_GAPPED_FINAL directly to allow alignments to push through
                        // low-similarity regions (NCBI BLAST uses 100 for nucleotide)
                        let (
                            final_qs,
                            final_qe,
                            final_ss,
                            final_se,
                            score,
                            matches,
                            mismatches,
                            gaps,
                            gap_letters,
                            dp_cells,
                        ) = extend_gapped_heuristic(
                            q_record.seq(),
                            search_seq,
                            qs,
                            ss,
                            qe - qs,
                            reward,
                            penalty,
                            gap_open,
                            gap_extend,
                            X_DROP_GAPPED_FINAL,
                            use_dp, // Use DP for blastn task, greedy for megablast
                        );

                        // Debug: show gapped extension results for window seeds
                        if in_window && debug_mode {
                            let aln_len = matches + mismatches + gap_letters;
                            eprintln!("[DEBUG WINDOW] Gapped extension result: q={}-{}, s={}-{}, score={}, len={}", final_qs, final_qe, final_ss, final_se, score, aln_len);
                        }

                        // Calculate statistics
                        let aln_len = matches + mismatches + gap_letters;
                        let identity = if aln_len > 0 {
                            (matches as f64) / (aln_len as f64)
                        } else {
                            0.0
                        };

                        // Calculate e-value and bit score using Karlin-Altschul statistics
                        let (bit_score, eval) = calculate_evalue(
                            score,
                            q_record.seq().len(),
                            db_len_total,
                            db_num_seqs,
                            &params,
                        );
                        
                        // Skip if e-value is too high
                        if eval > args.evalue {
                            continue;
                        }
                        
                        // Calculate alignment length and identity
                        let aln_len = matches + mismatches + gap_letters;
                        let identity = if aln_len > 0 {
                            ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                        } else {
                            0.0
                        };

                        // Store sequence data for post-processing (only for pairs we haven't seen)
                        let q_id = &query_ids[q_idx as usize];
                        if !seen_pairs.contains(&(q_id.clone(), s_id.clone())) {
                            seen_pairs.insert((q_id.clone(), s_id.clone()));
                            local_sequences.push((
                                q_id.clone(),
                                s_id.clone(),
                                q_record.seq().to_vec(),
                                s_seq.to_vec(),
                            ));
                        }

                        // Convert coordinates for minus strand hits
                        // For minus strand: positions in reverse complement need to be converted
                        // back to original subject coordinates, with s_start > s_end to indicate minus strand
                        let (hit_s_start, hit_s_end) = if is_minus_strand {
                            // Convert from reverse complement coordinates to original coordinates
                            // In reverse complement: position 0 corresponds to original position s_len-1
                            // final_ss and final_se are 0-based positions in the reverse complement
                            // We need to convert them to 1-based positions in the original sequence
                            // with s_start > s_end to indicate minus strand
                            let orig_s_start = s_len - final_ss; // 1-based, larger value
                            let orig_s_end = s_len - final_se + 1; // 1-based, smaller value
                            (orig_s_start, orig_s_end)
                        } else {
                            (final_ss + 1, final_se)
                        };

                        local_hits.push(Hit {
                            query_id: q_id.clone(),
                            subject_id: s_id.clone(),
                            identity,
                            length: aln_len,
                            mismatch: mismatches,
                            gapopen: gaps,
                            q_start: final_qs + 1,
                            q_end: final_qe,
                            s_start: hit_s_start,
                            s_end: hit_s_end,
                            e_value: eval,
                            bit_score,
                        });
                    }
                    } // end of s_pos loop for original lookup
                } // end of if/else for two-stage vs original lookup
            } // end of strand loop

            // Print debug summary for this subject
            if debug_mode {
                eprintln!(
                    "[DEBUG] Subject {}: positions={}, seeds_found={}, mask_skipped={}, ungapped_low={}, two_hit_failed={}, gapped_attempted={}, window_seeds={}, hits={}",
                    s_id, dbg_total_s_positions, dbg_seeds_found, dbg_mask_skipped, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds, local_hits.len()
                );
            }

            // Suppress unused variable warnings when not in debug mode
            let _ = (dbg_total_s_positions, dbg_ambiguous_skipped, dbg_no_lookup_match, dbg_seeds_found, dbg_mask_skipped, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds);

            if !local_hits.is_empty() {
                tx.send(Some((local_hits, local_sequences))).unwrap();
            }
            bar.inc(1);
        });

    bar.finish();
    if args.verbose {
        eprintln!("[INFO] Parallel processing complete, sending completion signal...");
    }

    // Send completion signal to writer thread
    // This ensures the writer thread exits even if sender-drop semantics are delayed
    tx_main.send(None).unwrap();
    drop(tx_main); // Explicitly drop to ensure channel closes

    writer_handle.join().unwrap()?;
    Ok(())
}
