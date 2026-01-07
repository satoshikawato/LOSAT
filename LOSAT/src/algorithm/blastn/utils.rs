use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::mpsc::channel;
use crate::common::{write_output, Hit};
use crate::stats::{lookup_nucl_params, KarlinParams};
use super::args::BlastnArgs;
use super::alignment::extend_gapped_heuristic;
use super::extension::{extend_hit_ungapped, type_of_word};
use super::constants::TWO_HIT_WINDOW;
use super::coordination::{configure_task, read_sequences, prepare_sequence_data, build_lookup_tables};
use super::lookup::{reverse_complement, pack_diag_key};
// SIMD comparison functions removed - now using type_of_word instead
use super::ncbi_cutoffs::{compute_blastn_cutoff_score, GAP_TRIGGER_BIT_SCORE_NUCL};
use super::query_info::QueryInfo;

// NCBI reference: blast_extend.h:57-60
// typedef struct DiagStruct {
//     Int4 last_hit;  /**< last hit position on this diagonal */
//     Uint4 flag;     /**< flag indicating if hit was saved/extended */
// } DiagStruct;
/// Diagonal structure for tracking hits (equivalent to NCBI's DiagStruct)
#[derive(Clone, Copy, Default)]
struct DiagStruct {
    /// Last hit position on this diagonal (with diag_offset added)
    /// NCBI reference: blast_extend.c:103: diag_struct_array[i].last_hit = -diag->window;
    /// LOSAT: Initialize to 0 (equivalent behavior for first hit)
    last_hit: usize,
    /// Flag indicating if hit was saved/extended (1 = extended, 0 = new hit)
    /// NCBI reference: na_ungapped.c:666: hit_saved = hit_level_array[real_diag].flag;
    flag: u8,
}

// Re-export for backward compatibility
pub use crate::algorithm::common::evalue::calculate_evalue_database_search as calculate_evalue;
/// Filter HSPs to remove redundant overlapping hits.
/// This function applies overlap filtering to match NCBI BLAST's behavior of outputting individual HSPs.
/// Chaining/clustering functionality has been removed to match NCBI BLAST blastn behavior.
pub fn filter_hsps(
    hits: Vec<Hit>,
    _sequences: &FxHashMap<(String, String), (Vec<u8>, Vec<u8>)>,
    _reward: i32,
    _penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    _db_len_total: usize,
    _db_num_seqs: usize,
    _params: &KarlinParams,
    _use_dp: bool,
    verbose: bool,
) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    let total_start = std::time::Instant::now();

    // Always output individual HSPs (NCBI BLAST behavior)
    // NCBI BLAST's blastn does not use chaining/clustering - all HSPs are saved individually
    let mut result_hits: Vec<Hit> = hits;

    // Final pass: remove redundant overlapping hits using NCBI BLAST's s_DominateTest algorithm
    // Reference: hspfilter_culling.c:s_DominateTest (lines 79-120)
    //
    // Key features:
    // - Uses query coordinates only (not subject coordinates)
    // - No diagonal gating (all HSPs are compared)
    // - Uses raw_score (not bit_score)
    // - Score/length tradeoff formula: d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
    // - 50% overlap check on candidate's length
    let filter_start = std::time::Instant::now();
    let result_hits_count = result_hits.len();

    // Sort by raw_score (descending) to match NCBI BLAST's sorting order
    // NCBI reference: blast_hits.c:ScoreCompareHSPs
    // Order: score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC
    result_hits.sort_by(|a, b| {
        b.raw_score
            .cmp(&a.raw_score)
            .then_with(|| a.s_start.cmp(&b.s_start))  // s_start ASC
            .then_with(|| b.s_end.cmp(&a.s_end))       // s_end DESC
            .then_with(|| a.q_start.cmp(&b.q_start))    // q_start ASC
            .then_with(|| b.q_end.cmp(&a.q_end))        // q_end DESC
    });

    // NCBI reference: hspfilter_culling.c:s_BlastHSPCullingRun
    // NCBI uses interval tree for O(n log n) culling, but for simplicity we use O(n²) comparison
    // which matches NCBI's s_DominateTest logic exactly
    // NO SPATIAL BINNING - NCBI does NOT use binning, it compares all HSPs
    let mut final_hits: Vec<Hit> = Vec::new();

    for hit in result_hits {
        // NCBI reference: hspfilter_culling.c:603-650
        // Check if this hit is dominated by any already-kept hit
        // NCBI: for each hsp in hsp_list, check if it's dominated by any kept hsp
        let mut dominated = false;
        for kept in &final_hits {
            // NCBI BLAST's s_DominateTest algorithm
            // Reference: hspfilter_culling.c:79-120
            
            // Must be same query-subject pair
            if hit.query_id != kept.query_id || hit.subject_id != kept.subject_id {
                continue;
            }

            // Query coordinates (already normalized to plus strand in blastn)
            // NCBI reference: hspfilter_culling.c:79-87
            // static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
            //     Int8 b1 = p->begin;  // p = already kept HSP
            //     Int8 b2 = y->begin;  // y = candidate HSP
            //     Int8 e1 = p->end;
            //     Int8 e2 = y->end;
            //     Int8 s1 = p->hsp->score;  // p's score
            //     Int8 s2 = y->hsp->score;  // y's score
            //     Int8 l1 = e1 - b1;  // p's length
            //     Int8 l2 = e2 - b2;  // y's length (candidate)
            // }
            // In LOSAT: kept = p (already kept), hit = y (candidate)
            // CRITICAL: Match NCBI variable mapping exactly
            let b1 = kept.q_start as i64;  // p->begin (already kept)
            let e1 = kept.q_end as i64;    // p->end (already kept)
            let b2 = hit.q_start as i64;   // y->begin (candidate)
            let e2 = hit.q_end as i64;     // y->end (candidate)
            
            // Query lengths
            // NCBI reference: hspfilter_culling.c:86-87
            // Int8 l1 = e1 - b1;  // p's length (already kept)
            // Int8 l2 = e2 - b2;  // y's length (candidate)
            let l1 = e1 - b1;  // kept's length (already kept)
            let l2 = e2 - b2;  // hit's length (candidate)
            
            // Calculate overlap (query coordinates only)
            // NCBI reference: hspfilter_culling.c:88: overlap = MIN(e1,e2) - MAX(b1,b2)
            // If overlap is negative (no overlap), it will fail the 50% check below
            let overlap = e1.min(e2) - b1.max(b2);
            
            // If not overlap by more than 50% of candidate's length, don't dominate
            // NCBI reference: hspfilter_culling.c:92: if(2 *overlap < l2) return FALSE;
            // l2 = y->end - y->begin (candidate's length)
            // Note: if overlap < 0 (no overlap), this condition is always true
            if overlap < 0 || 2 * overlap < l2 {
                continue;
            }
            
            // Raw scores (not bit_score)
            // NCBI reference: hspfilter_culling.c:84-85
            // Int8 s1 = p->hsp->score;  // already kept HSP's score
            // Int8 s2 = y->hsp->score;  // candidate HSP's score
            let s1 = kept.raw_score as i64;  // p's score (already kept)
            let s2 = hit.raw_score as i64;   // y's score (candidate)
            
            // Main criterion: 2 * (%diff in score) + 1 * (%diff in length)
            // Formula: d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
            // where s1, l1 = kept (already kept), s2, l2 = hit (candidate)
            // NCBI reference: hspfilter_culling.c:99
            let d = 4 * s1 * l1 + 2 * s1 * l2 - 2 * s2 * l1 - 4 * s2 * l2;
            
            // Tie-breaker for identical HSPs
            // NCBI reference: hspfilter_culling.c:101-113
            // if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
            //     if(s1 != s2) {
            //         return (s1>s2);  // p dominates y if p's score > y's score
            //     }
            //     if(p->sid != y->sid) {
            //         return (p->sid < y->sid);  // p dominates y if p->sid < y->sid
            //     }
            //     if(p->hsp->subject.offset > y->hsp->subject.offset) {
            //         return FALSE;  // p does NOT dominate y
            //     }
            //     return TRUE;  // p dominates y (default)
            // }
            if (s1 == s2 && b1 == b2 && l1 == l2) || d == 0 {
                if s1 != s2 {
                    // NCBI: return (s1>s2);  // p dominates y if p's score > y's score
                    // s1 = kept.raw_score (p's score), s2 = hit.raw_score (y's score)
                    if s1 > s2 {
                        // kept's score > hit's score → kept dominates hit
                        dominated = true;
                        break;
                    }
                    // kept's score <= hit's score → kept does NOT dominate hit
                    continue;
                }
                // Same score: use subject_id, then subject offset
                // NCBI: if(p->sid != y->sid) return (p->sid < y->sid);
                // p->sid = kept.subject_id, y->sid = hit.subject_id
                if kept.subject_id != hit.subject_id {
                    // NCBI: return (p->sid < y->sid);  // p dominates y if p->sid < y->sid
                    if kept.subject_id < hit.subject_id {
                        // kept.subject_id < hit.subject_id → kept dominates hit
                        dominated = true;
                        break;
                    }
                    // kept.subject_id >= hit.subject_id → kept does NOT dominate hit
                    continue;
                }
                // Same subject: use subject offset
                // NCBI: if(p->hsp->subject.offset > y->hsp->subject.offset) return FALSE;
                //       return TRUE;  (default: p dominates y)
                // p->hsp->subject.offset = kept.s_start, y->hsp->subject.offset = hit.s_start
                if kept.s_start > hit.s_start {
                    // NCBI: return FALSE;  // p does NOT dominate y
                    continue;
                }
                // NCBI: return TRUE;  // p dominates y (default)
                dominated = true;
                break;
            }
            
            // If d < 0: kept does NOT dominate hit (hit should be kept)
            // If d > 0: kept dominates hit (hit should be filtered)
            // NCBI reference: hspfilter_culling.c:115-119
            // if (d < 0) {
            //     return FALSE;  // p does NOT dominate y
            // }
            // return TRUE;  // p dominates y
            if d < 0 {
                // kept does NOT dominate hit
                continue;
            }
            
            // d >= 0: kept dominates hit (hit should be filtered)
            dominated = true;
            break;
        }

        if !dominated {
            // NCBI reference: hspfilter_culling.c:650
            // Add hit to final list if not dominated
            final_hits.push(hit);
        }
    }
    let filter_time = filter_start.elapsed();

        if verbose {
            eprintln!("[INFO] filter_hsps summary:");
            eprintln!(
                "[INFO]   Total time: {:.3}s",
                total_start.elapsed().as_secs_f64()
            );
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

    // Configure task-specific parameters
    let config = configure_task(&args);
    
    if args.verbose {
        eprintln!("[INFO] Scan stride optimization: scan_step={} (word_size={})", config.scan_step, config.effective_word_size);
    }

    // Read sequences
    let (queries, query_ids, subjects) = read_sequences(&args)?;
    if queries.is_empty() || subjects.is_empty() {
        return Ok(());
    }

    // Prepare sequence data (including DUST masking)
    let seq_data = prepare_sequence_data(&args, queries, query_ids, subjects);
    
    // Get Karlin-Altschul parameters
    use crate::config::NuclScoringSpec;
    let scoring_spec = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    let params = lookup_nucl_params(&scoring_spec);

    // Build lookup tables
    let (lookup_tables, scan_step) = build_lookup_tables(
        &config,
        &args,
        &seq_data.queries,
        &seq_data.query_masks,
    );

    // Initialize query info for context management
    // NCBI reference: blast_setup.c (query info initialization)
    // For blastn, each query has 2 contexts: forward strand (0) and reverse strand (1)
    let query_seqs: Vec<Vec<u8>> = seq_data.queries.iter().map(|r| r.seq().to_vec()).collect();
    let query_info = QueryInfo::new_blastn(&query_seqs);

    // Pre-compute cutoff scores per context
    // NCBI reference: blast_parameters.c:368-374
    // Per-context cutoff calculation in BlastInitialWordParametersUpdate
    // cutoffs[context].cutoff_score = MIN(gap_trigger * scale_factor, 
    //                                     cutoffs[context].cutoff_score_max);
    let evalue_threshold = args.evalue;
    let mut cutoff_scores: Vec<i32> = Vec::with_capacity(query_info.contexts.len());
    for context in &query_info.contexts {
        // For each context, compute cutoff_score using context-specific query_length
        // NCBI: cutoff_score calculation uses query_info->contexts[context].query_length
        let query_len = context.query_length as i64 * 2; // Both strands for blastn (already per-strand in context)
        // For blastn, we use the same query_length for both forward and reverse strands
        // The context already represents a single strand, so we don't multiply by 2
        let query_len_for_cutoff = context.query_length as i64;
        
        // Use a representative subject length for cutoff calculation
        // In practice, cutoff_score is computed per-subject, but we pre-compute
        // a baseline here. The actual cutoff may be updated per-subject.
        let subject_len = if !seq_data.subjects.is_empty() {
            seq_data.subjects[0].seq().len() as i64
        } else {
            1000 // Default fallback
        };
        
        let cutoff_score = compute_blastn_cutoff_score(
            query_len_for_cutoff,
            subject_len,
            evalue_threshold,
            GAP_TRIGGER_BIT_SCORE_NUCL,
            &params, // ungapped params (same as gapped for blastn)
            &params, // gapped params (same as ungapped for blastn)
            1.0, // scale_factor (typically 1.0 for standard scoring)
        );
        cutoff_scores.push(cutoff_score);
    }

    if args.verbose {
        eprintln!("Searching... ({} contexts, cutoff_scores: {:?})", query_info.contexts.len(), cutoff_scores);
    }

    let bar = ProgressBar::new(seq_data.subjects.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    // Channel for sending hits and sequence data
    // Use Option to signal completion: None means "all subjects processed"
    let (tx, rx) = channel::<Option<(Vec<Hit>, Vec<(String, String, Vec<u8>, Vec<u8>)>)>>();
    let out_path = args.out.clone();
    let params_clone = params.clone();

    // Keep a sender for the main thread to send the completion signal
    let tx_main = tx.clone();
    let verbose = args.verbose;

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

        // Filter HSPs to remove redundant overlapping hits
        // Always outputs individual HSPs (NCBI BLAST behavior - no chaining/clustering)
        let filtered_hits = filter_hsps(
            all_hits,
            &all_sequences,
            config.reward,
            config.penalty,
            config.gap_open,
            config.gap_extend,
            seq_data.db_len_total,
            seq_data.db_num_seqs,
            &params_clone,
            config.use_dp,
            verbose,
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

    // NCBI BLAST does NOT use mask_array/mask_hash for diagonal suppression
    // REMOVED: disable_mask variable (no longer needed)

    if debug_mode {
        // Build marker to verify correct code is running
        eprintln!(
            "[DEBUG] BLEMIR build: 2024-12-24-v7 (adaptive banding with MAX_WINDOW_SIZE=50000)"
        );
        eprintln!(
            "[DEBUG] Task: {}, Scoring: reward={}, penalty={}, gap_open={}, gap_extend={}",
            args.task, config.reward, config.penalty, config.gap_open, config.gap_extend
        );
        if let Some((q_start, q_end, s_start, s_end)) = debug_window {
            eprintln!(
                "[DEBUG] Focusing on window: query {}-{}, subject {}-{}",
                q_start, q_end, s_start, s_end
            );
        }
    }

    // Pass lookup tables and queries for use in the closure
    let two_stage_lookup_ref = lookup_tables.two_stage_lookup.as_ref();
    let pv_direct_lookup_ref = lookup_tables.pv_direct_lookup.as_ref();
    let hash_lookup_ref = lookup_tables.hash_lookup.as_ref();
    let queries_ref = &seq_data.queries;
    let query_ids_ref = &seq_data.query_ids;
    let query_masks_ref = &seq_data.query_masks;
    
    // Capture config values for use in closure
    let effective_word_size = config.effective_word_size;
    let min_ungapped_score = config.min_ungapped_score;
    let use_dp = config.use_dp;
    let use_direct_lookup = config.use_direct_lookup;
    let reward = config.reward;
    let penalty = config.penalty;
    let gap_open = config.gap_open;
    let gap_extend = config.gap_extend;
    let x_drop_gapped = config.x_drop_gapped;
    let scan_range = config.scan_range; // For off-diagonal hit detection
    let db_len_total = seq_data.db_len_total;
    let db_num_seqs = seq_data.db_num_seqs;
    let params_for_closure = params.clone();
    let evalue_threshold = args.evalue;
    let query_info_ref = &query_info;
    let cutoff_scores_ref = &cutoff_scores;
    
    seq_data.subjects
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

            // NCBI reference: blast_parameters.c:368-374
            // Per-subject cutoff_score update
            // NCBI: BlastInitialWordParametersUpdate recalculates cutoff_score per-subject
            // because cutoff_score depends on subject length (used in eff_searchsp calculation)
            // Reference: blast_parameters.c:348-374
            let subject_len = s_len as i64;
            let mut per_subject_cutoff_scores = cutoff_scores_ref.clone();
            
            // NCBI reference: blast_parameters.c:368-374
            // Recompute cutoff_score for each context with actual subject length
            // NCBI: cutoff_score calculation uses subject length in eff_searchsp
            for (ctx_idx, context) in query_info_ref.contexts.iter().enumerate() {
                if ctx_idx < per_subject_cutoff_scores.len() {
                    // NCBI reference: blast_parameters.c:368-374
                    // Recompute cutoff_score for this context with actual subject length
                    let query_len = context.query_length as i64;
                    let updated_cutoff = compute_blastn_cutoff_score(
                        query_len,
                        subject_len,
                        evalue_threshold,
                        GAP_TRIGGER_BIT_SCORE_NUCL,
                        &params_for_closure,
                        &params_for_closure,
                        1.0,
                    );
                    per_subject_cutoff_scores[ctx_idx] = updated_cutoff;
                }
            }

            // Debug counters for this subject
            let mut dbg_total_s_positions = 0usize;
            let dbg_ambiguous_skipped = 0usize;
            let dbg_no_lookup_match = 0usize;
            let mut dbg_seeds_found = 0usize;
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

            // NCBI reference: blast_extend.c:52-61, 141
            // s_BlastDiagTableNew(query_length, multiple_hits, word_params->options->window_size)
            // NCBI creates diag_table per query with query_length
            // diag_array_length is calculated as smallest power of 2 >= (qlen + window_size)
            // NCBI processes each query separately with its own diag_table
            // Reference: blast_extend.c:141
            // ewp->diag_table = diag_table = s_BlastDiagTableNew(query_length, multiple_hits, word_params->options->window_size);
            
            // Pre-compute diag_array_length for each query (NCBI creates diag_table per query)
            let query_diag_lengths: Vec<(usize, usize)> = queries.iter().map(|q| {
                let q_len = q.seq().len();
                let mut diag_array_length = 1usize;
                // NCBI reference: blast_extend.c:54-56
                // while (diag_array_length < (qlen+window_size)) { diag_array_length = diag_array_length << 1; }
                while diag_array_length < (q_len + TWO_HIT_WINDOW) {
                    diag_array_length <<= 1;
                }
                let diag_mask = diag_array_length - 1;
                (diag_array_length, diag_mask)
            }).collect();

            // Search both strands: plus (forward) and minus (reverse complement)
            // is_minus_strand: false = plus strand, true = minus strand
            for is_minus_strand in [false, true] {
                // Select the sequence to search based on strand
                let search_seq: &[u8] = if is_minus_strand {
                    &s_seq_rc
                } else {
                    s_seq
                };
                
                // NCBI reference: na_ungapped.c:663-664
                // NCBI uses real_diag (masked diag) directly as array index
                // Array size is diag_array_length (power of 2)
                // NCBI reference: blast_extend.c:141
                // NCBI creates diag_table per query, so each query has its own hit_level_array
                // For multiple queries, NCBI processes each query separately
                // LOSAT uses HashMap for multiple queries (NCBI would create separate diag_table for each)
                const MAX_ARRAY_DIAG_SIZE: usize = 2_000_000; // ~1M diagonals max
                // NCBI creates diag_table per query, so check if all queries can use array indexing
                let all_queries_small = queries.len() == 1 && query_diag_lengths[0].0 <= MAX_ARRAY_DIAG_SIZE;
                let use_array_indexing = all_queries_small;
                // NCBI does NOT use diag_offset - it uses real_diag (masked) directly
                let diag_offset = 0isize; // Not used in NCBI approach

                // NCBI BLAST does NOT use a separate mask array for diagonal suppression
                // NCBI only uses last_hit in hit_level_array (checked at line 753)
                // This mask_array/mask_hash was LOSAT-specific and caused excessive filtering
                // REMOVED to match NCBI BLAST behavior

                // NCBI reference: blast_extend.h:77-80, na_ungapped.c:660-666
                // Two-hit tracking: DiagStruct array for tracking last_hit and flag per diagonal
                // hit_level_array[diag] contains {last_hit, flag} (equivalent to NCBI's DiagStruct)
                // hit_len_array[diag] contains hit length (0 = no hit, >0 = hit length)
                // NCBI reference: blast_extend.c:141, 145-149
                // NCBI creates diag_table per query with query_length
                // diag_table->hit_level_array = calloc(diag_table->diag_array_length, sizeof(DiagStruct))
                // diag_table->hit_len_array = calloc(diag_table->diag_array_length, sizeof(Uint1))
                // NCBI processes each query separately, so each query has its own hit_level_array
                // LOSAT processes multiple queries together, so we use HashMap keyed by (q_idx, real_diag)
                // For single query: use Vec for O(1) access (matches NCBI), otherwise use HashMap
                let diag_array_size = if use_array_indexing {
                    query_diag_lengths[0].0 // NCBI uses diag_array_length directly (per query)
                } else {
                    0
                };
                let mut hit_level_array: Vec<DiagStruct> = if use_array_indexing {
                    vec![DiagStruct::default(); diag_array_size]
                } else {
                    Vec::new()
                };
                // NCBI creates diag_table per query, so for multiple queries we need per-query tracking
                // Use HashMap with (q_idx, real_diag) as key to match NCBI behavior
                let mut hit_level_hash: FxHashMap<u64, DiagStruct> = if !use_array_indexing {
                    FxHashMap::default()
                } else {
                    FxHashMap::default()
                };
                
                // NCBI reference: na_ungapped.c:696, 705: diag_table->hit_len_array[off_diag]
                // Hit length array: stores the length of the last hit on each diagonal
                // 0 = no hit, >0 = hit length
                // NCBI reference: blast_extend.c:148-149
                // diag_table->hit_len_array = calloc(diag_table->diag_array_length, sizeof(Uint1))
                // NCBI creates diag_table per query, so each query has its own hit_len_array
                let mut hit_len_array: Vec<usize> = if use_array_indexing {
                    vec![0; diag_array_size]
                } else {
                    Vec::new()
                };
                // NCBI creates diag_table per query, so for multiple queries we need per-query tracking
                // Use HashMap with (q_idx, real_diag) as key to match NCBI behavior
                let mut hit_len_hash: FxHashMap<u64, usize> = if !use_array_indexing {
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
                            let q_record = &queries[q_idx as usize];
                            let q_seq = q_record.seq();
                            let q_pos_usize = q_pos as usize;
                            
                            // NCBI reference: na_ungapped.c:730-731
                            // Int4 context = BSearchContextInfo(q_off, query_info);
                            // cutoffs = word_params->cutoffs + context;
                            // For blastn, we need to determine context from q_off
                            // However, q_off is relative to the query, not the concatenated super-query
                            // For blastn, each query has 2 contexts: 0 (forward) and 1 (reverse)
                            // We need to determine which strand we're on based on the query sequence
                            // For now, we'll use q_idx * 2 as the base context (forward strand)
                            // and determine if we need reverse strand based on the actual query
                            // In practice, blastn processes forward and reverse strands separately
                            // For simplicity, we'll use context 0 (forward) for now
                            // TODO: Determine strand from actual query processing
                            let context = (q_idx as usize) * 2; // Forward strand context
                            let cutoff_score = if context < cutoff_scores_ref.len() {
                                per_subject_cutoff_scores[context]
                            } else {
                                // Fallback: compute cutoff_score if context is out of bounds
                                let query_len = q_seq.len() as i64;
                                let subject_len = s_len as i64;
                                compute_blastn_cutoff_score(
                                    query_len,
                                    subject_len,
                                    evalue_threshold,
                                    GAP_TRIGGER_BIT_SCORE_NUCL,
                                    &params_for_closure,
                                    &params_for_closure,
                                    1.0,
                                )
                            };
                            
                            // NCBI reference: na_ungapped.c:1081-1140
                            // CRITICAL: In two-stage lookup, seed finding phase ONLY finds lut_word_length matches.
                            // The word_length matching happens LATER in s_BlastnExtendInitialHit (extension phase).
                            //
                            // NCBI BLAST flow:
                            // 1. Seed finding: Find all lut_word_length matches (this phase)
                            // 2. Extension: For each seed, call s_BlastnExtendInitialHit which:
                            //    - Does LEFT extension (backwards from lut_word_length start)
                            //    - Does RIGHT extension (forwards from lut_word_length end)
                            //    - Verifies word_length match (ext_left + ext_right >= ext_to)
                            //
                            // We should NOT filter seeds during seed finding - pass ALL seeds to extension.
                            // The extension phase will handle word_length matching.
                            
                            // Skip only if sequences are too short for lut_word_length match
                            // (word_length check happens in extension phase)
                            if q_pos_usize + lut_word_length > q_seq.len() || kmer_start + lut_word_length > s_len {
                                continue;
                            }
                            
                            // For two-stage lookup, we use the lut_word_length match position
                            // The actual word_length matching will happen in extension phase
                            // NCBI BLAST: s_BlastnExtendInitialHit does left+right extension to verify word_length
                            let extended = 0; // Will be calculated in extension phase
                            
                            // NCBI reference: na_ungapped.c:663-664
                            // diag = s_off + diag_table->diag_array_length - q_off;
                            // real_diag = diag & diag_table->diag_mask;
                            // CRITICAL: Must match NCBI calculation exactly
                            // NCBI creates diag_table per query, so use query-specific diag_array_length
                            let (q_diag_array_length, q_diag_mask) = query_diag_lengths[q_idx as usize];
                            let s_off = kmer_start as isize;
                            let q_off = q_pos_usize as isize;
                            let diag = s_off + q_diag_array_length as isize - q_off;
                            let real_diag = (diag as usize) & q_diag_mask;

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

                            // NCBI BLAST does NOT check mask_array/mask_hash at seed level
                            // NCBI reference: na_ungapped.c:671-672
                            // NCBI only checks: if (s_off_pos < last_hit) return 0;
                            // This check happens below (line 753) using hit_level_array[real_diag].last_hit
                            // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                            // REMOVED to match NCBI BLAST behavior

                            // NCBI reference: na_ungapped.c:665-666, 672
                            // last_hit = hit_level_array[real_diag].last_hit;
                            // hit_saved = hit_level_array[real_diag].flag;
                            // if (s_off_pos < last_hit) return 0;  // hit within explored area
                            // CRITICAL: Use real_diag (masked) directly as array index, matching NCBI exactly
                            // NCBI uses real_diag directly: hit_level_array[real_diag]
                            let diag_idx = real_diag; // NCBI uses real_diag directly
                            
                            // Get current diagonal state
                            let (last_hit, hit_saved) = if use_array_indexing && diag_idx < diag_array_size {
                                let diag_entry = &hit_level_array[diag_idx];
                                (diag_entry.last_hit, diag_entry.flag != 0)
                            } else if !use_array_indexing {
                                let diag_key = pack_diag_key(q_idx, diag);
                                let diag_entry = hit_level_hash.get(&diag_key).copied().unwrap_or_default();
                                (diag_entry.last_hit, diag_entry.flag != 0)
                            } else {
                                (0, false)
                            };
                            
                            // NCBI reference: na_ungapped.c:668
                            // s_off_pos = s_off + diag_table->offset;
                            // diag_table->offset = window_size (blast_extend.c:63)
                            // CRITICAL: Match NCBI exactly
                            let s_off_pos = kmer_start + TWO_HIT_WINDOW;
                            
                            // NCBI: if (s_off_pos < last_hit) return 0;
                            // Hit within explored area should be rejected
                            if s_off_pos < last_hit {
                                continue;
                            }

                            // NCBI reference: na_ungapped.c:1081-1140
                            // For two-stage lookup, verify word_length match BEFORE ungapped extension
                            // This is done by left+right extension from the lut_word_length match position
                            let (q_ext_start, s_ext_start) = if word_length > lut_word_length {
                                // NCBI BLAST: s_BlastnExtendInitialHit does left+right extension
                                // Reference: na_ungapped.c:1106-1120
                                let ext_to = word_length - lut_word_length;
                                
                                // LEFT extension (backwards from lut_word_length start)
                                // NCBI BLAST: na_ungapped.c:1101-1114
                                // Int4 ext_left = 0;
                                // Int4 s_off = s_offset;
                                // Uint1 *q = query->sequence + q_offset;
                                // Uint1 *s = subject->sequence + s_off / COMPRESSION_RATIO;
                                // for (; ext_left < MIN(ext_to, s_offset); ++ext_left) {
                                //     s_off--;
                                //     q--;
                                //     if (s_off % COMPRESSION_RATIO == 3)
                                //         s--;
                                //     if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
                                //         break;
                                // }
                                // CRITICAL: NCBI relies on MIN(ext_to, s_offset) limit, NO explicit boundary check
                                let mut ext_left = 0;
                                let mut q_left = q_pos_usize;
                                let mut s_left = kmer_start;
                                
                                // NCBI BLAST: MIN(ext_to, s_offset) - limit left extension by s_offset
                                // s_offset (kmer_start) is the maximum number of positions we can extend left
                                // NCBI reference: na_ungapped.c:1106-1114
                                // for (; ext_left < MIN(ext_to, s_offset); ++ext_left) {
                                //     s_off--;
                                //     q--;
                                //     if (s_off % COMPRESSION_RATIO == 3)
                                //         s--;
                                //     if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
                                //         break;
                                // }
                                // CRITICAL: NCBI only limits by MIN(ext_to, s_offset) - NO query bounds check
                                // NCBI doesn't check query bounds (undefined behavior if out of bounds)
                                // For Rust safety, we use unsafe indexing but trust sequences are long enough
                                // (lookup table only finds matches within valid ranges, so sequences should be long enough)
                                let max_ext_left = ext_to.min(kmer_start);
                                
                                // NCBI BLAST: Loop condition is ext_left < MIN(ext_to, s_offset)
                                // NO explicit boundary check inside loop (NCBI doesn't have it)
                                // SAFETY: We've limited by subject bounds. For query, we trust sequences are long enough
                                // (NCBI doesn't check query bounds either - it's undefined behavior if out of bounds)
                                while ext_left < max_ext_left {
                                    // NCBI BLAST: s_off--; q--;
                                    // Reference: na_ungapped.c:1107-1108
                                    q_left -= 1;
                                    s_left -= 1;
                                    // NCBI BLAST: if (base mismatch) break;
                                    // Reference: na_ungapped.c:1111-1114
                                    if unsafe { *q_seq.get_unchecked(q_left) } != unsafe { *search_seq.get_unchecked(s_left) } {
                                        break;
                                    }
                                    ext_left += 1;
                                }
                                
                                // RIGHT extension (forwards from lut_word_length end)
                                // NCBI BLAST: na_ungapped.c:1120-1136
                                // if (ext_left < ext_to) {
                                //     Int4 ext_right = 0;
                                //     s_off = s_offset + lut_word_length;
                                //     if (s_off + ext_to - ext_left > s_range) 
                                //         continue;
                                //     q = query->sequence + q_offset + lut_word_length;
                                //     s = subject->sequence + s_off / COMPRESSION_RATIO;
                                //     for (; ext_right < ext_to - ext_left; ++ext_right) {
                                //         if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
                                //             break;
                                //         s_off++;
                                //         q++;
                                //         if (s_off % COMPRESSION_RATIO == 0)
                                //             s++;
                                //     }
                                // }
                                // Only do right extension if left didn't get all bases
                                // NCBI reference: na_ungapped.c:1120-1136
                                // if (ext_left < ext_to) {
                                //     Int4 ext_right = 0;
                                //     s_off = s_offset + lut_word_length;
                                //     if (s_off + ext_to - ext_left > s_range) 
                                //         continue;
                                //     q = query->sequence + q_offset + lut_word_length;
                                //     s = subject->sequence + s_off / COMPRESSION_RATIO;
                                //     for (; ext_right < ext_to - ext_left; ++ext_right) {
                                //         if (mismatch) break;
                                //         s_off++;
                                //         q++;
                                //     }
                                // }
                                if ext_left < ext_to {
                                    // NCBI BLAST: s_off = s_offset + lut_word_length;
                                    // NCBI BLAST: if (s_off + ext_to - ext_left > s_range) continue;
                                    // Reference: na_ungapped.c:1122-1124
                                    // CRITICAL: NCBI only checks subject bounds, NOT query bounds
                                    let s_off = kmer_start + lut_word_length;
                                    if s_off + (ext_to - ext_left) > s_len {
                                        // Not enough room in subject, skip this seed (NCBI behavior)
                                        continue;
                                    }
                                    
                                    let mut ext_right = 0;
                                    let q_off = q_pos_usize + lut_word_length;
                                    let mut q_right = q_off;
                                    let mut s_right = s_off;
                                    
                                    // NCBI BLAST: for (; ext_right < ext_to - ext_left; ++ext_right)
                                    // Reference: na_ungapped.c:1128-1136
                                    // NO bounds checks in loop condition (NCBI relies on pre-loop subject check only)
                                    // NCBI doesn't check query bounds - it would access out of bounds (undefined behavior)
                                    // For Rust safety, we use unsafe indexing but trust sequences are long enough
                                    // (lookup table only finds matches within valid ranges, so sequences should be long enough)
                                    while ext_right < (ext_to - ext_left) {
                                        // NCBI BLAST: if (base mismatch) break;
                                        // Reference: na_ungapped.c:1129-1131
                                        // SAFETY: We've checked subject bounds above. For query, we trust sequences are long enough
                                        // (NCBI doesn't check query bounds either - it's undefined behavior if out of bounds)
                                        if unsafe { *q_seq.get_unchecked(q_right) } != unsafe { *search_seq.get_unchecked(s_right) } {
                                            break;
                                        }
                                        ext_right += 1;
                                        q_right += 1;
                                        s_right += 1;
                                    }
                                    
                                    // NCBI BLAST: if (ext_left + ext_right < ext_to) continue;
                                    // Reference: na_ungapped.c:1125
                                    if ext_left + ext_right < ext_to {
                                        // Word_length match failed, skip this seed
                                        continue;
                                    }
                                }
                                
                                // Adjust positions for type_of_word and ungapped extension
                                // NCBI BLAST: q_offset -= ext_left; s_offset -= ext_left;
                                // Reference: na_ungapped.c:1143-1144
                                (q_pos_usize - ext_left, kmer_start - ext_left)
                            } else {
                                // word_length == lut_word_length, no extension needed
                                (q_pos_usize, kmer_start)
                            };
                            
                            // NCBI reference: na_ungapped.c:674-683
                            // After word_length verification, call type_of_word for two-hit mode
                            // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
                            //     word_type = s_TypeOfWord(...);
                            //     if (!word_type) return 0;
                            //     s_end += extended;
                            //     s_end_pos += extended;
                            // }
                            // NCBI reference: na_ungapped.c:656-658
                            // Boolean two_hits = (window_size > 0);
                            // Boolean off_found = FALSE;
                            // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                            // CRITICAL: Delta must be calculated at function start, BEFORE two_hits block
                            // This matches NCBI implementation exactly (na_ungapped.c:658)
                            let two_hits = TWO_HIT_WINDOW > 0;
                            let window_size = TWO_HIT_WINDOW;
                            // NCBI: Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                            let mut delta = scan_range.min(window_size.saturating_sub(word_length)) as isize;
                            let q_off = q_ext_start;
                            let s_off = s_ext_start;
                            // NCBI reference: na_ungapped.c:667-669
                            // s_end = s_off + word_length;
                            // s_off_pos = s_off + diag_table->offset;
                            // s_end_pos = s_end + diag_table->offset;
                            // diag_table->offset = window_size (blast_extend.c:63)
                            // CRITICAL: Match NCBI exactly
                            let mut s_end = s_ext_start + word_length;
                            // NCBI reference: na_ungapped.c:668
                            // s_off_pos = s_off + diag_table->offset;
                            // diag_table->offset = window_size (blast_extend.c:63)
                            // CRITICAL: Match NCBI exactly - s_off_pos is required for off-diagonal search
                            let s_off_pos = s_off + TWO_HIT_WINDOW;
                            let mut s_end_pos = s_end + TWO_HIT_WINDOW;
                            let mut word_type = 1u8; // Default: single word (when word_length == lut_word_length)
                            let mut extended = 0usize;
                            let mut off_found = false;
                            let mut hit_ready = true;
                            
                            // NCBI: if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
                            if two_hits && (hit_saved || s_end_pos > last_hit + TWO_HIT_WINDOW) {
                                // NCBI reference: na_ungapped.c:677-680
                                // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
                                //                          query_mask, query_info, s_range, 
                                //                          word_length, lut_word_length, lut, TRUE, &extended);
                                let query_mask = &query_masks_ref[q_idx as usize];
                                let (wt, ext) = type_of_word(
                                    q_seq,
                                    search_seq,
                                    q_off,
                                    s_off,
                                    query_mask,
                                    word_length,
                                    lut_word_length,
                                    true, // check_double = TRUE
                                );
                                word_type = wt;
                                extended = ext;
                                
                                // NCBI: if (!word_type) return 0;
                                if word_type == 0 {
                                    // Non-word, skip this hit
                                    continue;
                                }
                                
                                // NCBI: s_end += extended;
                                // NCBI: s_end_pos += extended;
                                s_end += extended;
                                s_end_pos += extended;
                                
                                // NCBI reference: na_ungapped.c:685-717
                                // for single word, also try off diagonals
                                // if (word_type == 1) {
                                //     /* try off-diagonals */
                                //     Int4 orig_diag = real_diag + diag_table->diag_array_length;
                                //     Int4 s_a = s_off_pos + word_length - window_size;
                                //     Int4 s_b = s_end_pos - 2 * word_length;
                                //     Int4 delta;
                                //     if (Delta < 0) Delta = 0;
                                //     for (delta = 1; delta <= Delta ; ++delta) {
                                //         // Check diag + delta and diag - delta
                                //         ...
                                //     }
                                //     if (!off_found) {
                                //         hit_ready = 0;
                                //     }
                                // }
                                if word_type == 1 {
                                    // NCBI reference: na_ungapped.c:692
                                    // if (Delta < 0) Delta = 0;
                                    // CRITICAL: Delta was calculated at function start (line 658)
                                    // Here we only check if it's negative and set to 0
                                    if delta < 0 {
                                        delta = 0;
                                    }
                                    let delta_max = delta;
                                    
                                    // NCBI reference: na_ungapped.c:689-690
                                    // Int4 s_a = s_off_pos + word_length - window_size;
                                    // Int4 s_b = s_end_pos - 2 * word_length;
                                    // CRITICAL: NCBI uses signed arithmetic (Int4), so s_a and s_b can be negative
                                    // LOSAT must use signed arithmetic to match NCBI behavior
                                    let s_a = s_off_pos as isize + word_length as isize - window_size as isize;
                                    let s_b = s_end_pos as isize - 2 * word_length as isize;
                                    
                                    // NCBI reference: na_ungapped.c:688, 694
                                    // Int4 orig_diag = real_diag + diag_table->diag_array_length;
                                    // Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
                                    // CRITICAL: Match NCBI exactly
                                    // NCBI creates diag_table per query, so use query-specific diag_array_length
                                    let (q_diag_array_length, q_diag_mask) = query_diag_lengths[q_idx as usize];
                                    let orig_diag = real_diag as isize + q_diag_array_length as isize;
                                    
                                    // NCBI reference: na_ungapped.c:693
                                    // for (delta = 1; delta <= Delta ; ++delta) {
                                    for delta in 1..=delta_max {
                                        // NCBI reference: na_ungapped.c:694-702
                                        // Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
                                        // Int4 off_s_end = hit_level_array[off_diag].last_hit;
                                        // Int4 off_s_l   = diag_table->hit_len_array[off_diag];
                                        // if ( off_s_l
                                        //  && off_s_end - delta >= s_a 
                                        //  && off_s_end - off_s_l <= s_b) {
                                        //     off_found = TRUE;
                                        //     break;
                                        // }
                                        // CRITICAL: Match NCBI exactly - off_diag is masked and used directly as array index
                                        // NCBI reference: na_ungapped.c:694
                                        // Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
                                        // Int4 off_s_end = hit_level_array[off_diag].last_hit;
                                        // CRITICAL: off_diag is already masked, use directly as array index (no diag_offset)
                                        let off_diag = ((orig_diag + delta as isize) as usize) & q_diag_mask;
                                        // NCBI reference: na_ungapped.c:695-696
                                        // Int4 off_s_end = hit_level_array[off_diag].last_hit;
                                        // Int4 off_s_l   = diag_table->hit_len_array[off_diag];
                                        // CRITICAL: off_diag is already masked, use directly as array index
                                        let (off_s_end, off_s_l) = if use_array_indexing && off_diag < diag_array_size {
                                            let off_entry = &hit_level_array[off_diag];
                                            (off_entry.last_hit, hit_len_array[off_diag])
                                        } else if !use_array_indexing {
                                            // For hash indexing, we need to use the original diag value
                                            // But NCBI doesn't use hash indexing - it always uses array
                                            // This is a fallback for LOSAT when array is too large
                                            let off_diag_key = pack_diag_key(q_idx, off_diag as isize);
                                            let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                            let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                            (off_entry.last_hit, off_len)
                                        } else {
                                            (0, 0)
                                        };
                                        
                                        // NCBI: off_s_end - delta >= s_a (signed comparison)
                                        // Convert to signed for comparison to match NCBI behavior
                                        if off_s_l > 0 
                                            && (off_s_end as isize - delta as isize) >= s_a
                                            && (off_s_end as isize - off_s_l as isize) <= s_b {
                                            off_found = true;
                                            break;
                                        }
                                        
                                        // NCBI reference: na_ungapped.c:703-711
                                        // off_diag  = (orig_diag - delta) & diag_table->diag_mask;
                                        // off_s_end = hit_level_array[off_diag].last_hit;
                                        // off_s_l   = diag_table->hit_len_array[off_diag];
                                        // if ( off_s_l
                                        //  && off_s_end >= s_a 
                                        //  && off_s_end - off_s_l + delta <= s_b) {
                                        //     off_found = TRUE;
                                        //     break;
                                        // }
                                        // CRITICAL: Match NCBI exactly - off_diag is masked and used directly as array index
                                        // NCBI creates diag_table per query, so use query-specific diag_mask
                                        let off_diag = ((orig_diag - delta as isize) as usize) & q_diag_mask;
                                        // NCBI reference: na_ungapped.c:704-705
                                        // Int4 off_s_end = hit_level_array[off_diag].last_hit;
                                        // Int4 off_s_l   = diag_table->hit_len_array[off_diag];
                                        // CRITICAL: off_diag is already masked, use directly as array index
                                        let (off_s_end, off_s_l) = if use_array_indexing && off_diag < diag_array_size {
                                            let off_entry = &hit_level_array[off_diag];
                                            (off_entry.last_hit, hit_len_array[off_diag])
                                        } else if !use_array_indexing {
                                            // For hash indexing, we need to use the original diag value
                                            // But NCBI doesn't use hash indexing - it always uses array
                                            // This is a fallback for LOSAT when array is too large
                                            let off_diag_key = pack_diag_key(q_idx, off_diag as isize);
                                            let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                            let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                            (off_entry.last_hit, off_len)
                                        } else {
                                            (0, 0)
                                        };
                                        
                                        // NCBI: off_s_end >= s_a (signed comparison)
                                        // Convert to signed for comparison to match NCBI behavior
                                        if off_s_l > 0 
                                            && (off_s_end as isize) >= s_a
                                            && (off_s_end as isize - off_s_l as isize + delta as isize) <= s_b {
                                            off_found = true;
                                            break;
                                        }
                                    }
                                    
                                    // NCBI: if (!off_found) {
                                    //     hit_ready = 0;
                                    // }
                                    if !off_found {
                                        // NCBI reference: na_ungapped.c:713-716
                                        // if (!off_found) {
                                        //     /* This is a new hit */
                                        //     hit_ready = 0;
                                        // }
                                        hit_ready = false;
                                    }
                                }
                            } else {
                                // NCBI reference: na_ungapped.c:718-726
                                // else if (check_masks) {
                                //     /* check the masks for the word */
                                //     if(!s_TypeOfWord(query, subject, &q_off, &s_off,
                                //                     query_mask, query_info, s_range, 
                                //                     word_length, lut_word_length, lut, FALSE, &extended)) return 0;
                                //     /* update the right end*/
                                //     s_end += extended;
                                //     s_end_pos += extended;
                                // }
                                // In NCBI, check_masks is TRUE by default (only FALSE when lut->stride is true)
                                // For LOSAT, we always check masks (equivalent to check_masks = TRUE)
                                let query_mask = &query_masks_ref[q_idx as usize];
                                let (wt, ext) = type_of_word(
                                    q_seq,
                                    search_seq,
                                    q_off,
                                    s_off,
                                    query_mask,
                                    word_length,
                                    lut_word_length,
                                    false, // check_double = FALSE (not in two-hit block)
                                );
                                if wt == 0 {
                                    // Non-word, skip this hit
                                    continue;
                                }
                                // NCBI: s_end += extended;
                                // NCBI: s_end_pos += extended;
                                s_end += ext;
                                s_end_pos += ext;
                                // hit_ready remains true (default) - extension will proceed
                            }
                            
                            // Now do ungapped extension from the adjusted position
                            // NCBI BLAST: s_BlastnExtendInitialHit calls ungapped extension after word_length verification
                            // Use q_off and s_off (adjusted by type_of_word if called)
                            let (qs, qe, ss, _se, ungapped_score) = extend_hit_ungapped(
                                q_seq,
                                search_seq,
                                q_off,
                                s_off,
                                reward,
                                penalty,
                                None, // Use default X-drop
                            );

                            // NCBI reference: na_ungapped.c:752-772
                            // if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
                            //     s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                            // } else {
                            //     hit_ready = 0;
                            // }
                            // hit_level_array[real_diag].last_hit = s_end_pos;
                            // hit_level_array[real_diag].flag = hit_ready;
                            // diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;
                            // CRITICAL: NCBI updates hit_level_array ALWAYS at the END, after all processing
                            // If extension fails, hit_ready = 0 but s_end_pos is NOT updated (remains at pre-extension value)
                            // If extension succeeds, s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset
                            let mut final_s_end_pos = s_end_pos; // Default: pre-extension position
                            if !(off_found || ungapped_score >= cutoff_score) {
                                // NCBI reference: na_ungapped.c:759-761
                                // If extension is skipped, hit_ready = 0, but s_end_pos is NOT updated
                                hit_ready = false;
                                dbg_ungapped_low += 1;
                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos_usize, kmer_start, ungapped_score, min_ungapped_score);
                                }
                                // final_s_end_pos remains at s_end_pos (pre-extension position)
                                // hit_level_array will be updated immediately after this block (before gapped extension)
                            } else {
                                // NCBI reference: na_ungapped.c:757-758
                                // s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                                // ungapped_data->s_start = ss (subject start after ungapped extension)
                                // ungapped_data->length = (qe - qs) (ungapped extension length)
                                // diag_table->offset = window_size (blast_extend.c:63)
                                final_s_end_pos = ss + (qe - qs) + TWO_HIT_WINDOW;
                                hit_ready = true; // Extension was successful

                                dbg_gapped_attempted += 1;

                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> GAPPED EXTENSION (ungapped_score={}, seed_len={})", q_pos_usize, kmer_start, ungapped_score, qe - qs);
                                }
                            }

                            // NCBI reference: na_ungapped.c:768-772
                            // hit_level_array[real_diag].last_hit = s_end_pos;
                            // hit_level_array[real_diag].flag = hit_ready;
                            // if (two_hits) {
                            //     diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;
                            // }
                            // CRITICAL: NCBI updates hit_level_array IMMEDIATELY after ungapped extension,
                            // BEFORE gapped extension (na_ungapped.c:768 is executed before function returns)
                            // This is different from LOSAT's previous implementation where it was updated after gapped extension
                            if use_array_indexing && diag_idx < diag_array_size {
                                hit_level_array[diag_idx].last_hit = final_s_end_pos;
                                hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                                if TWO_HIT_WINDOW > 0 {
                                    // NCBI reference: na_ungapped.c:770-771
                                    // Only update hit_len_array if two_hits (window_size > 0)
                                    hit_len_array[diag_idx] = if hit_ready { 0 } else { final_s_end_pos - s_off_pos };
                                }
                            } else if !use_array_indexing {
                                let diag_key = pack_diag_key(q_idx, diag);
                                hit_level_hash.insert(
                                    diag_key,
                                    DiagStruct {
                                        last_hit: final_s_end_pos,
                                        flag: if hit_ready { 1 } else { 0 },
                                    },
                                );
                                if TWO_HIT_WINDOW > 0 {
                                    hit_len_hash.insert(diag_key, if hit_ready { 0 } else { final_s_end_pos - s_off_pos });
                                }
                            }

                            // NCBI reference: na_ungapped.c:756-758
                            // BLAST_SaveInitialHit saves ungapped hit to init_hitlist
                            // Gapped extension happens later in blast_gapalign.c (separate function call)
                            // NCBI reference: na_ungapped.c:768-772
                            // hit_level_array is updated AFTER ungapped extension, BEFORE gapped extension
                            // LOSAT performs gapped extension inline (structural difference from NCBI),
                            // but execution order matches NCBI: hit_level_array update → gapped extension
                            if hit_ready {
                                let q_record = &queries[q_idx as usize];
                                // Gapped extension with task-specific X-drop
                                // NCBI BLAST: blastn uses 30 (non-greedy), megablast uses 25 (greedy)
                                // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:122-148
                                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c: gap_align->gap_x_dropoff = ext_params->gap_x_dropoff;
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
                                    x_drop_gapped, // Task-specific: blastn=30, megablast=25
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

                                let (bit_score, eval) = calculate_evalue(
                                    score,
                                    q_record.seq().len(),
                                    db_len_total,
                                    db_num_seqs,
                                    &params_for_closure,
                                );

                                if eval <= evalue_threshold {
                                    // Identity is matches / alignment_length, capped at 100%
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
                                    let (hit_s_start, hit_s_end) = if is_minus_strand {
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
                                        q_idx: 0,
                                        s_idx: 0,
                                        raw_score: score,
                                    });
                                }
                            }
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
                            pv_direct_lookup_ref.map(|pv_dl| pv_dl.get_hits_checked(current_kmer)).unwrap_or(&[])
                        } else {
                            // Use hash-based lookup for larger word sizes
                            hash_lookup_ref.and_then(|hl| hl.get(&current_kmer).map(|v| v.as_slice())).unwrap_or(&[])
                        };

                        for &(q_idx, q_pos) in matches_slice {
                            dbg_seeds_found += 1;
                            
                            let q_record = &queries[q_idx as usize];
                            
                            // NCBI reference: na_ungapped.c:730-731
                            // Int4 context = BSearchContextInfo(q_off, query_info);
                            // cutoffs = word_params->cutoffs + context;
                            // For blastn, we need to determine context from q_idx
                            // Each query has 2 contexts: 0 (forward) and 1 (reverse)
                            let context = (q_idx as usize) * 2; // Forward strand context
                            let cutoff_score = if context < cutoff_scores_ref.len() {
                                per_subject_cutoff_scores[context]
                            } else {
                                // Fallback: compute cutoff_score if context is out of bounds
                                let query_len = q_record.seq().len() as i64;
                                let subject_len = s_len as i64;
                                compute_blastn_cutoff_score(
                                    query_len,
                                    subject_len,
                                    evalue_threshold,
                                    GAP_TRIGGER_BIT_SCORE_NUCL,
                                    &params_for_closure,
                                    &params_for_closure,
                                    1.0,
                                )
                            };
                            
                        // NCBI reference: na_ungapped.c:663-664
                        // diag = s_off + diag_table->diag_array_length - q_off;
                        // real_diag = diag & diag_table->diag_mask;
                        // CRITICAL: Must match NCBI calculation exactly
                        // NCBI creates diag_table per query, so use query-specific diag_array_length
                        let (q_diag_array_length, q_diag_mask) = query_diag_lengths[q_idx as usize];
                        let s_off = kmer_start as isize;
                        let q_off = q_pos as isize;
                        let diag = s_off + q_diag_array_length as isize - q_off;
                        let real_diag = (diag as usize) & q_diag_mask;

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

                        // NCBI BLAST does NOT check mask_array/mask_hash at seed level
                        // NCBI reference: na_ungapped.c:671-672
                        // NCBI only checks: if (s_off_pos < last_hit) return 0;
                        // This check happens below using hit_level_array[real_diag].last_hit
                        // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                        // REMOVED to match NCBI BLAST behavior

                        // NCBI reference: na_ungapped.c:665-666, 672
                        // last_hit = hit_level_array[real_diag].last_hit;
                        // hit_saved = hit_level_array[real_diag].flag;
                        // if (s_off_pos < last_hit) return 0;  // hit within explored area
                        // CRITICAL: Use real_diag (masked) directly as array index, matching NCBI exactly
                        // NCBI uses real_diag directly: hit_level_array[real_diag]
                        let diag_idx = real_diag; // NCBI uses real_diag directly
                        
                        // Get current diagonal state
                        let (last_hit, hit_saved) = if use_array_indexing && diag_idx < diag_array_size {
                            let diag_entry = &hit_level_array[diag_idx];
                            (diag_entry.last_hit, diag_entry.flag != 0)
                        } else if !use_array_indexing {
                            let diag_key = pack_diag_key(q_idx, diag);
                            let diag_entry = hit_level_hash.get(&diag_key).copied().unwrap_or_default();
                            (diag_entry.last_hit, diag_entry.flag != 0)
                        } else {
                            (0, false)
                        };
                        
                        // NCBI reference: na_ungapped.c:668
                        // s_off_pos = s_off + diag_table->offset;
                        // diag_table->offset = window_size (blast_extend.c:63)
                        // CRITICAL: Match NCBI exactly
                        let s_off_pos = kmer_start + TWO_HIT_WINDOW;
                        
                        // NCBI: if (s_off_pos < last_hit) return 0;
                        // Hit within explored area should be rejected
                        if s_off_pos < last_hit {
                            continue;
                        }
                        
                        // NCBI reference: na_ungapped.c:674-683
                        // After two-hit check, call type_of_word for two-hit mode
                        // if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
                        //     word_type = s_TypeOfWord(...);
                        //     if (!word_type) return 0;
                        //     s_end += extended;
                        //     s_end_pos += extended;
                        // }
                        // NCBI reference: na_ungapped.c:667-669
                        // s_end = s_off + word_length;
                        // s_off_pos = s_off + diag_table->offset;
                        // s_end_pos = s_end + diag_table->offset;
                        // diag_table->offset = window_size (blast_extend.c:63)
                        // CRITICAL: Match NCBI exactly
                        // NCBI reference: na_ungapped.c:656-658
                        // Boolean two_hits = (window_size > 0);
                        // Boolean off_found = FALSE;
                        // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                        // CRITICAL: Delta must be calculated at function start, BEFORE two_hits block
                        // This matches NCBI implementation exactly (na_ungapped.c:658)
                        let two_hits = TWO_HIT_WINDOW > 0;
                        let window_size = TWO_HIT_WINDOW;
                        // NCBI: Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                        let mut delta = scan_range.min(window_size.saturating_sub(safe_k)) as isize;
                        let q_off = q_pos as usize;
                        let s_off = kmer_start;
                        let mut s_end = kmer_start + safe_k;
                        // NCBI reference: na_ungapped.c:668
                        // s_off_pos = s_off + diag_table->offset;
                        // diag_table->offset = window_size (blast_extend.c:63)
                        // CRITICAL: Match NCBI exactly - s_off_pos is required for off-diagonal search
                        let s_off_pos = s_off + TWO_HIT_WINDOW;
                        let mut s_end_pos = s_end + TWO_HIT_WINDOW;
                        let mut word_type = 1u8; // Default: single word (when word_length == lut_word_length)
                        let mut extended = 0usize;
                        let mut off_found = false;
                        let mut hit_ready = true;
                        
                        let q_record = &queries[q_idx as usize];
                        
                        // NCBI: if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) {
                        if two_hits && (hit_saved || s_end_pos > last_hit + TWO_HIT_WINDOW) {
                            // NCBI reference: na_ungapped.c:677-680
                            // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
                            //                          query_mask, query_info, s_range, 
                            //                          word_length, lut_word_length, lut, TRUE, &extended);
                            // For non-two-stage lookup, word_length == lut_word_length, so type_of_word returns (1, 0)
                            let query_mask = &query_masks_ref[q_idx as usize];
                            let (wt, ext) = type_of_word(
                                q_record.seq(),
                                search_seq,
                                q_off,
                                s_off,
                                query_mask,
                                safe_k, // word_length
                                safe_k, // lut_word_length (same for non-two-stage)
                                true, // check_double = TRUE
                            );
                            word_type = wt;
                            extended = ext;
                            
                            // NCBI: if (!word_type) return 0;
                            if word_type == 0 {
                                // Non-word, skip this hit
                                continue;
                            }
                            
                            // NCBI: s_end += extended;
                            // NCBI: s_end_pos += extended;
                            s_end += extended;
                            s_end_pos += extended;
                            
                            // NCBI reference: na_ungapped.c:852-881 - off-diagonal search (DiagHash version)
                            if word_type == 1 {
                                // NCBI reference: na_ungapped.c:858
                                // if (Delta < 0) Delta = 0;
                                // CRITICAL: Delta was calculated at function start (line 825 for Hash version, line 658 for Array version)
                                // Here we only check if it's negative and set to 0
                                if delta < 0 {
                                    delta = 0;
                                }
                                let delta_max = delta;
                                
                                // NCBI reference: na_ungapped.c:855-856
                                // Int4 s_a = s_off_pos + word_length - window_size;
                                // Int4 s_b = s_end_pos - 2 * word_length;
                                // CRITICAL: NCBI uses signed arithmetic (Int4), so s_a and s_b can be negative
                                // LOSAT must use signed arithmetic to match NCBI behavior
                                let s_a = s_off_pos as isize + safe_k as isize - window_size as isize;
                                let s_b = s_end_pos as isize - 2 * safe_k as isize;
                                
                                // NCBI reference: na_ungapped.c:688, 694
                                // Int4 orig_diag = real_diag + diag_table->diag_array_length;
                                // Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
                                // CRITICAL: Match NCBI exactly
                                // NCBI creates diag_table per query, so use query-specific diag_array_length
                                let (q_diag_array_length, q_diag_mask) = query_diag_lengths[q_idx as usize];
                                let orig_diag = real_diag as isize + q_diag_array_length as isize;
                                
                                // NCBI reference: na_ungapped.c:859
                                // for (delta = 1; delta <= Delta; ++delta) {
                                for delta in 1..=delta_max {
                                    // NCBI reference: na_ungapped.c:860-871
                                    // Int4 off_s_end = 0;
                                    // Int4 off_s_l = 0;
                                    // Int4 off_hit_saved = 0;
                                    // Int4 off_rc = s_BlastDiagHashRetrieve(hash_table, diag + delta, 
                                    //               &off_s_end, &off_s_l, &off_hit_saved);
                                    // if ( off_rc
                                    //   && off_s_l
                                    //   && off_s_end - delta >= s_a
                                    //   && off_s_end - off_s_l <= s_b) {
                                    //      off_found = TRUE;
                                    //      break;
                                    // }
                                    // CRITICAL: Match NCBI exactly - off_diag is masked and used directly as array index
                                    // NCBI creates diag_table per query, so use query-specific diag_mask
                                    let off_diag = ((orig_diag + delta as isize) as usize) & q_diag_mask;
                                    // NCBI reference: na_ungapped.c:695-696
                                    // Int4 off_s_end = hit_level_array[off_diag].last_hit;
                                    // Int4 off_s_l   = diag_table->hit_len_array[off_diag];
                                    // CRITICAL: off_diag is already masked, use directly as array index
                                    let (off_s_end, off_s_l) = if use_array_indexing && off_diag < diag_array_size {
                                        let off_entry = &hit_level_array[off_diag];
                                        (off_entry.last_hit, hit_len_array[off_diag])
                                    } else if !use_array_indexing {
                                        // For hash indexing, we need to use the original diag value
                                        // But NCBI doesn't use hash indexing - it always uses array
                                        // This is a fallback for LOSAT when array is too large
                                        let off_diag_key = pack_diag_key(q_idx, off_diag as isize);
                                        let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                        let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                        (off_entry.last_hit, off_len)
                                    } else {
                                        (0, 0)
                                    };
                                    
                                    // NCBI: off_s_end - delta >= s_a (signed comparison)
                                    // Convert to signed for comparison to match NCBI behavior
                                    if off_s_l > 0 
                                        && (off_s_end as isize - delta as isize) >= s_a
                                        && (off_s_end as isize - off_s_l as isize) <= s_b {
                                        off_found = true;
                                        break;
                                    }
                                    
                                    // NCBI reference: na_ungapped.c:872-880
                                    // off_rc = s_BlastDiagHashRetrieve(hash_table, diag - delta, 
                                    //               &off_s_end, &off_s_l, &off_hit_saved);
                                    // if ( off_rc
                                    //   && off_s_l
                                    //   && off_s_end >= s_a
                                    //   && off_s_end - off_s_l + delta <= s_b) {
                                    //      off_found = TRUE;
                                    //      break;
                                    // }
                                    // CRITICAL: Match NCBI exactly - off_diag is masked and used directly as array index
                                    // NCBI creates diag_table per query, so use query-specific diag_mask
                                    let off_diag = ((orig_diag - delta as isize) as usize) & q_diag_mask;
                                    // NCBI reference: na_ungapped.c:704-705
                                    // Int4 off_s_end = hit_level_array[off_diag].last_hit;
                                    // Int4 off_s_l   = diag_table->hit_len_array[off_diag];
                                    // CRITICAL: off_diag is already masked, use directly as array index
                                    let (off_s_end, off_s_l) = if use_array_indexing && off_diag < diag_array_size {
                                        let off_entry = &hit_level_array[off_diag];
                                        (off_entry.last_hit, hit_len_array[off_diag])
                                    } else if !use_array_indexing {
                                        // For hash indexing, we need to use the original diag value
                                        // But NCBI doesn't use hash indexing - it always uses array
                                        // This is a fallback for LOSAT when array is too large
                                        let off_diag_key = pack_diag_key(q_idx, off_diag as isize);
                                        let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                        let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                        (off_entry.last_hit, off_len)
                                    } else {
                                        (0, 0)
                                    };
                                    
                                    // NCBI: off_s_end >= s_a (signed comparison)
                                    // Convert to signed for comparison to match NCBI behavior
                                    if off_s_l > 0 
                                        && (off_s_end as isize) >= s_a
                                        && (off_s_end as isize - off_s_l as isize + delta as isize) <= s_b {
                                        off_found = true;
                                        break;
                                    }
                                }
                                
                                if !off_found {
                                    // NCBI reference: na_ungapped.c:713-716
                                    // if (!off_found) {
                                    //     /* This is a new hit */
                                    //     hit_ready = 0;
                                    // }
                                    hit_ready = false;
                                }
                            }
                        } else {
                            // NCBI reference: na_ungapped.c:718-726
                            // else if (check_masks) {
                            //     /* check the masks for the word */
                            //     if(!s_TypeOfWord(query, subject, &q_off, &s_off,
                            //                     query_mask, query_info, s_range, 
                            //                     word_length, lut_word_length, lut, FALSE, &extended)) return 0;
                            //     /* update the right end*/
                            //     s_end += extended;
                            //     s_end_pos += extended;
                            // }
                            // In NCBI, check_masks is TRUE by default (only FALSE when lut->stride is true)
                            // For LOSAT, we always check masks (equivalent to check_masks = TRUE)
                            let q_record = &queries[q_idx as usize];
                            let query_mask = &query_masks_ref[q_idx as usize];
                            let (wt, ext) = type_of_word(
                                q_record.seq(),
                                search_seq,
                                q_off,
                                s_off,
                                query_mask,
                                safe_k, // word_length
                                safe_k, // lut_word_length (same for non-two-stage)
                                false, // check_double = FALSE (not in two-hit block)
                            );
                            if wt == 0 {
                                // Non-word, skip this hit
                                continue;
                            }
                            // NCBI: s_end += extended;
                            // NCBI: s_end_pos += extended;
                            s_end += ext;
                            s_end_pos += ext;
                            // hit_ready remains true (default) - extension will proceed
                        }

                        // Now do ungapped extension from the adjusted position
                        // Use q_off and s_off (adjusted by type_of_word if called)
                        let (qs, qe, ss, _se, ungapped_score) = extend_hit_ungapped(
                            q_record.seq(),
                            search_seq,
                            q_off,
                            s_off,
                            reward,
                            penalty,
                            None, // Use default X-drop
                        );

                        // NCBI reference: na_ungapped.c:752-772
                        // if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
                        //     s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                        // } else {
                        //     hit_ready = 0;
                        // }
                        // hit_level_array[real_diag].last_hit = s_end_pos;
                        // hit_level_array[real_diag].flag = hit_ready;
                        // diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;
                        // CRITICAL: NCBI updates hit_level_array IMMEDIATELY after ungapped extension,
                        // BEFORE gapped extension (na_ungapped.c:768 is executed before function returns)
                        let mut final_s_end_pos = s_end_pos; // Default: pre-extension position
                        if !(off_found || ungapped_score >= cutoff_score) {
                            // NCBI reference: na_ungapped.c:759-761
                            // If extension is skipped, hit_ready = 0, but s_end_pos is NOT updated
                            hit_ready = false;
                            dbg_ungapped_low += 1;
                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos, kmer_start, ungapped_score, min_ungapped_score);
                            }
                            // final_s_end_pos remains at s_end_pos (pre-extension position)
                            // hit_level_array will be updated immediately after this block (before gapped extension)
                        } else {
                            // NCBI reference: na_ungapped.c:757-758
                            // s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                            // ungapped_data->s_start = ss (subject start after ungapped extension)
                            // ungapped_data->length = (qe - qs) (ungapped extension length)
                            // diag_table->offset = window_size (blast_extend.c:63)
                            final_s_end_pos = ss + (qe - qs) + TWO_HIT_WINDOW;
                            hit_ready = true; // Extension was successful

                            dbg_gapped_attempted += 1;

                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> GAPPED EXTENSION (ungapped_score={}, seed_len={})", q_pos, kmer_start, ungapped_score, qe - qs);
                            }
                        }

                        // NCBI reference: na_ungapped.c:768-772
                        // hit_level_array[real_diag].last_hit = s_end_pos;
                        // hit_level_array[real_diag].flag = hit_ready;
                        // if (two_hits) {
                        //     diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;
                        // }
                        // CRITICAL: NCBI updates hit_level_array IMMEDIATELY after ungapped extension,
                        // BEFORE gapped extension (na_ungapped.c:768 is executed before function returns)
                        // This is different from LOSAT's previous implementation where it was updated after gapped extension
                        if use_array_indexing && diag_idx < diag_array_size {
                            hit_level_array[diag_idx].last_hit = final_s_end_pos;
                            hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                            if TWO_HIT_WINDOW > 0 {
                                // NCBI reference: na_ungapped.c:770-771
                                // Only update hit_len_array if two_hits (window_size > 0)
                                hit_len_array[diag_idx] = if hit_ready { 0 } else { final_s_end_pos - s_off_pos };
                            }
                        } else if !use_array_indexing {
                            let diag_key = pack_diag_key(q_idx, diag);
                            hit_level_hash.insert(
                                diag_key,
                                DiagStruct {
                                    last_hit: final_s_end_pos,
                                    flag: if hit_ready { 1 } else { 0 },
                                },
                            );
                            if TWO_HIT_WINDOW > 0 {
                                hit_len_hash.insert(diag_key, if hit_ready { 0 } else { final_s_end_pos - s_off_pos });
                            }
                        }

                        // NCBI reference: na_ungapped.c:756-758
                        // BLAST_SaveInitialHit saves ungapped hit to init_hitlist
                        // Gapped extension happens later in blast_gapalign.c (separate function call)
                        // NCBI reference: na_ungapped.c:768-772
                        // hit_level_array is updated AFTER ungapped extension, BEFORE gapped extension
                        // LOSAT performs gapped extension inline (structural difference from NCBI),
                        // but execution order matches NCBI: hit_level_array update → gapped extension
                        if hit_ready {
                            let q_record = &queries[q_idx as usize];
                            // Gapped extension with task-specific X-drop
                            // NCBI BLAST: blastn uses 30 (non-greedy), megablast uses 25 (greedy)
                            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:122-148
                            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c: gap_align->gap_x_dropoff = ext_params->gap_x_dropoff;
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
                                x_drop_gapped, // Task-specific: blastn=30, megablast=25
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
                                &params_for_closure,
                            );
                            
                            // Skip if e-value is too high
                            if eval > evalue_threshold {
                                // hit_level_array was already updated above, regardless of e-value filtering
                            } else {
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
                                    q_idx: 0,
                                    s_idx: 0,
                                    raw_score: score,
                                });
                            }
                        }
                    } // end of s_pos loop for original lookup
                } // end of if two_stage_lookup_ref.is_none()
                } // end of if let Some(two_stage)
            } // end of strand loop

            // Print debug summary for this subject
            if debug_mode {
                eprintln!(
                    "[DEBUG] Subject {}: positions={}, seeds_found={}, ungapped_low={}, two_hit_failed={}, gapped_attempted={}, window_seeds={}, hits={}",
                    s_id, dbg_total_s_positions, dbg_seeds_found, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds, local_hits.len()
                );
            }

            // Suppress unused variable warnings when not in debug mode
            let _ = (dbg_total_s_positions, dbg_ambiguous_skipped, dbg_no_lookup_match, dbg_seeds_found, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds);

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
