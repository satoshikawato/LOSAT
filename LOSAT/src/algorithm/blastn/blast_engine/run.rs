//! Main BLASTN run function
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process.

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::mpsc::channel;

use crate::common::{write_output_ncbi_order, Hit};
use crate::stats::lookup_nucl_params;
use crate::core::blast_encoding::{encode_iupac_to_blastna, encode_iupac_to_ncbi2na};

// Import from parent blastn module (super::super:: because we're in blast_engine/)
use super::super::args::BlastnArgs;
use super::super::alignment::{
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    extend_final_traceback,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_traceback_with_scratch,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    greedy_gapped_alignment_score_only,
    greedy_gapped_alignment_with_traceback,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
};
use super::super::extension::{extend_hit_ungapped, type_of_word};
use crate::utils::dust::MaskedInterval;
use super::super::constants::TWO_HIT_WINDOW;
use super::super::coordination::{configure_task, finalize_task_config, read_sequences, prepare_sequence_data, build_lookup_tables};
use super::super::lookup::{reverse_complement, pack_diag_key};
use super::super::ncbi_cutoffs::{compute_blastn_cutoff_score, GAP_TRIGGER_BIT_SCORE_NUCL};
use super::super::blast_extend::DiagStruct;
use super::super::interval_tree::{BlastIntervalTree, TreeHsp, IndexMethod};
use super::super::filtering::{
    purge_hsps_with_common_endpoints,
    purge_hsps_with_common_endpoints_ex,
    reevaluate_hsp_with_ambiguities_gapped_ex,
    ReevalParams,
};

// Import from this module (blast_engine)
use super::calculate_evalue;

/// Structure to hold ungapped hit data for batch processing
/// NCBI reference: blast_gapalign.c - init_hsp_array is sorted by score descending
/// before gapped extension
#[derive(Clone)]
struct UngappedHit {
    q_idx: u32,
    // Ungapped extension results (0-based coordinates)
    qs: usize,        // query start
    qe: usize,        // query end (exclusive)
    ss: usize,        // subject start in search_seq coordinates
    se: usize,        // subject end in search_seq coordinates
    score: i32,       // ungapped score
    is_minus_strand: bool,
    // Original sequence data for gapped extension
    q_seq_idx: usize, // index into queries array
}

/// Internal 0-based coordinates for Phase 2 interval tree processing
/// NCBI reference: blast_traceback.c:679-692 - Phase 2 uses internal coordinates
#[derive(Clone)]
struct InternalHitData {
    q_offset_0: usize,  // 0-based query start
    q_end_0: usize,     // 0-based query end
    s_offset_0: usize,  // 0-based subject start (canonical, always < s_end_0)
    s_end_0: usize,     // 0-based subject end (canonical, always > s_offset_0)
    is_minus_strand: bool,
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

    // Configure task-specific parameters (initial configuration)
    let mut config = configure_task(&args);

    // Read sequences
    let (queries, query_ids, subjects) = read_sequences(&args)?;
    if queries.is_empty() || subjects.is_empty() {
        return Ok(());
    }

    // Finalize configuration with query-dependent parameters (adaptive lookup table selection)
    // NCBI reference: blast_nalookup.c (BlastChooseNaLookupTable)
    let total_query_length: usize = queries.iter().map(|r| r.seq().len()).sum();
    finalize_task_config(&mut config, total_query_length);

    if args.verbose {
        eprintln!("[INFO] Adaptive lookup: approx_entries={}, lut_word_length={}, two_stage={}, scan_step={}",
                  total_query_length, config.lut_word_length, config.use_two_stage, config.scan_step);
    }

    // Prepare sequence data (including DUST masking)
    let seq_data = prepare_sequence_data(&args, queries, query_ids, subjects);

    // Get Karlin-Altschul parameters
    // CRITICAL: NCBI uses DIFFERENT params for ungapped and gapped calculations!
    // - Ungapped params (kbp_std): Used for gap_trigger calculation
    // - Gapped params (kbp_gap): Used for cutoff_score_max and length adjustment
    // Reference: blast_parameters.c:343-344 uses kbp_std for gap_trigger
    use crate::config::NuclScoringSpec;

    // Gapped params (for length adjustment and cutoff_score_max)
    let scoring_spec_gapped = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    let params_gapped = lookup_nucl_params(&scoring_spec_gapped);

    // Ungapped params (for gap_trigger calculation)
    // NCBI uses gap_open=0, gap_extend=0 to get ungapped params
    let scoring_spec_ungapped = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: 0,
        gap_extend: 0,
    };
    let params_ungapped = lookup_nucl_params(&scoring_spec_ungapped);

    // Use gapped params for E-value calculations (same as before)
    let params = params_gapped.clone();

    // Build lookup tables
    let (lookup_tables, scan_step) = build_lookup_tables(
        &config,
        &args,
        &seq_data.queries,
        &seq_data.query_masks,
    );

    // NCBI reference: blast_traceback.c:654 (cutoff_score_min computed per subject)
    // Cutoff scores are computed per subject in the main loop; no global map needed.

    if args.verbose {
        eprintln!("Searching...");
    }

    let bar = ProgressBar::new(seq_data.subjects.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    // Channel for sending hits
    // Use Option to signal completion: None means "all subjects processed"
    let (tx, rx) = channel::<Option<Vec<Hit>>>();
    let out_path = args.out.clone();

    // Keep a sender for the main thread to send the completion signal
    let tx_main = tx.clone();
    let verbose = args.verbose;

    let writer_handle = std::thread::spawn(move || -> Result<()> {
        if verbose {
            eprintln!("[INFO] Writer thread started, waiting for hits...");
        }
        let mut all_hits = Vec::new();
        let mut messages_received = 0usize;

        while let Ok(msg) = rx.recv() {
            match msg {
                Some(hits) => {
                    messages_received += 1;
                    if verbose && (messages_received == 1 || messages_received % 100 == 0) {
                        eprintln!(
                            "[INFO] Received message #{}, {} hits so far",
                            messages_received,
                            all_hits.len() + hits.len()
                        );
                    }
                    all_hits.extend(hits);
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

        // NCBI reference: blast_traceback.c:633-692 (post-gapped processing is per-subject)
        // Per-subject traceback already applies purge/reevaluation; skip duplicate global pass.
        let filtered_hits = all_hits;

        if verbose {
            eprintln!(
                "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
                chain_start.elapsed().as_secs_f64(),
                filtered_hits.len()
            );
        }
        let write_start = std::time::Instant::now();

        write_output_ncbi_order(filtered_hits, out_path.as_ref())?;

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
    // BLASTN-specific debug mode: set LOSAT_DEBUG_BLASTN=1 to enable detailed hit loss diagnostics
    let blastn_debug = std::env::var("LOSAT_DEBUG_BLASTN").is_ok();
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
    let x_drop_final = config.x_drop_final;  // Final traceback X-drop (100)
    let scan_range = config.scan_range; // For off-diagonal hit detection
    let min_diag_separation = config.min_diag_separation; // For MB_HSP_CLOSE containment check
    let db_len_total = seq_data.db_len_total;
    let db_num_seqs = seq_data.db_num_seqs;
    let params_for_closure = params.clone();  // Gapped params for E-value
    let params_ungapped_for_closure = params_ungapped.clone();  // Ungapped params for gap_trigger
    let params_gapped_for_closure = params_gapped.clone();  // Gapped params for cutoff_score_max
    let evalue_threshold = args.evalue;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:866-869 (query blastna, subject ncbi2na)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
    let encoded_queries_blastna: Vec<Vec<u8>> = queries_ref
        .iter()
        .map(|q_record| encode_iupac_to_blastna(q_record.seq()))
        .collect();

    seq_data.subjects
        .par_iter()
        .enumerate()
        .for_each_init(
            || {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
                // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct per thread)
                (tx.clone(), GapAlignScratch::new())
            },
            |state, (_s_idx, s_record)| {
            let (tx, gap_scratch) = state;
            let queries = queries_ref;
            let query_ids = query_ids_ref;
            // NCBI reference: blast_traceback.c:679-692
            // Use hits_with_internal to preserve 0-based coordinates for Phase 2
            let mut hits_with_internal: Vec<(Hit, InternalHitData)> = Vec::new();

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

            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:866-869 (query blastna, subject ncbi2na)
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:476-489 (BlastCompressBlastnaSequence old_seq[i] & 3)
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:503-507 (greedy traceback uses uncompressed subject)
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:806-833 (GetReverseNuclSequence uses ncbi4na/blastna)
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
            let s_seq_blastna = encode_iupac_to_blastna(s_seq);
            let s_seq_rc_blastna = encode_iupac_to_blastna(&s_seq_rc);
            let (s_seq_ncbi2na, s_seq_rc_ncbi2na) = if use_dp {
                (None, None)
            } else {
                (
                    Some(encode_iupac_to_ncbi2na(s_seq)),
                    Some(encode_iupac_to_ncbi2na(&s_seq_rc)),
                )
            };

            // PERFORMANCE OPTIMIZATION: Pre-compute cutoff scores per query
            // This avoids redundant compute_blastn_cutoff_score calls in the inner loop
            // (cutoff only depends on query_len and subject_len, not on seed position)
            //
            // NCBI reference: blast_parameters.c:340-344
            // NCBI uses UNGAPPED params (kbp_std) for gap_trigger calculation
            // and GAPPED params (kbp_gap) for cutoff_score_max calculation
            let subject_len = s_len as i64;
            let cutoff_scores: Vec<i32> = queries.iter().map(|q_record| {
                let query_len = q_record.seq().len() as i64 * 2; // Both strands for blastn
                compute_blastn_cutoff_score(
                    query_len,
                    subject_len,
                    evalue_threshold,
                    GAP_TRIGGER_BIT_SCORE_NUCL,
                    &params_ungapped_for_closure,  // UNGAPPED for gap_trigger (NCBI: kbp_std)
                    &params_gapped_for_closure,    // GAPPED for cutoff_score_max (NCBI: kbp_gap)
                    1.0,
                )
            }).collect();

            // BLASTN debug: Log cutoff scores for this query-subject pair
            if blastn_debug {
                eprintln!(
                    "[BLASTN_DEBUG] Subject {} (len={}): cutoff_scores={:?}, params_ungapped=(lambda={:.4}, K={:.4}), params_gapped=(lambda={:.4}, K={:.4})",
                    s_id, s_len, cutoff_scores,
                    params_ungapped_for_closure.lambda, params_ungapped_for_closure.k,
                    params_gapped_for_closure.lambda, params_gapped_for_closure.k
                );
            }

            // NCBI reference: blast_gapalign.c:3826-3831 Blast_IntervalTreeInit
            // Initialize interval tree for HSP containment checking
            // Tree is indexed by query offsets (primary) and subject offsets (midpoint subtrees)
            let max_query_len = queries.iter().map(|q| q.seq().len()).max().unwrap_or(1);
            let mut interval_tree = BlastIntervalTree::new(
                0,                           // q_min
                (max_query_len + 1) as i32,  // q_max
                0,                           // s_min
                (s_len + 1) as i32,          // s_max
            );

            // NCBI architecture: Collect all ungapped hits first, then process in score-descending order
            // NCBI reference: blast_gapalign.c:3824 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))
            // This is critical for correct containment checking: high-score HSPs must be processed first
            let mut ungapped_hits: Vec<UngappedHit> = Vec::new();

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
                // For very large subject sequences (>6MB), fall back to HashMap
                // NCBI reference: blast_extend.c uses array-based diagonal tracking
                const MAX_ARRAY_DIAG_SIZE: usize = 12_000_000; // ~6M diagonals max (for ~6MB subject)
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

                // NCBI BLAST does NOT use a separate mask array for diagonal suppression
                // NCBI only uses last_hit in hit_level_array (checked at line 753)
                // This mask_array/mask_hash was LOSAT-specific and caused excessive filtering
                // REMOVED to match NCBI BLAST behavior

                // NCBI reference: blast_extend.h:77-80, na_ungapped.c:660-666
                // Two-hit tracking: DiagStruct array for tracking last_hit and flag per diagonal
                // hit_level_array[diag] contains {last_hit, flag} (equivalent to NCBI's DiagStruct)
                // hit_len_array[diag] contains hit length (0 = no hit, >0 = hit length)
                // For single query: use Vec for O(1) access, otherwise use HashMap
                let mut hit_level_array: Vec<DiagStruct> = if use_array_indexing {
                    vec![DiagStruct::default(); diag_array_size]
                } else {
                    Vec::new()
                };
                let mut hit_level_hash: FxHashMap<u64, DiagStruct> = if !use_array_indexing {
                    FxHashMap::default()
                } else {
                    FxHashMap::default()
                };

                // NCBI reference: na_ungapped.c:696, 705: diag_table->hit_len_array[off_diag]
                // Hit length array: stores the length of the last hit on each diagonal
                // 0 = no hit, >0 = hit length
                let mut hit_len_array: Vec<usize> = if use_array_indexing {
                    vec![0; diag_array_size]
                } else {
                    Vec::new()
                };
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

                    // DEBUG: Count loop iterations
                    let mut dbg_left_ext_iters = 0usize;
                    let mut dbg_right_ext_iters = 0usize;
                    let mut dbg_ungapped_ext_calls = 0usize;
                    let dbg_start_time = std::time::Instant::now();

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

                        // DEBUG: Track processing time
                        let _dbg_start = std::time::Instant::now();

                        // Lookup using lut_word_length k-mer
                        let matches_slice = two_stage.get_hits(current_lut_kmer);

                        // DEBUG: Log max matches
                        if matches_slice.len() > 1000 {
                            eprintln!("[WARN] Large matches_slice: len={} for kmer at position {}", matches_slice.len(), kmer_start);
                        }

                        // For each match, check if word_length match exists
                        for &(q_idx, q_pos) in matches_slice {
                            dbg_seeds_found += 1;

                            // Check if word_length match exists starting at these positions
                            // Need to verify that subject[kmer_start..kmer_start+word_length]
                            // matches query[q_pos..q_pos+word_length]
                            let q_record = &queries[q_idx as usize];
                            let q_seq = q_record.seq();
                            let q_pos_usize = q_pos as usize;

                            // Use pre-computed cutoff score (computed once per query-subject pair)
                            let cutoff_score = cutoff_scores[q_idx as usize];

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

                            // NCBI BLAST does NOT check mask_array/mask_hash at seed level
                            // NCBI reference: na_ungapped.c:671-672
                            // NCBI only checks: if (s_off_pos < last_hit) return 0;
                            // This check happens below (line 753) using hit_level_array[real_diag].last_hit
                            // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                            // REMOVED to match NCBI BLAST behavior

                            // NCBI reference: na_ungapped.c:656-672
                            // Two-hit filter: check if there's a previous hit within window_size
                            // NCBI: Boolean two_hits = (window_size > 0);
                            // NCBI: last_hit = hit_level_array[real_diag].last_hit;
                            // NCBI: hit_saved = hit_level_array[real_diag].flag;
                            // NCBI: if (s_off_pos < last_hit) return 0;  // hit within explored area
                            let diag_idx = if use_array_indexing {
                                (diag + diag_offset) as usize
                            } else {
                                0 // Not used for hash indexing
                            };

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

                            // NCBI: s_off_pos = s_off + diag_table->offset;
                            // In LOSAT, we use kmer_start directly (0-based), but need to add offset for comparison
                            let s_off_pos = kmer_start + diag_offset as usize;

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
                                // But LOSAT must also limit by q_pos to prevent underflow
                                // NCBI assumes sequences are long enough, but we need explicit bounds
                                let max_ext_left = ext_to.min(kmer_start).min(q_pos_usize);

                                // NCBI BLAST: Loop condition is ext_left < MIN(ext_to, s_offset)
                                // NO explicit boundary check inside loop (NCBI doesn't have it)
                                // SAFETY: We've limited by subject bounds. For query, we trust sequences are long enough
                                // (NCBI doesn't check query bounds either - it's undefined behavior if out of bounds)
                                while ext_left < max_ext_left {
                                    dbg_left_ext_iters += 1;
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
                                    let s_off = kmer_start + lut_word_length;
                                    if s_off + (ext_to - ext_left) > s_len {
                                        // Not enough room in subject, skip this seed (NCBI behavior)
                                        continue;
                                    }

                                    let mut ext_right = 0;
                                    let q_off = q_pos_usize + lut_word_length;
                                    // LOSAT-specific: Also check query bounds (NCBI doesn't, causing UB)
                                    if q_off + (ext_to - ext_left) > q_seq.len() {
                                        // Not enough room in query, skip this seed
                                        continue;
                                    }
                                    let mut q_right = q_off;
                                    let mut s_right = s_off;

                                    // NCBI BLAST: for (; ext_right < ext_to - ext_left; ++ext_right)
                                    // Reference: na_ungapped.c:1128-1136
                                    // NO bounds checks in loop condition (NCBI relies on pre-loop subject check only)
                                    // NCBI doesn't check query bounds - it would access out of bounds (undefined behavior)
                                    // For Rust safety, we use unsafe indexing but trust sequences are long enough
                                    // (lookup table only finds matches within valid ranges, so sequences should be long enough)
                                    while ext_right < (ext_to - ext_left) {
                                        dbg_right_ext_iters += 1;
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
                            let two_hits = TWO_HIT_WINDOW > 0;
                            let q_off = q_ext_start;
                            let s_off = s_ext_start;
                            let mut s_end = s_ext_start + word_length;
                            let mut s_end_pos = s_end + diag_offset as usize;
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
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
                                // s_TypeOfWord uses query_mask to skip masked seeds
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
                                    // NCBI reference: na_ungapped.c:658, 692
                                    // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                                    // if (Delta < 0) Delta = 0;
                                    let window_size = TWO_HIT_WINDOW;
                                    let delta_calc = window_size as isize - word_length as isize;
                                    let delta_max = if delta_calc < 0 {
                                        0
                                    } else {
                                        scan_range.min(delta_calc as usize) as isize
                                    };

                                    // NCBI reference: na_ungapped.c:689-690
                                    // Int4 s_a = s_off_pos + word_length - window_size;
                                    // Int4 s_b = s_end_pos - 2 * word_length;
                                    // CRITICAL: NCBI uses signed arithmetic (Int4), so s_a and s_b can be negative
                                    // LOSAT must use signed arithmetic to match NCBI behavior
                                    let s_a = s_off_pos as isize + word_length as isize - window_size as isize;
                                    let s_b = s_end_pos as isize - 2 * word_length as isize;

                                    // NCBI reference: na_ungapped.c:688
                                    // Int4 orig_diag = real_diag + diag_table->diag_array_length;
                                    // In LOSAT, we use diag directly (not masked)
                                    let orig_diag = diag;

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
                                        let off_diag = orig_diag + delta;
                                        let (off_s_end, off_s_l) = if use_array_indexing {
                                            let off_diag_idx = (off_diag + diag_offset) as usize;
                                            if off_diag_idx < diag_array_size {
                                                let off_entry = &hit_level_array[off_diag_idx];
                                                (off_entry.last_hit, hit_len_array[off_diag_idx])
                                            } else {
                                                (0, 0)
                                            }
                                        } else {
                                            let off_diag_key = pack_diag_key(q_idx, off_diag);
                                            let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                            let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                            (off_entry.last_hit, off_len)
                                        };

                                        // NCBI: off_s_end - delta >= s_a (signed comparison)
                                        // Convert to signed for comparison to match NCBI behavior
                                        if off_s_l > 0
                                            && (off_s_end as isize - delta) >= s_a
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
                                        let off_diag = orig_diag - delta;
                                        let (off_s_end, off_s_l) = if use_array_indexing {
                                            let off_diag_idx = (off_diag + diag_offset) as usize;
                                            if off_diag_idx < diag_array_size {
                                                let off_entry = &hit_level_array[off_diag_idx];
                                                (off_entry.last_hit, hit_len_array[off_diag_idx])
                                            } else {
                                                (0, 0)
                                            }
                                        } else {
                                            let off_diag_key = pack_diag_key(q_idx, off_diag);
                                            let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                            let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                            (off_entry.last_hit, off_len)
                                        };

                                        // NCBI: off_s_end >= s_a (signed comparison)
                                        // Convert to signed for comparison to match NCBI behavior
                                        if off_s_l > 0
                                            && (off_s_end as isize) >= s_a
                                            && (off_s_end as isize - off_s_l as isize + delta) <= s_b {
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
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:718-725
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
                            dbg_ungapped_ext_calls += 1;
                            let (qs, qe, ss, ungapped_se, ungapped_score) = extend_hit_ungapped(
                                q_seq,
                                search_seq,
                                q_off,
                                s_off,
                                reward,
                                penalty,
                                None, // Use default X-drop
                            );

                            // NCBI reference: na_ungapped.c:757-758
                            // s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                            // This is the END of the UNGAPPED extension, used for last_hit update
                            let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

                            // Skip if ungapped score is too low
                            // NCBI reference: na_ungapped.c:752
                            // if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
                            // Use dynamically calculated cutoff_score instead of fixed threshold
                            // off_found is set by off-diagonal search (Step 3)
                            // NCBI: if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
                            if !(off_found || ungapped_score >= cutoff_score) {
                                dbg_ungapped_low += 1;
                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos_usize, kmer_start, ungapped_score, min_ungapped_score);
                                }
                                // NCBI reference: na_ungapped.c:768-771
                                // Update hit_level_array and hit_len_array even if extension is skipped
                                // hit_level_array[real_diag].last_hit = s_end_pos;
                                // hit_level_array[real_diag].flag = hit_ready;
                                // diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;
                                if use_array_indexing && diag_idx < diag_array_size {
                                    hit_level_array[diag_idx].last_hit = s_end_pos;
                                    hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                                    hit_len_array[diag_idx] = if hit_ready { 0 } else { s_end_pos - s_off_pos };
                                } else if !use_array_indexing {
                                    let diag_key = pack_diag_key(q_idx, diag);
                                    hit_level_hash.insert(diag_key, DiagStruct {
                                        last_hit: s_end_pos,
                                        flag: if hit_ready { 1 } else { 0 },
                                    });
                                    hit_len_hash.insert(diag_key, if hit_ready { 0 } else { s_end_pos - s_off_pos });
                                }
                                continue;
                            }

                            dbg_gapped_attempted += 1;

                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> COLLECT UNGAPPED (score={}, len={})", q_pos_usize, kmer_start, ungapped_score, qe - qs);
                            }

                            // NCBI architecture: Collect ungapped hit for batch processing
                            // NCBI reference: blast_gapalign.c:3824 - init_hsp_array is sorted by score DESCENDING
                            // Gapped extension will be done later in score order with containment check
                            ungapped_hits.push(UngappedHit {
                                q_idx,
                                qs,
                                qe,
                                ss,
                                se: ungapped_se,
                                score: ungapped_score,
                                is_minus_strand,
                                q_seq_idx: q_idx as usize,
                            });

                            // NCBI reference: na_ungapped.c:768-771
                            // Update hit_level_array after collecting ungapped hit
                            hit_ready = true;
                            if use_array_indexing && diag_idx < diag_array_size {
                                hit_level_array[diag_idx].last_hit = ungapped_s_end_pos;
                                hit_level_array[diag_idx].flag = 1;
                                hit_len_array[diag_idx] = 0;
                            } else if !use_array_indexing {
                                let diag_key = pack_diag_key(q_idx, diag);
                                hit_level_hash.insert(diag_key, DiagStruct {
                                    last_hit: ungapped_s_end_pos,
                                    flag: 1,
                                });
                                hit_len_hash.insert(diag_key, 0);
                            }
                            // Gapped extension deferred to batch processing phase
                        } // end of for matches_slice in two-stage lookup
                    } // end of s_pos loop for two-stage lookup

                    // DEBUG: Log stats for this subject
                    let elapsed = dbg_start_time.elapsed();
                    eprintln!("[PERF] Subject scan took {:?}: left_ext={}, right_ext={}, ungapped={}, seeds={}, valid_pos={}",
                        elapsed, dbg_left_ext_iters, dbg_right_ext_iters, dbg_ungapped_ext_calls,
                        dbg_seeds_found, dbg_total_s_positions);
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

                            // Use pre-computed cutoff score (computed once per query-subject pair)
                            let cutoff_score = cutoff_scores[q_idx as usize];
                            let q_record = &queries[q_idx as usize];

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

                        // NCBI BLAST does NOT check mask_array/mask_hash at seed level
                        // NCBI reference: na_ungapped.c:671-672
                        // NCBI only checks: if (s_off_pos < last_hit) return 0;
                        // This check happens below using hit_level_array[real_diag].last_hit
                        // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                        // REMOVED to match NCBI BLAST behavior

                        // NCBI reference: na_ungapped.c:656-672
                        // Two-hit filter: check if there's a previous hit within window_size
                        // NCBI: Boolean two_hits = (window_size > 0);
                        // NCBI: last_hit = hit_level_array[real_diag].last_hit;
                        // NCBI: hit_saved = hit_level_array[real_diag].flag;
                        // NCBI: if (s_off_pos < last_hit) return 0;  // hit within explored area
                        // NCBI: if (two_hits && (hit_saved || s_end_pos > last_hit + window_size)) { ... }
                        let diag_idx = if use_array_indexing {
                            (diag + diag_offset) as usize
                        } else {
                            0 // Not used for hash indexing
                        };

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

                        // NCBI: s_off_pos = s_off + diag_table->offset;
                        // In LOSAT, we use kmer_start directly (0-based), but need to add offset for comparison
                        let s_off_pos = kmer_start + diag_offset as usize;

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
                        let two_hits = TWO_HIT_WINDOW > 0;
                        let q_off = q_pos as usize;
                        let s_off = kmer_start;
                        let mut s_end = kmer_start + safe_k;
                        let mut s_end_pos = s_end + diag_offset as usize;
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
                            // For non-two-stage lookup, word_length == lut_word_length, so type_of_word returns (1, 0)
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
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
                                // Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
                                // if (Delta < 0) Delta = 0;
                                let window_size = TWO_HIT_WINDOW;
                                let delta_calc = window_size as isize - safe_k as isize;
                                let delta_max = if delta_calc < 0 {
                                    0
                                } else {
                                    scan_range.min(delta_calc as usize) as isize
                                };

                                // NCBI reference: na_ungapped.c:855-856
                                // Int4 s_a = s_off_pos + word_length - window_size;
                                // Int4 s_b = s_end_pos - 2 * word_length;
                                // CRITICAL: NCBI uses signed arithmetic (Int4), so s_a and s_b can be negative
                                // LOSAT must use signed arithmetic to match NCBI behavior
                                let s_a = s_off_pos as isize + safe_k as isize - window_size as isize;
                                let s_b = s_end_pos as isize - 2 * safe_k as isize;

                                let orig_diag = diag;

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
                                    let off_diag = orig_diag + delta;
                                    let (off_s_end, off_s_l) = if use_array_indexing {
                                        let off_diag_idx = (off_diag + diag_offset) as usize;
                                        if off_diag_idx < diag_array_size {
                                            let off_entry = &hit_level_array[off_diag_idx];
                                            (off_entry.last_hit, hit_len_array[off_diag_idx])
                                        } else {
                                            (0, 0)
                                        }
                                    } else {
                                        let off_diag_key = pack_diag_key(q_idx, off_diag);
                                        let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                        let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                        (off_entry.last_hit, off_len)
                                    };

                                    // NCBI: off_s_end - delta >= s_a (signed comparison)
                                    // Convert to signed for comparison to match NCBI behavior
                                    if off_s_l > 0
                                        && (off_s_end as isize - delta) >= s_a
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
                                    let off_diag = orig_diag - delta;
                                    let (off_s_end, off_s_l) = if use_array_indexing {
                                        let off_diag_idx = (off_diag + diag_offset) as usize;
                                        if off_diag_idx < diag_array_size {
                                            let off_entry = &hit_level_array[off_diag_idx];
                                            (off_entry.last_hit, hit_len_array[off_diag_idx])
                                        } else {
                                            (0, 0)
                                        }
                                    } else {
                                        let off_diag_key = pack_diag_key(q_idx, off_diag);
                                        let off_entry = hit_level_hash.get(&off_diag_key).copied().unwrap_or_default();
                                        let off_len = hit_len_hash.get(&off_diag_key).copied().unwrap_or(0);
                                        (off_entry.last_hit, off_len)
                                    };

                                    // NCBI: off_s_end >= s_a (signed comparison)
                                    // Convert to signed for comparison to match NCBI behavior
                                    if off_s_l > 0
                                        && (off_s_end as isize) >= s_a
                                        && (off_s_end as isize - off_s_l as isize + delta) <= s_b {
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
                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:718-725
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
                        let (qs, qe, ss, ungapped_se, ungapped_score) = extend_hit_ungapped(
                            q_record.seq(),
                            search_seq,
                            q_off,
                            s_off,
                            reward,
                            penalty,
                            None, // Use default X-drop
                        );

                        // NCBI reference: na_ungapped.c:757-758
                        // s_end_pos = ungapped_data->length + ungapped_data->s_start + diag_table->offset;
                        // This is the END of the UNGAPPED extension, used for last_hit update
                        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

                        // Skip if ungapped score is too low
                        // NCBI reference: na_ungapped.c:752
                        // if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
                        // Use dynamically calculated cutoff_score instead of fixed threshold
                        // off_found is set by off-diagonal search
                        // NCBI: if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
                        if !(off_found || ungapped_score >= cutoff_score) {
                            // NCBI reference: na_ungapped.c:768-771
                            // Update hit_level_array and hit_len_array even if extension is skipped
                            if use_array_indexing && diag_idx < diag_array_size {
                                hit_level_array[diag_idx].last_hit = s_end_pos;
                                hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                                hit_len_array[diag_idx] = if hit_ready { 0 } else { s_end_pos - s_off_pos };
                            } else if !use_array_indexing {
                                let diag_key = pack_diag_key(q_idx, diag);
                                hit_level_hash.insert(diag_key, DiagStruct {
                                    last_hit: s_end_pos,
                                    flag: if hit_ready { 1 } else { 0 },
                                });
                                hit_len_hash.insert(diag_key, if hit_ready { 0 } else { s_end_pos - s_off_pos });
                            }
                            dbg_ungapped_low += 1;
                            if in_window && debug_mode {
                                eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos, kmer_start, ungapped_score, min_ungapped_score);
                            }
                            continue;
                        }

                        dbg_gapped_attempted += 1;

                        if in_window && debug_mode {
                            eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> COLLECT UNGAPPED (score={}, len={})", q_pos, kmer_start, ungapped_score, qe - qs);
                        }

                        // NCBI architecture: Collect ungapped hit for batch processing
                        // NCBI reference: blast_gapalign.c:3824 - init_hsp_array is sorted by score DESCENDING
                        ungapped_hits.push(UngappedHit {
                            q_idx,
                            qs,
                            qe,
                            ss,
                            se: ungapped_se,
                            score: ungapped_score,
                            is_minus_strand,
                            q_seq_idx: q_idx as usize,
                        });

                        // NCBI reference: na_ungapped.c:768-771
                        // Update hit_level_array after collecting ungapped hit
                        hit_ready = true;
                        if use_array_indexing && diag_idx < diag_array_size {
                            hit_level_array[diag_idx].last_hit = ungapped_s_end_pos;
                            hit_level_array[diag_idx].flag = 1;
                            hit_len_array[diag_idx] = 0;
                        } else if !use_array_indexing {
                            let diag_key = pack_diag_key(q_idx, diag);
                            hit_level_hash.insert(diag_key, DiagStruct {
                                last_hit: ungapped_s_end_pos,
                                flag: 1,
                            });
                            hit_len_hash.insert(diag_key, 0);
                        }
                    }
                    } // end of s_pos loop for original lookup
                } // end of if/else for two-stage vs original lookup
            } // end of strand loop

            // ============================================================
            // NCBI ARCHITECTURE: Batch process ungapped hits in score order
            // ============================================================
            // NCBI reference: blast_gapalign.c:3824 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))
            // NCBI processes init_hsp_array sorted by score DESCENDING
            // High-score HSPs are processed first, their GAPPED HSPs added to interval tree
            // Lower-score HSPs are then checked for containment against these GAPPED HSPs

            // Sort by score DESCENDING (highest first)
            ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

            // Debug counters for containment analysis
            let mut dbg_containment_skipped = 0usize;
            let total_ungapped = ungapped_hits.len();
            let gapped_start_time = std::time::Instant::now();
            let mut dbg_gapped_calls = 0usize;

            // Process each ungapped hit in score order
            for (idx, uh) in ungapped_hits.iter().enumerate() {
                // NCBI reference: blast_gapalign.c:3886-4090 (no per-HSP logging)
                if verbose && idx % 100 == 0 {
                    eprintln!("[INFO] Gapped extension: {}/{}", idx, total_ungapped);
                }
                // Get query sequence for this hit
                let q_record = &queries[uh.q_seq_idx];
                let q_seq = q_record.seq();

                // Get subject sequence (forward or reverse complement)
                let search_seq: &[u8] = if uh.is_minus_strand {
                    &s_seq_rc
                } else {
                    s_seq
                };

                // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:866-869 (query blastna, subject ncbi2na)
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:503-507 (greedy traceback uses uncompressed subject)
                let (q_seq_greedy, s_seq_greedy_score, s_seq_greedy_trace) = if use_dp {
                    (q_seq, search_seq, search_seq)
                } else {
                    let encoded_queries = &encoded_queries_blastna;
                    let q_seq_blastna = encoded_queries[uh.q_seq_idx].as_slice();
                    let s_seq_score = if uh.is_minus_strand {
                        s_seq_rc_ncbi2na.as_ref().unwrap().as_slice()
                    } else {
                        s_seq_ncbi2na.as_ref().unwrap().as_slice()
                    };
                    let s_seq_trace = if uh.is_minus_strand {
                        s_seq_rc_blastna.as_slice()
                    } else {
                        s_seq_blastna.as_slice()
                    };
                    (q_seq_blastna, s_seq_score, s_seq_trace)
                };

                // NCBI reference: blast_gapalign.c:3908-3915
                // Create tmp_hsp from ungapped data for containment checking
                // NCBI uses 0-based coordinates internally for the tree
                let subject_frame_sign = if uh.is_minus_strand { -1i32 } else { 1i32 };
                let ungapped_tree_hsp = TreeHsp {
                    query_offset: uh.qs as i32,
                    query_end: uh.qe as i32,
                    subject_offset: uh.ss as i32,
                    subject_end: uh.se as i32,
                    score: uh.score,
                    query_context_offset: 0, // Single query context in LOSAT
                    subject_frame_sign,
                };

                // NCBI reference: blast_gapalign.c:3918 BlastIntervalTreeContainsHSP
                // Check if UNGAPPED HSP is contained in existing GAPPED HSPs
                let is_contained = interval_tree.contains_hsp(&ungapped_tree_hsp, 0, min_diag_separation);

                if is_contained {
                    // NCBI: Skip gapped extension if ungapped HSP is contained
                    dbg_containment_skipped += 1;
                    continue;
                }

                // Select gapped-start seed within the ungapped HSP.
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3248-3389
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4031
                let cutoff_score = cutoff_scores[uh.q_idx as usize];
                let (
                    prelim_qs,
                    prelim_qe,
                    prelim_ss,
                    prelim_se,
                    prelim_score,
                    final_seed_qs,
                    final_seed_ss,
                ) = if use_dp {
                    // DP seed selection (blastn)
                    let (seed_qs, seed_ss) = match blast_get_offsets_for_gapped_alignment(
                        q_seq,
                        search_seq,
                        uh.qs,
                        uh.qe,
                        uh.ss,
                        uh.se,
                        reward,
                        penalty,
                    ) {
                        Some((q_start, s_start)) => blast_get_start_for_gapped_alignment_nucl(
                            q_seq,
                            search_seq,
                            uh.qs,
                            uh.qe,
                            uh.ss,
                            uh.se,
                            q_start,
                            s_start,
                        ),
                        None => {
                            continue;
                        }
                    };

                    // Preliminary DP gapped extension (score-only)
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:176-183
                    let (p_qs, p_qe, p_ss, p_se, p_score, _, _, _, _, _) = extend_gapped_heuristic_with_scratch(
                        q_seq,
                        search_seq,
                        seed_qs,
                        seed_ss,
                        1,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop_gapped,
                        gap_scratch,
                        true,
                    );
                    (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss)
                } else {
                    // Greedy seed selection (megablast)
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4017
                    let ungapped_len = uh.qe.saturating_sub(uh.qs);
                    let seed_qs = uh.qs + ungapped_len / 2;
                    let seed_ss = uh.ss + ungapped_len / 2;

                    let prelim = match greedy_gapped_alignment_score_only(
                        q_seq_greedy,
                        s_seq_greedy_score,
                        seed_qs,
                        seed_ss,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop_gapped,
                    ) {
                        Some(value) => value,
                        None => {
                            continue;
                        }
                    };

                    let (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss) = prelim;
                    (p_qs, p_qe, p_ss, p_se, p_score, seed_qs, seed_ss)
                };

                dbg_gapped_calls += 1;
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4091
                if prelim_score < cutoff_score {
                    continue;
                }

                // Add preliminary GAPPED HSP to interval tree for containment checks
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4089
                let gapped_tree_hsp = TreeHsp {
                    query_offset: prelim_qs as i32,
                    query_end: prelim_qe as i32,
                    subject_offset: prelim_ss as i32,
                    subject_end: prelim_se as i32,
                    score: prelim_score,
                    query_context_offset: 0,
                    subject_frame_sign: if uh.is_minus_strand { -1 } else { 1 },
                };
                interval_tree.add_hsp(gapped_tree_hsp, 0, IndexMethod::QueryAndSubject);

                // Final gapped extension with traceback
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:503-512
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:1422-1431
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
                    edit_ops,
                ) = if use_dp {
                    extend_gapped_heuristic_with_traceback_with_scratch(
                        q_seq,
                        search_seq,
                        final_seed_qs,
                        final_seed_ss,
                        1,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop_final,
                        gap_scratch,
                    )
                } else {
                    match greedy_gapped_alignment_with_traceback(
                        q_seq_greedy,
                        s_seq_greedy_trace,
                        final_seed_qs,
                        final_seed_ss,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        x_drop_final,
                    ) {
                        Some(value) => value,
                        None => {
                            continue;
                        }
                    }
                };

                // Calculate statistics
                let aln_len = matches + mismatches + gap_letters;

                // Calculate e-value
                let (bit_score, eval) = calculate_evalue(
                    score,
                    q_seq.len(),
                    db_len_total,
                    db_num_seqs,
                    &params_for_closure,
                );

                if eval > evalue_threshold {
                    continue;
                }

                let identity = if aln_len > 0 {
                    ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                } else {
                    0.0
                };

                let q_id = query_ids[uh.q_idx as usize].clone();

                // Convert coordinates for minus strand
                let (hit_s_start, hit_s_end) = if uh.is_minus_strand {
                    let orig_s_start = s_len - final_ss;
                    let orig_s_end = s_len - final_se + 1;
                    (orig_s_start, orig_s_end)
                } else {
                    (final_ss + 1, final_se)
                };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2567 (gap_align->edit_script stored in HSP)
                let gap_info = if edit_ops.is_empty() {
                    None
                } else {
                    Some(edit_ops)
                };

                // NCBI reference: blast_traceback.c:679-692
                // Store internal 0-based coordinates for Phase 2 interval tree processing
                let internal = InternalHitData {
                    q_offset_0: final_qs,
                    q_end_0: final_qe,
                    s_offset_0: final_ss,  // Always canonical (0-based, search_seq coords)
                    s_end_0: final_se,     // Always canonical (0-based, search_seq coords)
                    is_minus_strand: uh.is_minus_strand,
                };

                // Add hit with internal coordinates to results
                hits_with_internal.push((Hit {
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
                    q_idx: uh.q_idx,  // Store query index for re-evaluation
                    s_idx: 0,
                    raw_score: score,
                    gap_info,
                }, internal));

            }

            // DEBUG: Log gapped extension stats
            let gapped_elapsed = gapped_start_time.elapsed();
            if verbose {
                eprintln!(
                    "[INFO] Gapped extension took {:?}: calls={}, skipped={}, total_ungapped={}",
                    gapped_elapsed,
                    dbg_gapped_calls,
                    dbg_containment_skipped,
                    total_ungapped
                );
            }

            // Print debug summary for this subject (before post-processing)
            let hits_before_phase2 = hits_with_internal.len();
            if debug_mode || blastn_debug {
                // Get the cutoff_score for the first query (typically only one query in blastn)
                let cutoff_score_0 = cutoff_scores.first().copied().unwrap_or(0);
                eprintln!(
                    "[DEBUG] Subject {}: positions={}, seeds_found={}, ungapped_low={}, two_hit_failed={}, gapped_attempted={}, hits={}, cutoff_score={}",
                    s_id, dbg_total_s_positions, dbg_seeds_found, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, hits_before_phase2, cutoff_score_0
                );
                // Show containment statistics
                eprintln!(
                    "[DEBUG] Interval tree: ungapped_collected={}, containment_skipped={}, gapped_added={}, tree_nodes={}",
                    total_ungapped, dbg_containment_skipped, interval_tree.node_count(), hits_before_phase2
                );
            }

            // Suppress unused variable warnings when not in debug mode
            let _ = (dbg_total_s_positions, dbg_ambiguous_skipped, dbg_no_lookup_match, dbg_seeds_found, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds, total_ungapped, dbg_containment_skipped);

            // =================================================================
            // NCBI POST-GAPPED PROCESSING
            // Reference: blast_traceback.c:633-692
            // =================================================================

            // Step 1: Extract hits and internal coordinates
            let (local_hits, internals): (Vec<Hit>, Vec<InternalHitData>) =
                hits_with_internal.into_iter().unzip();

            // Step 2: Endpoint purging pass 1 (trim, purge=false)
            // NCBI reference: blast_traceback.c:637-638
            // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
            let hits_step1 = local_hits.len();
            let (mut local_hits, extra_start) = purge_hsps_with_common_endpoints_ex(local_hits, false);
            let hits_step2 = local_hits.len();

            // Step 3: Re-evaluate trimmed HSPs
            // NCBI reference: blast_traceback.c:647-665
            // The remaining part of the hsp may be extended further

            // Create ReevalParams for E-value/bit_score recalculation
            // NCBI reference: blast_traceback.c:234-250 Blast_HSPListGetEvalues, Blast_HSPListGetBitScores
            let reeval_params = ReevalParams {
                lambda: params_for_closure.lambda,
                k: params_for_closure.k,
                eff_searchsp: (db_len_total as i64) * (queries[0].seq().len() as i64), // Approximate
                db_len: db_len_total,
                db_num_seqs,
            };

            for hit in local_hits.iter_mut().skip(extra_start) {
                // Get sequences for re-evaluation
                let q_seq = queries[hit.q_idx as usize].seq();
                let cutoff = cutoff_scores.get(hit.q_idx as usize).copied().unwrap_or(0);

                // NCBI reference: blast_traceback.c:653-665 (reevaluate with blastna sequences)
                let q_seq_blastna = encoded_queries_blastna[hit.q_idx as usize].as_slice();
                let s_seq_eval = if hit.s_start > hit.s_end {
                    s_seq_rc_blastna.as_slice()
                } else {
                    s_seq_blastna.as_slice()
                };
                let delete = reevaluate_hsp_with_ambiguities_gapped_ex(
                    hit,
                    q_seq_blastna,
                    s_seq_eval,
                    reward,
                    penalty,
                    gap_open,
                    gap_extend,
                    cutoff,
                    Some(&reeval_params),
                );
                if delete {
                    hit.raw_score = i32::MIN; // Mark for removal
                }
            }
            local_hits.retain(|h| h.raw_score != i32::MIN);
            let hits_step3 = local_hits.len();

            // Step 4: Endpoint purging pass 2 (delete, purge=true) - BLASTN only
            // NCBI reference: blast_traceback.c:667-668
            // if(program_number == eBlastTypeBlastn) {
            //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
            // }
            let mut local_hits = purge_hsps_with_common_endpoints(local_hits);
            let hits_step4 = local_hits.len();

            // Step 5: Re-sort by gapped score
            // NCBI reference: blast_traceback.c:672
            // Blast_HSPListSortByScore(hsp_list);
            local_hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));

            // Step 6: Phase 2 - Tree reset and second containment pass
            // NCBI reference: blast_traceback.c:678
            // Blast_IntervalTreeReset(tree);
            interval_tree.reset();

            // NCBI reference: blast_traceback.c:679-692
            // Remove any HSPs that are contained within other HSPs.
            // Since the list is sorted by score already, any HSP
            // contained by a previous HSP is guaranteed to have a
            // lower score, and may be purged.
            let mut final_hits = Vec::with_capacity(local_hits.len());
            for hit in local_hits {
                // Reconstruct TreeHsp from hit (convert 1-based back to 0-based)
                // NCBI uses canonical coordinates: subject.offset < subject.end always
                let tree_hsp = TreeHsp {
                    query_offset: (hit.q_start.saturating_sub(1)) as i32,
                    query_end: hit.q_end as i32,
                    subject_offset: (hit.s_start.min(hit.s_end).saturating_sub(1)) as i32,
                    subject_end: hit.s_start.max(hit.s_end) as i32,
                    score: hit.raw_score,
                    query_context_offset: 0,
                    subject_frame_sign: if hit.s_start > hit.s_end { -1 } else { 1 },
                };

                // NCBI reference: blast_traceback.c:682-688
                if !interval_tree.contains_hsp(&tree_hsp, 0, min_diag_separation) {
                    interval_tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);
                    final_hits.push(hit);
                }
                // else: HSP is contained within another, skip (implicit delete)
            }

            // Phase 2 debug output
            if debug_mode || blastn_debug {
                let phase2_tree_removed = hits_step4 - final_hits.len();
                eprintln!(
                    "[DEBUG] Phase 2 breakdown: step1(extract)={}, step2(purge1)={} (-{}), step3(reeval)={} (-{}), step4(purge2)={} (-{}), step6(tree)={} (-{})",
                    hits_step1,
                    hits_step2, hits_step1 - hits_step2,
                    hits_step3, hits_step2 - hits_step3,
                    hits_step4, hits_step3 - hits_step4,
                    final_hits.len(), phase2_tree_removed
                );
            }

            // Suppress unused warning
            let _ = internals;

            if !final_hits.is_empty() {
                // NCBI reference: blast_traceback.c:633-692 (post-gapped processing is per-subject)
                tx.send(Some(final_hits)).unwrap();
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
