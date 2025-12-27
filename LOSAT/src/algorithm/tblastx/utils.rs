//! Main execution logic for TBLASTX
//!
//! This module contains the main `run` function that orchestrates the TBLASTX search.

use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::atomic::Ordering as AtomicOrdering;
use std::sync::mpsc::channel;

use crate::common::{write_output, Hit};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::{lookup_protein_params, search_space::SearchSpace};
use crate::utils::dust::{DustMasker, MaskedInterval};
use crate::utils::genetic_code::GeneticCode;

use super::args::TblastxArgs;
use super::constants::{
    CUTOFF_E_TBLASTX, GAP_EXTEND, GAP_OPEN, TWO_HIT_WINDOW, X_DROP_GAPPED_FINAL,
    X_DROP_GAPPED_PRELIM,
};
use super::chaining::{ExtendedHit, SequenceData, SequenceKey};
use super::diagnostics::{diagnostics_enabled, DiagnosticCounters, print_summary as print_diagnostics_summary};
use super::extension::{
    convert_coords, extend_gapped_protein, extend_hit_two_hit,
    get_score,
};
use crate::stats::karlin::{bit_score as calc_bit_score, evalue as calc_evalue, raw_score_from_evalue};
use super::lookup::build_direct_lookup;
use super::translation::{generate_frames, QueryFrame};

pub fn run(args: TblastxArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);

    // Initialize diagnostic counters if enabled
    let diag_enabled = diagnostics_enabled();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    eprintln!("Reading queries...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();

    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| {
            r.id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string()
        })
        .collect();

    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        queries_raw
            .iter()
            .map(|r| masker.mask_sequence(r.seq()))
            .collect()
    } else {
        queries_raw.iter().map(|_| Vec::new()).collect()
    };

    let query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();

    eprintln!("Building lookup table...");
    let lookup = build_direct_lookup(&query_frames, &query_masks);

    eprintln!("Reading subject file...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> =
        subject_reader.records().filter_map(|r| r.ok()).collect();

    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    let db_len_total_bp: usize = subjects_raw.iter().map(|r| r.seq().len()).sum();
    let db_len_aa_total = db_len_total_bp / 3;

    let scoring_spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: GAP_OPEN.abs(),
        gap_extend: GAP_EXTEND.abs(),
    };
    let params = lookup_protein_params(&scoring_spec);

    eprintln!(
        "Searching {} Queries vs {} Subjects...",
        queries_raw.len(),
        subjects_raw.len()
    );

    let bar = ProgressBar::new(subjects_raw.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap(),
    );

    // Channel now sends ExtendedHit with sequence data for HSP chaining
    type HitData = (Vec<ExtendedHit>, Vec<(SequenceKey, SequenceData)>);
    let (tx, rx) = channel::<HitData>();

    let out_path = args.out.clone();
    let params_clone = params.clone();
    let evalue_threshold = args.evalue;
    let diagnostics_clone = diagnostics.clone();
    let diag_enabled_clone = diag_enabled;
    let writer_handle = std::thread::spawn(move || -> Result<()> {
        let mut all_hits: Vec<ExtendedHit> = Vec::new();

        while let Ok((hits, _seq_data)) = rx.recv() {
            all_hits.extend(hits);
            // Note: seq_data is not used since we skip HSP chaining for TBLASTX
            // (NCBI BLAST does not call Blast_LinkHsps for TBLASTX)
        }

        // NCBI BLAST does not call Blast_LinkHsps for TBLASTX
        // Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1129
        // ASSERT(program_number != eBlastTypeTblastx);
        // For TBLASTX, individual HSPs are output directly without chaining
        // Only apply e-value threshold filtering
        let diag_ref = if diag_enabled_clone { Some(diagnostics_clone.as_ref()) } else { None };
        
        // Record HSPs before filtering
        let hsps_before = all_hits.len();
        if let Some(diag) = diag_ref {
            diag.base.hsps_before_chain.store(hsps_before, AtomicOrdering::Relaxed);
            diag.base.hsps_after_chain.store(hsps_before, AtomicOrdering::Relaxed);
        }

        // Sort by bit score (highest first) for consistent output order
        all_hits.sort_by(|a, b| {
            b.hit.bit_score
                .partial_cmp(&a.hit.bit_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // NCBI BLAST does not apply HSP culling for TBLASTX by default
        // (see blast_options.c: "Gapped search is not allowed for tblastx")
        // Skip domination filter to match NCBI BLAST behavior
        if let Some(diag) = diag_ref {
            diag.base.hsps_culled_dominated.store(0, AtomicOrdering::Relaxed);
        }

        // Convert to Hit and filter by e-value threshold
        let mut output_from_ungapped = 0usize;
        let mut output_from_gapped = 0usize;
        let filtered_hits: Vec<Hit> = all_hits
            .into_iter()
            .filter(|ext_hit| ext_hit.hit.e_value <= evalue_threshold)
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
        if let Some(diag) = diag_ref {
            diag.output_from_ungapped.store(output_from_ungapped, AtomicOrdering::Relaxed);
            diag.output_from_gapped.store(output_from_gapped, AtomicOrdering::Relaxed);
            diag.base.hsps_after_overlap_filter.store(filtered_hits.len(), AtomicOrdering::Relaxed);
        }

        write_output(&filtered_hits, out_path.as_ref())?;
        Ok(())
    });

    let diagnostics_for_parallel = diagnostics.clone();
    subjects_raw
        .par_iter()
        .enumerate()
        .for_each_with((tx, diagnostics_for_parallel), |(tx, diag), (_s_idx, s_record)| {
            let s_frames = generate_frames(s_record.seq(), &db_code);
            let s_id = s_record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            let s_len = s_record.seq().len();

            // Parallelize over subject frames for intra-subject parallelism
            // Each frame is independent since mask_key includes s_f_idx
            let diag_inner = diag.clone();
            let frame_results: Vec<(Vec<ExtendedHit>, Vec<(SequenceKey, SequenceData)>)> =
                s_frames
                    .par_iter()
                    .enumerate()
                    .map(|(s_f_idx, s_frame)| {
                        let mut local_hits: Vec<ExtendedHit> = Vec::new();
                        let mut local_sequences: Vec<(SequenceKey, SequenceData)> = Vec::new();
                        let mut seen_pairs: std::collections::HashSet<SequenceKey> =
                            std::collections::HashSet::new();

                        // Mask to track already-extended regions on each diagonal (per frame)
                        // NCBI BLAST style: use diag_coord = (query_offset - subject_offset) & diag_mask
                        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:351
                        //            ncbi-blast/c++/src/algo/blast/core/blast_extend.c:42-61
                        let mut mask: FxHashMap<(u32, u8, isize), usize> = FxHashMap::default();
                        // Two-hit tracking: stores the last seed position on each diagonal for each query/frame combination
                        // NCBI BLAST style: use diag_coord instead of raw diagonal
                        let mut last_seed: FxHashMap<(u32, u8, isize), usize> =
                            FxHashMap::default();
                        let s_aa = &s_frame.aa_seq;
                        let s_aa_len = s_aa.len();
                        
                        // Calculate diag_mask for each query frame (NCBI BLAST style)
                        // diag_array_length = smallest power of 2 >= (query_length + window_size)
                        // diag_mask = diag_array_length - 1
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
                        let mut diag_mask_cache: FxHashMap<(u32, u8), isize> = FxHashMap::default();
                        for (q_idx, q_frames) in query_frames.iter().enumerate() {
                            for (q_f_idx, q_frame) in q_frames.iter().enumerate() {
                                let q_aa_len = q_frame.aa_seq.len();
                                // Calculate diag_array_length as smallest power of 2 >= (q_aa_len + TWO_HIT_WINDOW)
                                let mut diag_array_length = 1usize;
                                while diag_array_length < (q_aa_len + TWO_HIT_WINDOW) {
                                    diag_array_length <<= 1;
                                }
                                let diag_mask = (diag_array_length - 1) as isize;
                                diag_mask_cache.insert((q_idx as u32, q_f_idx as u8), diag_mask);
                            }
                        }
                        
                        // Cache effective search space per query-subject-frame combination
                        // NCBI BLAST calculates this once per context (query frame + subject frame)
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846
                        // IMPORTANT: context_key must include subject frame index to match NCBI BLAST behavior
                        let mut search_space_cache: FxHashMap<(u32, u8, u8), SearchSpace> = FxHashMap::default();
                        // Cache cutoff_score per context (calculated from CUTOFF_E_TBLASTX = 1e-300)
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:349-360
                        let mut cutoff_score_cache: FxHashMap<(u32, u8, u8), i32> = FxHashMap::default();
                        if s_aa.len() < 3 {
                            return (local_hits, local_sequences);
                        }

                        for s_pos in 0..=(s_aa.len() - 3) {
                            let c1 = unsafe { *s_aa.get_unchecked(s_pos) } as usize;
                            let c2 = unsafe { *s_aa.get_unchecked(s_pos + 1) } as usize;
                            let c3 = unsafe { *s_aa.get_unchecked(s_pos + 2) } as usize;

                            // シードに終止コドンが含まれる場合はスキップ
                            // NCBI BLAST uses 25 amino acids: ARNDCQEGHILKMFPSTWYVBJZX* (indices 0-24)
                            if c1 >= 25 || c2 >= 25 || c3 >= 25 {
                                continue;
                            }
                            // 25^3 = 15,625 possible k-mers
                            let kmer = c1 * 625 + c2 * 25 + c3;

                            let matches = unsafe { lookup.get_unchecked(kmer) };

                            // Count k-mer matches for diagnostics
                            if diag_enabled && !matches.is_empty() {
                                diag_inner.base.kmer_matches.fetch_add(matches.len(), AtomicOrdering::Relaxed);
                            }

                            for &(q_idx, q_f_idx, q_pos) in matches {
                                // NCBI BLAST style: diag_coord = (query_offset - subject_offset) & diag_mask
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:351
                                // Note: NCBI BLAST uses query_offset - subject_offset, but we use q_pos - s_pos
                                // which is the same diagonal but with opposite sign convention
                                let diag = q_pos as isize - s_pos as isize;
                                // Get diag_mask for this query frame
                                let diag_mask = *diag_mask_cache.get(&(q_idx, q_f_idx))
                                    .unwrap_or(&((1isize << 20) - 1)); // Fallback: 2^20 - 1
                                // Calculate diag_coord using NCBI BLAST style masking
                                let diag_coord = diag & diag_mask;
                                // Use diag_coord instead of raw diagonal for mask_key (NCBI BLAST compatible)
                                let mask_key = (q_idx, q_f_idx, diag_coord);

                                // Check if this region was already extended
                                if let Some(&last_end) = mask.get(&mask_key) {
                                    if s_pos <= last_end {
                                        continue;
                                    }
                                }

                                let q_frame = &query_frames[q_idx as usize][q_f_idx as usize];
                                let q_aa = &q_frame.aa_seq;
                                let q_aa_len = q_aa.len();

                                let seed_score = get_score(c1 as u8, q_aa[q_pos as usize])
                                    + get_score(c2 as u8, q_aa[q_pos as usize + 1])
                                    + get_score(c3 as u8, q_aa[q_pos as usize + 2]);

                                // NCBI BLAST uses BLAST_WORD_THRESHOLD_TBLASTX = 13 for seed score threshold
                                // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:115
                                // #define BLAST_WORD_THRESHOLD_TBLASTX 13
                                // Also see: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1204
                                // For tblastx, it adds 2 to the base BLOSUM62 threshold (11) -> 13
                                if seed_score < 13 {
                                    if diag_enabled {
                                        diag_inner.base.seeds_low_score.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // Two-hit requirement: check if there's a previous seed on this diagonal
                                // within the window distance. NCBI BLAST style: when two hits occur,
                                // extend from the second hit (R) to the left and require it to reach
                                // the first hit (L) before doing the right extension.
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:535-551
                                let k_size = 3usize; // Word size for TBLASTX
                                let prev_s_pos_opt = last_seed.get(&mask_key).copied();
                                let two_hit_info = if let Some(prev_s_pos) = prev_s_pos_opt {
                                    let diff = s_pos.saturating_sub(prev_s_pos);
                                    // Check if distance is within window (diff <= window)
                                    // Reference: aa_ungapped.c:538
                                    if diff <= TWO_HIT_WINDOW {
                                        // Check if hits overlap (diff < wordsize) - skip if they do
                                        // Reference: aa_ungapped.c:549-551
                                        if diff < k_size {
                                            if diag_enabled {
                                                diag_inner.base.seeds_second_hit_overlap.fetch_add(1, AtomicOrdering::Relaxed);
                                            }
                                            None
                                        } else {
                                            if diag_enabled {
                                                diag_inner.base.seeds_second_hit_window.fetch_add(1, AtomicOrdering::Relaxed);
                                            }
                                            Some(prev_s_pos)
                                        }
                                    } else {
                                        if diag_enabled {
                                            diag_inner.base.seeds_second_hit_too_far.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        None
                                    }
                                } else {
                                    if diag_enabled {
                                        diag_inner.base.seeds_first_hit.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    None
                                };

                                // Update the last seed position for this diagonal
                                last_seed.insert(mask_key, s_pos);

                                // NCBI BLAST behavior: Two-hit requirement controls right extension execution,
                                // but left extension is always performed. If left extension score exceeds
                                // cutoff_score, the HSP is saved even when two-hit requirement is not met.
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:576-591
                                // - Two-hit satisfied: left + right extension → higher score HSP
                                // - Two-hit not satisfied: left extension only → left_score based HSP
                                // Both cases: if score >= cutoff_score, save as HSP
                                
                                // Count seeds that passed to extension
                                // Note: Even when two-hit requirement is not satisfied, extension is still performed
                                // (left extension only, which may still produce valid HSPs if score >= cutoff_score)
                                if diag_enabled {
                                    diag_inner.base.seeds_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                    diag_inner.base.ungapped_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                    if two_hit_info.is_none() {
                                        // Track one-hit extensions (left-only, two-hit not satisfied)
                                        diag_inner.base.ungapped_one_hit_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                    } else {
                                        // Track two-hit extensions (left+right, two-hit satisfied)
                                        diag_inner.base.ungapped_two_hit_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                }

                                let k_size = 3usize;
                                // Two-hit extension: extend from current seed (R) to the left,
                                // only extend right if left extension reaches the first hit (L)
                                // NCBI BLAST passes s_left_off = last_hit + wordsize (end of first word)
                                // When two-hit requirement is not met, use current seed position as s_left_off
                                // so that left extension cannot reach it, resulting in left-only extension
                                let (qs, qe, ss, se_ungapped, ungapped_score, right_extended, s_last_off) = {
                                    let s_left_off = if let Some(prev_s_pos) = two_hit_info {
                                        // Two-hit requirement satisfied: use first hit position
                                        prev_s_pos + k_size  // End of first hit (L + wordsize), NCBI BLAST style
                                    } else {
                                        // Two-hit requirement not satisfied: use current seed position
                                        // This ensures left extension cannot reach it, so only left extension is performed
                                        s_pos + k_size
                                    };
                                    
                                    let (qs, qe, ss, se, score, right_ext, s_last) = extend_hit_two_hit(
                                        q_aa,
                                        s_aa,
                                        s_left_off,
                                        s_pos,                          // Second hit position (R)
                                        q_pos as usize,
                                    );
                                    (qs, qe, ss, se, score, right_ext, s_last)
                                };

                                // NCBI BLAST does not use a fixed MIN_UNGAPPED_SCORE threshold.
                                // Instead, it uses cutoffs->cutoff_score which is dynamically calculated
                                // from E-value (CUTOFF_E_TBLASTX = 1e-300) via BLAST_Cutoffs function.
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:588
                                //            ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:360
                                //            ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:80
                                // Since we filter by E-value later, we don't need a fixed score threshold here.
                                // All hits will be filtered by E-value threshold, matching NCBI BLAST behavior.

                                // NCBI BLAST does not allow gapped search for TBLASTX
                                // Reference: ncbi-blast/c++/src/algo/blast/api/tblastx_options.cpp:73
                                // m_Opts->SetGappedMode(false);
                                // Also, BLAST_GAP_X_DROPOFF_TBLASTX = 0 and BLAST_GAP_X_DROPOFF_FINAL_TBLASTX = 0
                                // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:134,148
                                // Always use ungapped-only path to match NCBI BLAST behavior
                                if true {
                                    if diag_enabled {
                                        diag_inner.ungapped_only_hits.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    // Record the ungapped hit if it passes e-value threshold
                                    // NCBI BLAST uses query length and subject length (amino acid) for effective search space
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:842-843
                                    // effective_search_space = effective_db_length * (query_length - length_adjustment)
                                    // For TBLASTX, lengths are in amino acids (nucleotide length / 3)
                                    // NCBI BLAST calculates effective search space once per context (query frame + subject frame)
                                    // and reuses it for all HSPs in that context
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846
                                    // IMPORTANT: context_key must include subject frame index (s_f_idx) to match NCBI BLAST
                                    // Each context (query frame + subject frame combination) has its own effective search space
                                    let context_key = (q_idx, q_f_idx, s_f_idx as u8);
                                    let search_space = *search_space_cache.entry(context_key).or_insert_with(|| {
                                        // Use NCBI-style length adjustment for single-sequence comparison
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-824
                                        let ss = SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, &params);
                                        // Log effective search space for debugging (first time only per context)
                                        if diag_enabled {
                                            eprintln!("[DEBUG] Effective search space for context (q_frame={}, s_frame={}): effective_query_len={:.1}, effective_db_len={:.1}, effective_space={:.1}, length_adj={}", 
                                                q_frame.frame, s_frame.frame, ss.effective_query_len, ss.effective_db_len, ss.effective_space, ss.length_adjustment);
                                        }
                                        ss
                                    });
                                    
                                    // Calculate cutoff_score from CUTOFF_E_TBLASTX = 1e-300
                                    // NCBI BLAST uses this to filter extension results before E-value check
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:349-360
                                    //            ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:588
                                    // The search space used for cutoff_score calculation is MIN(q_len, s_len) * s_len
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:361-362
                                    // IMPORTANT: cutoff_score is capped by gap_trigger (BLAST_GAP_TRIGGER_PROT = 22.0 bit score)
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:366-367
                                    let cutoff_score = *cutoff_score_cache.entry(context_key).or_insert_with(|| {
                                        // Use MIN(q_len, s_len) * s_len for cutoff_score calculation
                                        // This matches NCBI BLAST's implementation
                                        let cutoff_search_space = SearchSpace {
                                            effective_query_len: (q_aa_len.min(s_aa_len)) as f64,
                                            effective_db_len: s_aa_len as f64,
                                            effective_space: ((q_aa_len.min(s_aa_len)) * s_aa_len) as f64,
                                            length_adjustment: 0,
                                        };
                                        let cutoff_from_evalue = raw_score_from_evalue(CUTOFF_E_TBLASTX, &params, &cutoff_search_space);
                                        
                                        // Calculate gap_trigger from bit score: gap_trigger = (gap_trigger_bit * ln(2) + ln(K)) / lambda
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:343-344
                                        use super::constants::GAP_TRIGGER_BIT_SCORE;
                                        use crate::stats::karlin::raw_score_from_bit_score;
                                        let gap_trigger = raw_score_from_bit_score(GAP_TRIGGER_BIT_SCORE, &params);
                                        
                                        // Use MIN(cutoff_from_evalue, gap_trigger) as NCBI BLAST does
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:367
                                        let cutoff = cutoff_from_evalue.min(gap_trigger);
                                        
                                        if diag_enabled {
                                            eprintln!("[DEBUG] Cutoff score for context (q_frame={}, s_frame={}): cutoff_e={}, cutoff_from_evalue={}, gap_trigger={}, final_cutoff_score={}", 
                                                q_frame.frame, s_frame.frame, CUTOFF_E_TBLASTX, cutoff_from_evalue, gap_trigger, cutoff);
                                        }
                                        cutoff
                                    });
                                    
                                    // NCBI BLAST checks cutoff_score AFTER extension, before E-value check
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:588
                                    if ungapped_score < cutoff_score {
                                        if diag_enabled {
                                            diag_inner.ungapped_cutoff_failed.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        // Update mask even if cutoff_score check fails (NCBI BLAST behavior)
                                        let mask_end = if right_extended {
                                            s_last_off.saturating_sub(k_size - 1)
                                        } else {
                                            s_pos
                                        };
                                        mask.insert(mask_key, mask_end);
                                        continue;
                                    }
                                    
                                    let bit_score = calc_bit_score(ungapped_score, &params);
                                    let e_val = calc_evalue(bit_score, &search_space);

                                    if e_val <= args.evalue {
                                        if diag_enabled {
                                            diag_inner.ungapped_evalue_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        let len = qe - qs;
                                        let mut match_count = 0;
                                        for k in 0..len {
                                            if q_aa[qs + k] == s_aa[ss + k] {
                                                match_count += 1;
                                            }
                                        }
                                        let identity = (match_count as f64 / len as f64) * 100.0;

                                        let (q_start_bp, q_end_bp) =
                                            convert_coords(qs, qe, q_frame.frame, q_frame.orig_len);
                                        let (s_start_bp, s_end_bp) =
                                            convert_coords(ss, se_ungapped, s_frame.frame, s_len);

                                        // Store sequence data for HSP chaining
                                        let seq_key: SequenceKey = (
                                            query_ids[q_idx as usize].clone(),
                                            s_id.clone(),
                                            q_frame.frame,
                                            s_frame.frame,
                                        );
                                        if !seen_pairs.contains(&seq_key) {
                                            seen_pairs.insert(seq_key.clone());
                                            local_sequences.push((
                                                seq_key.clone(),
                                                (q_aa.clone(), s_aa.clone()),
                                            ));
                                        }

                                        local_hits.push(ExtendedHit {
                                            hit: Hit {
                                                query_id: query_ids[q_idx as usize].clone(),
                                                subject_id: s_id.clone(),
                                                identity,
                                                length: len,
                                                mismatch: len - match_count,
                                                gapopen: 0,
                                                q_start: q_start_bp,
                                                q_end: q_end_bp,
                                                s_start: s_start_bp,
                                                s_end: s_end_bp,
                                                e_value: e_val,
                                                bit_score,
                                            },
                                            q_frame: q_frame.frame,
                                            s_frame: s_frame.frame,
                                            q_aa_start: qs,
                                            q_aa_end: qe,
                                            s_aa_start: ss,
                                            s_aa_end: se_ungapped,
                                            q_orig_len: q_frame.orig_len,
                                            s_orig_len: s_len,
                                            from_gapped: false,
                                        });
                                    } else if diag_enabled {
                                        diag_inner.ungapped_evalue_failed.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    
                                    // NCBI BLAST style diagonal suppression: update mask AFTER E-value check
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:593-606
                                    // IMPORTANT: NCBI BLAST updates diag_array[diag_coord].last_hit AFTER
                                    // the score check (line 588) and HSP saving (line 589-591).
                                    // This means mask is updated regardless of whether the HSP was saved or passed E-value.
                                    // The mask update happens after extension and score/E-value check, but before
                                    // the next iteration, ensuring that even low-scoring extensions prevent
                                    // re-extending in already-scanned regions.
                                    let mask_end = if right_extended {
                                        // Right extension happened: use s_last_off - (wordsize - 1)
                                        // Reference: aa_ungapped.c:598-599
                                        s_last_off.saturating_sub(k_size - 1)
                                    } else {
                                        // No right extension: use current seed position (subject_offset)
                                        // Reference: aa_ungapped.c:604-605
                                        // Note: In LOSAT, we don't have diag_offset separately, so we use s_pos directly
                                        s_pos
                                    };
                                    mask.insert(mask_key, mask_end);
                                    
                                    continue;
                                }

                                // Count gapped extensions
                                if diag_enabled {
                                    diag_inner.gapped_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                // Two-phase gapped extension (NCBI-style):
                                // Phase 1: Preliminary extension with lower X-drop
                                let (
                                    prelim_qs,
                                    prelim_qe,
                                    prelim_ss,
                                    prelim_se,
                                    prelim_score,
                                    _,
                                    _,
                                    _,
                                    _,
                                ) = extend_gapped_protein(
                                    q_aa,
                                    s_aa,
                                    qs,
                                    ss,
                                    qe - qs,
                                    X_DROP_GAPPED_PRELIM,
                                );

                                // Phase 2: Final extension with higher X-drop for longer alignments
                                let (
                                    final_qs,
                                    final_qe,
                                    final_ss,
                                    final_se,
                                    gapped_score,
                                    matches,
                                    mismatches,
                                    gap_opens,
                                    gap_letters,
                                ) = if prelim_score > 0 {
                                    extend_gapped_protein(
                                        q_aa,
                                        s_aa,
                                        prelim_qs,
                                        prelim_ss,
                                        prelim_qe - prelim_qs,
                                        X_DROP_GAPPED_FINAL,
                                    )
                                } else {
                                    (
                                        prelim_qs,
                                        prelim_qe,
                                        prelim_ss,
                                        prelim_se,
                                        prelim_score,
                                        0,
                                        0,
                                        0,
                                        0,
                                    )
                                };

                                mask.insert(mask_key, final_se);

                                // Calculate alignment length as matches + mismatches + gap_letters
                                // This is the correct BLAST definition of alignment length
                                let aln_len = matches + mismatches + gap_letters;

                                // NCBI BLAST uses query length and subject length (amino acid) for effective search space
                                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:842-843
                                // effective_search_space = effective_db_length * (query_length - length_adjustment)
                                // For TBLASTX, lengths are in amino acids (nucleotide length / 3)
                                // NCBI BLAST calculates effective search space once per context (query frame + subject frame)
                                // and reuses it for all HSPs in that context
                                // IMPORTANT: context_key must include subject frame index (s_f_idx) to match NCBI BLAST
                                let context_key = (q_idx, q_f_idx, s_f_idx as u8);
                                let search_space = *search_space_cache.entry(context_key).or_insert_with(|| {
                                    // Use NCBI-style length adjustment for single-sequence comparison
                                    SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, &params)
                                });
                                let bit_score = calc_bit_score(gapped_score, &params);
                                let e_val = calc_evalue(bit_score, &search_space);

                                if e_val <= args.evalue {
                                    if diag_enabled {
                                        diag_inner.gapped_evalue_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                    }

                                    // Identity is matches / alignment_length, capped at 100%
                                    let identity = if aln_len > 0 {
                                        ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                                    } else {
                                        0.0
                                    };

                                    let (q_start_bp, q_end_bp) = convert_coords(
                                        final_qs,
                                        final_qe,
                                        q_frame.frame,
                                        q_frame.orig_len,
                                    );
                                    let (s_start_bp, s_end_bp) =
                                        convert_coords(final_ss, final_se, s_frame.frame, s_len);

                                    // Store sequence data for HSP chaining
                                    let seq_key: SequenceKey = (
                                        query_ids[q_idx as usize].clone(),
                                        s_id.clone(),
                                        q_frame.frame,
                                        s_frame.frame,
                                    );
                                    if !seen_pairs.contains(&seq_key) {
                                        seen_pairs.insert(seq_key.clone());
                                        local_sequences
                                            .push((seq_key.clone(), (q_aa.clone(), s_aa.clone())));
                                    }

                                    local_hits.push(ExtendedHit {
                                        hit: Hit {
                                            query_id: query_ids[q_idx as usize].clone(),
                                            subject_id: s_id.clone(),
                                            identity,
                                            length: aln_len,
                                            mismatch: mismatches,
                                            gapopen: gap_opens,
                                            q_start: q_start_bp,
                                            q_end: q_end_bp,
                                            s_start: s_start_bp,
                                            s_end: s_end_bp,
                                            e_value: e_val,
                                            bit_score,
                                        },
                                        q_frame: q_frame.frame,
                                        s_frame: s_frame.frame,
                                        q_aa_start: final_qs,
                                        q_aa_end: final_qe,
                                        s_aa_start: final_ss,
                                        s_aa_end: final_se,
                                        q_orig_len: q_frame.orig_len,
                                        s_orig_len: s_len,
                                        from_gapped: true,
                                    });
                                }
                            }
                        }
                        (local_hits, local_sequences)
                    })
                    .collect();

            // Flatten frame results and send
            let mut all_hits: Vec<ExtendedHit> = Vec::new();
            let mut all_sequences: Vec<(SequenceKey, SequenceData)> = Vec::new();
            for (hits, seqs) in frame_results {
                all_hits.extend(hits);
                all_sequences.extend(seqs);
            }
            if !all_hits.is_empty() {
                tx.send((all_hits, all_sequences)).unwrap();
            }
            bar.inc(1);
        });

    bar.finish();
    writer_handle.join().unwrap()?;

    // Print diagnostic summary if enabled
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }

    Ok(())
}
