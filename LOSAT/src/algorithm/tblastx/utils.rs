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
use crate::utils::seg::SegMasker;

use super::args::TblastxArgs;
use super::constants::{
    CUTOFF_E_TBLASTX, GAP_EXTEND, GAP_OPEN, TWO_HIT_WINDOW, X_DROP_GAPPED_FINAL,
    X_DROP_GAPPED_PRELIM, X_DROP_UNGAPPED,
};
use super::chaining::ExtendedHit;
use super::diagnostics::{diagnostics_enabled, DiagnosticCounters, print_summary as print_diagnostics_summary};
use super::sum_stats_linking::apply_sum_stats_even_gap_linking;
use super::extension::{
    convert_coords, extend_gapped_protein, extend_hit_two_hit, extend_hit_ungapped,
    get_score,
};
use crate::stats::karlin::{bit_score as calc_bit_score, evalue as calc_evalue, raw_score_from_evalue};
use super::lookup::build_direct_lookup_with_threshold;
use super::translation::{generate_frames, QueryFrame};

pub fn run(args: TblastxArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let only_qframe = args.only_qframe;
    let only_sframe = args.only_sframe;
    let valid_frame = |f: i8| matches!(f, 1 | 2 | 3 | -1 | -2 | -3);
    if let Some(f) = only_qframe {
        if !valid_frame(f) {
            anyhow::bail!("--only-qframe must be one of 1,2,3,-1,-2,-3 (got {f})");
        }
    }
    if let Some(f) = only_sframe {
        if !valid_frame(f) {
            anyhow::bail!("--only-sframe must be one of 1,2,3,-1,-2,-3 (got {f})");
        }
    }
    
    // Determine two-hit window size based on arguments
    // Use --window_size to set (default: 40)
    // window_size = 0 enables one-hit mode (like NCBI BLAST's -window_size 0)
    let two_hit_window: usize = args.window_size;
    let one_hit_mode = two_hit_window == 0;
    
    // X-drop value (default is now NCBI-compatible: 7)
    let x_drop_ungapped: i32 = X_DROP_UNGAPPED;

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

    // Apply DUST filter to DNA sequences (if enabled)
    let mut query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        queries_raw
            .iter()
            .map(|r| masker.mask_sequence(r.seq()))
            .collect()
    } else {
        queries_raw.iter().map(|_| Vec::new()).collect()
    };

    // Generate query frames
    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| {
            let mut frames = generate_frames(r.seq(), &query_code);
            if let Some(f) = only_qframe {
                frames.retain(|x| x.frame == f);
            }
            frames
        })
        .collect();
    if let Some(_f) = only_qframe {
        // Safety: if the query is too short, a frame translation may be empty; prevent silent no-op searches.
        if query_frames.iter().all(|v| v.is_empty()) {
            anyhow::bail!("--only-qframe filtered out all query frames (query too short?)");
        }
    }

    // Apply SEG filter to translated amino acid sequences (if enabled)
    // SEG masks low-complexity regions by replacing amino acids with 'X' (code 23)
    // This causes extensions to naturally terminate at masked regions (X has low scores)
    // NCBI BLAST Reference: blast_filter.c - SEG replaces residues with 'X'
    const AA_X: u8 = 23; // NCBI encoding for 'X' (unknown amino acid)
    
    if args.seg {
        eprintln!(
            "Applying SEG filter (window={}, locut={}, hicut={})...",
            args.seg_window, args.seg_locut, args.seg_hicut
        );
        let seg_masker = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        
        let mut total_aa_masked = 0usize;
        let mut total_aa = 0usize;
        
        for (q_idx, frames) in query_frames.iter_mut().enumerate() {
            let orig_len = queries_raw[q_idx].seq().len();
            let mut seg_masks_dna: Vec<MaskedInterval> = Vec::new();
            
            for frame in frames.iter_mut() {
                // Apply SEG to amino acid sequence (skip sentinels at positions 0 and len-1)
                // aa_seq layout: [SENTINEL, aa0, aa1, ..., aaN-1, SENTINEL]
                // Extract actual amino acids (positions 1 to len-2)
                if frame.aa_seq.len() < 3 {
                    continue; // Need at least sentinel + 1 AA + sentinel
                }
                
                total_aa += frame.aa_len;
                
                // Apply SEG filter to get masked intervals
                let aa_seq_slice = &frame.aa_seq[1..frame.aa_seq.len() - 1];
                let aa_masks = seg_masker.mask_sequence(aa_seq_slice);
                
                // Replace masked amino acids with 'X' in the aa_seq
                // This causes extensions to naturally terminate at masked regions
                for aa_mask in &aa_masks {
                    // aa_mask coordinates are 0-indexed logical positions
                    // In aa_seq, position i corresponds to aa_seq[i + 1] (due to leading sentinel)
                    for pos in aa_mask.start..aa_mask.end {
                        if pos + 1 < frame.aa_seq.len() - 1 {
                            frame.aa_seq[pos + 1] = AA_X;
                            total_aa_masked += 1;
                        }
                    }
                }
                
                // Also convert AA mask coordinates to DNA coordinates for lookup table filtering
                for aa_mask in aa_masks {
                    let (dna_start, dna_end) = if frame.frame > 0 {
                        // Forward frame: frame 1,2,3 have offset 0,1,2
                        let offset = (frame.frame - 1) as usize;
                        let dna_start = aa_mask.start * 3 + offset;
                        let dna_end = (aa_mask.end * 3 + offset).min(orig_len);
                        (dna_start, dna_end)
                    } else {
                        // Reverse frame: -1,-2,-3
                        let offset = (-frame.frame - 1) as usize;
                        let dna_end = orig_len - (aa_mask.start * 3 + offset);
                        let dna_start = orig_len.saturating_sub(aa_mask.end * 3 + offset);
                        (dna_start.max(0), dna_end)
                    };
                    
                    if dna_start < dna_end {
                        seg_masks_dna.push(MaskedInterval::new(dna_start, dna_end));
                    }
                }
            }
            
            // Merge SEG masks with DUST masks for lookup table filtering
            if !seg_masks_dna.is_empty() {
                // Sort by start position
                seg_masks_dna.sort_by_key(|m| m.start);
                
                // Merge overlapping or adjacent intervals
                let mut merged: Vec<MaskedInterval> = Vec::new();
                for mask in seg_masks_dna {
                    if let Some(last) = merged.last_mut() {
                        if mask.start <= last.end + 1 {
                            last.end = last.end.max(mask.end);
                        } else {
                            merged.push(mask);
                        }
                    } else {
                        merged.push(mask);
                    }
                }
                
                // Merge with existing DUST masks
                let mut all_masks = query_masks[q_idx].clone();
                all_masks.extend(merged);
                
                // Sort and merge all masks together
                all_masks.sort_by_key(|m| m.start);
                let mut final_masks: Vec<MaskedInterval> = Vec::new();
                for mask in all_masks {
                    if let Some(last) = final_masks.last_mut() {
                        if mask.start <= last.end + 1 {
                            last.end = last.end.max(mask.end);
                        } else {
                            final_masks.push(mask);
                        }
                    } else {
                        final_masks.push(mask);
                    }
                }
                
                query_masks[q_idx] = final_masks;
            }
        }
        
        // Print statistics
        let total_masked_dna: usize = query_masks.iter().map(|m| m.iter().map(|i| i.end - i.start).sum::<usize>()).sum();
        let total_bases: usize = queries_raw.iter().map(|r| r.seq().len()).sum();
        if total_bases > 0 {
            eprintln!(
                "SEG masked {} aa ({:.2}%) of {} total aa across all frames",
                total_aa_masked,
                100.0 * total_aa_masked as f64 / total_aa as f64,
                total_aa
            );
            eprintln!(
                "Combined DNA masks: {} bases ({:.2}%)",
                total_masked_dna,
                100.0 * total_masked_dna as f64 / total_bases as f64
            );
        }
    }

    eprintln!("Building lookup table...");
    // NCBI BLAST+ tblastx default: Neighboring words threshold = 13
    // (confirmed via pairwise output header: "Neighboring words threshold: 13")
    let lookup = build_direct_lookup_with_threshold(&query_frames, &query_masks, args.threshold);

    eprintln!("Reading subject file...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> =
        subject_reader.records().filter_map(|r| r.ok()).collect();

    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    let db_len_total_bp: usize = subjects_raw.iter().map(|r| r.seq().len()).sum();
    let db_len_aa_total = db_len_total_bp / 3;

    // TBLASTX uses UNGAPPED alignments, so we need ungapped Karlin-Altschul parameters
    // Reference: ncbi-blast/c++/src/algo/blast/api/tblastx_options.cpp:73
    // m_Opts->SetGappedMode(false);
    // For ungapped BLOSUM62, use gap_open = i32::MAX, gap_extend = i32::MAX
    // which corresponds to lambda=0.3176, K=0.134 (from blast_stat.c)
    let scoring_spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: i32::MAX, // Ungapped parameters
        gap_extend: i32::MAX, // Ungapped parameters
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

    // Channel sends ExtendedHit records for post-processing (sum-statistics linking, filtering, output)
    type HitData = Vec<ExtendedHit>;
    let (tx, rx) = channel::<HitData>();

    let out_path = args.out.clone();
    let params_clone = params.clone();
    let evalue_threshold = args.evalue;
    let diagnostics_clone = diagnostics.clone();
    let diag_enabled_clone = diag_enabled;
    let writer_handle = std::thread::spawn(move || -> Result<()> {
        let mut all_hits: Vec<ExtendedHit> = Vec::new();

        while let Ok(hits) = rx.recv() {
            all_hits.extend(hits);
        }

        let diag_ref = if diag_enabled_clone { Some(diagnostics_clone.as_ref()) } else { None };
        
        // Record HSPs before filtering
        let hsps_before = all_hits.len();
        if let Some(diag) = diag_ref {
            diag.base.hsps_before_chain.store(hsps_before, AtomicOrdering::Relaxed);
        }

        // === NCBI-like sum-statistics behavior ===
        // NCBI tblastx assigns a *set-level* E-value to each HSP after linking (BLAST_LinkHsps).
        // We apply an even-gap (small/large) linking step and overwrite `hit.e_value` with the
        // resulting set-level E-value.
        let linked_hits: Vec<ExtendedHit> =
            apply_sum_stats_even_gap_linking(all_hits, &params_clone);

        // Record HSPs after "linking" (before threshold)
        if let Some(diag) = diag_ref {
            diag.base.hsps_after_chain.store(linked_hits.len(), AtomicOrdering::Relaxed);
            diag.base.hsps_culled_dominated.store(0, AtomicOrdering::Relaxed);
        }

        // Filter by E-value threshold (now using set-level e-values)
        let mut output_from_ungapped = 0usize;
        let mut output_from_gapped = 0usize;
        let mut ungapped_evalue_passed = 0usize;
        let mut ungapped_evalue_failed = 0usize;
        let mut filtered_hits: Vec<Hit> = Vec::new();
        filtered_hits.reserve(linked_hits.len());
        for ext_hit in linked_hits.into_iter() {
            if ext_hit.hit.e_value <= evalue_threshold {
                ungapped_evalue_passed += 1;
                if ext_hit.from_gapped {
                    output_from_gapped += 1;
                } else {
                    output_from_ungapped += 1;
                }
                filtered_hits.push(ext_hit.hit);
            } else {
                ungapped_evalue_failed += 1;
            }
        }

        if let Some(diag) = diag_ref {
            diag.ungapped_evalue_passed
                .store(ungapped_evalue_passed, AtomicOrdering::Relaxed);
            diag.ungapped_evalue_failed
                .store(ungapped_evalue_failed, AtomicOrdering::Relaxed);
        }

        // Sort by bit score (highest first) for consistent output order
        filtered_hits.sort_by(|a, b| {
            b.bit_score
                .partial_cmp(&a.bit_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

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
            let mut s_frames = generate_frames(s_record.seq(), &db_code);
            if let Some(f) = only_sframe {
                s_frames.retain(|x| x.frame == f);
            }
            let s_id = s_record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            let s_len = s_record.seq().len();
            
            // Detect self-comparison (query == subject) for diagnostics
            // This is a common source of unexpectedly long alignments
            let is_self_comparison = query_ids.iter().any(|q_id| q_id == &s_id);
            if is_self_comparison && diag_enabled {
                diag.base.self_comparisons.fetch_add(1, AtomicOrdering::Relaxed);
            }

            // Parallelize over subject frames for intra-subject parallelism
            // Each frame is independent since mask_key includes s_f_idx
            let diag_inner = diag.clone();
            let frame_results: Vec<Vec<ExtendedHit>> =
                s_frames
                    .par_iter()
                    .enumerate()
                    .map(|(s_f_idx, s_frame)| {
                        let mut local_hits: Vec<ExtendedHit> = Vec::new();

                        // Mask to track already-extended regions on each diagonal (per frame)
                        // NCBI BLAST style: use diag_coord = (query_offset - subject_offset) & diag_mask
                        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:351
                        //            ncbi-blast/c++/src/algo/blast/core/blast_extend.c:42-61
                        let mut mask: FxHashMap<(u32, u8, isize), usize> = FxHashMap::default();
                        // Two-hit tracking: stores the last seed position on each diagonal for each query/frame combination
                        // NCBI BLAST style: use diag_coord instead of raw diagonal
                        let mut last_seed: FxHashMap<(u32, u8, isize), usize> =
                            FxHashMap::default();
                        // NCBI BLAST flag mechanism: tracks whether an extension just happened on a diagonal
                        // Reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:57-60
                        // When flag = true, the next hit on this diagonal that is PAST the extension end
                        // becomes a new first hit (not a two-hit pair).
                        let mut diag_extended: FxHashMap<(u32, u8, isize), bool> =
                            FxHashMap::default();
                        let s_aa = &s_frame.aa_seq;
                        // Use aa_len (actual amino acid count without sentinels) for calculations
                        let s_aa_len = s_frame.aa_len;
                        
                        // Calculate diag_mask for each query frame (NCBI BLAST style)
                        // diag_array_length = smallest power of 2 >= (query_length + window_size)
                        // diag_mask = diag_array_length - 1
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
                        let mut diag_mask_cache: FxHashMap<(u32, u8), isize> = FxHashMap::default();
                        for (q_idx, q_frames) in query_frames.iter().enumerate() {
                            for (q_f_idx, q_frame) in q_frames.iter().enumerate() {
                                // Use aa_len (actual amino acid count without sentinels) for calculations
                                let q_aa_len = q_frame.aa_len;
                                // Calculate diag_array_length as smallest power of 2 >= (q_aa_len + two_hit_window)
                                let mut diag_array_length = 1usize;
                                while diag_array_length < (q_aa_len + two_hit_window) {
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
                        // aa_seq layout: [SENTINEL, aa0, aa1, ..., aaN-1, SENTINEL]
                        // Need at least 5 bytes: sentinel + 3 AAs + sentinel
                        if s_aa.len() < 5 {
                            return local_hits;
                        }

                        // Scan from raw position 1 to len-4 (skip sentinels)
                        // s_pos is in RAW coordinates (index into aa_seq with sentinels)
                        for s_pos in 1..=(s_aa.len() - 4) {
                            let c1 = unsafe { *s_aa.get_unchecked(s_pos) } as usize;
                            let c2 = unsafe { *s_aa.get_unchecked(s_pos + 1) } as usize;
                            let c3 = unsafe { *s_aa.get_unchecked(s_pos + 2) } as usize;

                            // Skip k-mers containing stop codons (24) or sentinels (255)
                            // NCBI matrix order: 0-23 for amino acids, 24 for stop codon (*)
                            if c1 >= 24 || c2 >= 24 || c3 >= 24 {
                                continue;
                            }
                            // 24^3 = 13,824 possible k-mers (matching lookup.rs encoding)
                            let kmer = c1 * 576 + c2 * 24 + c3;

                            let matches = unsafe { lookup.get_unchecked(kmer) };

                            // Count k-mer matches for diagnostics
                            if diag_enabled && !matches.is_empty() {
                                diag_inner.base.kmer_matches.fetch_add(matches.len(), AtomicOrdering::Relaxed);
                            }

                            for &(q_idx, q_f_idx, q_pos) in matches {
                                // q_pos is LOGICAL position (0-indexed, not counting sentinels)
                                // s_pos is RAW position (index into aa_seq with sentinels)
                                // Convert q_pos to raw for consistent diagonal calculation
                                let q_pos_raw = q_pos as usize + 1;
                                
                                // NCBI BLAST style: diag_coord = (query_offset - subject_offset) & diag_mask
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:351
                                // Both positions are now in RAW coordinates for consistent diagonal
                                let diag = q_pos_raw as isize - s_pos as isize;
                                // Get diag_mask for this query frame
                                let diag_mask = *diag_mask_cache.get(&(q_idx, q_f_idx))
                                    .unwrap_or(&((1isize << 20) - 1)); // Fallback: 2^20 - 1
                                // Calculate diag_coord using NCBI BLAST style masking
                                let diag_coord = diag & diag_mask;
                                // Use diag_coord instead of raw diagonal for mask_key (NCBI BLAST compatible)
                                let mask_key = (q_idx, q_f_idx, diag_coord);

                                // Check if this region was already extended (mask check)
                                // NCBI BLAST style diagonal suppression
                                // Note: NCBI BLAST uses strict < comparison (subject_offset + diag_offset < last_hit)
                                // This allows a hit exactly AT the extension end to start a new sequence
                                if let Some(&last_end) = mask.get(&mask_key) {
                                    if s_pos < last_end {
                                        if diag_enabled {
                                            diag_inner.base.seeds_masked.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        continue;
                                    }
                                }

                                let q_frame = &query_frames[q_idx as usize][q_f_idx as usize];
                                let q_aa = &q_frame.aa_seq;
                                // Use aa_len (actual amino acid count without sentinels) for calculations
                                let q_aa_len = q_frame.aa_len;

                                // Access q_aa using RAW position (q_pos_raw = q_pos + 1)
                                let seed_score = get_score(c1 as u8, q_aa[q_pos_raw])
                                    + get_score(c2 as u8, q_aa[q_pos_raw + 1])
                                    + get_score(c3 as u8, q_aa[q_pos_raw + 2]);

                                // Early filtering: skip very low-scoring seeds (f593910 optimization)
                                if seed_score < 11 {
                                    if diag_enabled {
                                        diag_inner.base.seeds_low_score.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // f593910-style two-hit mechanism (simpler and faster)
                                let k_size = 3usize;
                                
                                // NCBI BLAST flag check: if an extension just happened on this diagonal,
                                // skip seeds that are before the extension end
                                // Reference: aa_ungapped.c:354-365, 519-530
                                if let Some(&extended) = diag_extended.get(&mask_key) {
                                    if extended {
                                        // An extension happened on this diagonal
                                        if let Some(&last_hit) = last_seed.get(&mask_key) {
                                            if s_pos < last_hit {
                                                // We've already extended past this hit, skip it
                                                if diag_enabled {
                                                    diag_inner.base.seeds_masked.fetch_add(1, AtomicOrdering::Relaxed);
                                                }
                                                continue;
                                            } else {
                                                // Past the extension end - start a new hit
                                                // Reset flag and update last_hit
                                                diag_extended.insert(mask_key, false);
                                                last_seed.insert(mask_key, s_pos);
                                            }
                                        }
                                    }
                                }
                                
                                let prev_s_pos_opt = last_seed.get(&mask_key).copied();
                                
                                // NCBI BLAST style two-hit logic
                                // Reference: aa_ungapped.c:532-551
                                let (two_hit_info, should_update_last_seed) = if one_hit_mode {
                                    // One-hit mode: bypass two-hit check
                                    (None, true)
                                } else if let Some(prev_s_pos) = prev_s_pos_opt {
                                    let diff = s_pos.saturating_sub(prev_s_pos);
                                    
                                    if diff >= two_hit_window {
                                        // Beyond the window - start a new hit
                                        // Reference: aa_ungapped.c:538-543
                                        if diag_enabled {
                                            diag_inner.base.seeds_second_hit_too_far.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        (None, true)
                                    } else if diff < k_size {
                                        // Overlapping hits (diff < wordsize) - give up
                                        // Reference: aa_ungapped.c:546-551
                                        // IMPORTANT: Do NOT update last_seed in this case
                                        if diag_enabled {
                                            diag_inner.base.seeds_second_hit_overlap.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        continue;
                                    } else {
                                        // Within two-hit window and not overlapping
                                        // k_size <= diff < two_hit_window
                                        // Proceed to extension
                                        if diag_enabled {
                                            diag_inner.base.seeds_second_hit_window.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        (Some(prev_s_pos), true)
                                    }
                                } else {
                                    // No previous hit on this diagonal - first hit
                                    if diag_enabled {
                                        diag_inner.base.seeds_first_hit.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    (None, true)
                                };
                                
                                // Update the last seed position for this diagonal
                                // (only if should_update_last_seed is true)
                                if should_update_last_seed {
                                    last_seed.insert(mask_key, s_pos);
                                }
                                
                                // NCBI BLAST two-hit mode: only extend when two-hit condition is met
                                // Unlike the previous "f593910 optimization", NCBI does NOT extend
                                // high-scoring seeds without a two-hit pair in two-hit mode.
                                // Reference: aa_ungapped.c - no one-hit extension in s_BlastAaWordFinder_TwoHit
                                if two_hit_info.is_none() && !one_hit_mode {
                                    if diag_enabled {
                                        diag_inner.base.seeds_two_hit_failed.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // Count seeds that passed to extension
                                if diag_enabled {
                                    diag_inner.base.seeds_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                    diag_inner.base.ungapped_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                    if one_hit_mode {
                                        diag_inner.base.ungapped_one_hit_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                    } else {
                                        diag_inner.base.ungapped_two_hit_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                }
                                
                                // Extension: f593910 style - use two_hit_info to decide mode
                                // All positions passed to extension functions are RAW coordinates
                                // Results (qs_raw, qe_raw, ss_raw, se_raw) are also RAW coordinates
                                let (qs_raw, qe_raw, ss_raw, se_ungapped_raw, ungapped_score, right_extended, s_last_off) = if let Some(prev_s_pos) = two_hit_info {
                                    // Two-hit extension: extend from current seed (R) to the left,
                                    // only extend right if left extension reaches the first hit (L)
                                    let (qs, qe, ss, se, score, right_ext, s_last) = extend_hit_two_hit(
                                        q_aa,
                                        s_aa,
                                        prev_s_pos + k_size, // End of first hit (L + wordsize) - RAW
                                        s_pos,               // Second hit position (R) - RAW
                                        q_pos_raw,           // Query position - RAW
                                        x_drop_ungapped,
                                    );
                                    (qs, qe, ss, se, score, right_ext, s_last)
                                } else {
                                    // One-hit extension for high-scoring seeds (or one_hit_mode)
                                    let (qs, qe, ss, se, score, s_last) = extend_hit_ungapped(
                                        q_aa,
                                        s_aa,
                                        q_pos_raw,           // Query position - RAW
                                        s_pos,               // Subject position - RAW
                                        seed_score,
                                        x_drop_ungapped,
                                    );
                                    (qs, qe, ss, se, score, true, s_last)
                                };
                                
                                // Convert RAW coordinates to LOGICAL coordinates (subtract 1 for sentinel offset)
                                // Keep both RAW and LOGICAL for different uses:
                                // - RAW: mask updates, array access (q_aa, s_aa)
                                // - LOGICAL: coordinate conversion, length calculation
                                let qs = qs_raw - 1;
                                let qe = qe_raw - 1;
                                let ss = ss_raw - 1;
                                let se_ungapped = se_ungapped_raw - 1;

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
                                        
                                        // Calculate cutoff_score_max FIRST (before gap_trigger)
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:372-373, 943-946
                                        // For TBLASTX, cutoff_score_max is set in BlastHitSavingParametersUpdate
                                        // and is calculated from expect_value (default 10.0) using eff_searchsp
                                        //
                                        // IMPORTANT: NCBI BLAST uses eff_searchsp (with length adjustment) for cutoff_score_max,
                                        // NOT the same search space as cutoff_score calculation.
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:926
                                        // searchsp = query_info->contexts[context].eff_searchsp;
                                        // This is the key difference - eff_searchsp is smaller, making cutoff_score_max smaller
                                        // and thus allowing lower-scoring hits to pass.
                                        // Also, cutoff_score_max is scaled by scale_factor after calculation.
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:985
                                        // params->cutoffs[context].cutoff_score_max *= (Int4) scale_factor;
                                        //
                                        // IMPORTANT: The user's insight is that cutoff_score_max should allow hits with
                                        // "evalue < 10 OR bitscore > threshold" (where threshold is gap_trigger = 22.0 bit score).
                                        // This means cutoff_score_max should be calculated to allow hits with bit score >= 22.0,
                                        // which corresponds to gap_trigger raw score (46). So cutoff_score_max should be
                                        // at least gap_trigger, or calculated from an evalue that corresponds to gap_trigger.
                                        let scale_factor = 1.0; // Default for TBLASTX (RPS-BLAST uses > 1.0)
                                        
                                        // Calculate gap_trigger first to determine the minimum cutoff_score_max
                                        use super::constants::GAP_TRIGGER_BIT_SCORE;
                                        use crate::stats::karlin::raw_score_from_bit_score;
                                        let gap_trigger = raw_score_from_bit_score(GAP_TRIGGER_BIT_SCORE, &params);
                                        
                                        // Calculate cutoff_score_max from args.evalue (default 10.0)
                                        let cutoff_score_max_from_evalue = (raw_score_from_evalue(args.evalue, &params, &search_space) as f64 * scale_factor) as i32;
                                        
                                        // NCBI BLAST uses gap_trigger (22 bit score ≈ 46 raw score) as minimum cutoff
                                        // This ensures short hits (like 3 aa) with low scores are filtered out
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:366-373
                                        let cutoff_score_max = gap_trigger as i32;
                                        
                                        // gap_trigger is already calculated above
                                        
                                        // IMPORTANT: NCBI BLAST applies gap_trigger BEFORE cutoff_score_max, but
                                        // if cutoff_score_max < gap_trigger, cutoff_score_max will override gap_trigger.
                                        // However, if cutoff_score_max >= gap_trigger, gap_trigger is applied first.
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:103-110, 372-373
                                        // The actual NCBI BLAST code always applies gap_trigger first, then cutoff_score_max:
                                        //   new_cutoff = MIN(new_cutoff, gap_trigger);
                                        //   new_cutoff *= scale_factor;
                                        //   new_cutoff = MIN(new_cutoff, cutoff_score_max);
                                        // This means if cutoff_score_max >= gap_trigger, gap_trigger is the limiting factor.
                                        // But the user wants cutoff_score_max to take precedence when it's larger.
                                        // IMPORTANT: Since cutoff_score_max = gap_trigger (to allow "evalue < 10 OR bitscore > threshold"),
                                        // we should use cutoff_score_max directly, not gap_trigger.
                                        // This ensures that hits with bit score >= 22.0 (gap_trigger) can pass even if evalue >= 10.0.
                                        let mut cutoff = cutoff_from_evalue.min(cutoff_score_max);
                                        
                                        // Apply scale_factor (usually 1.0 for TBLASTX)
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:371
                                        // Note: scale_factor is already applied to cutoff_score_max above
                                        cutoff = (cutoff as f64 * scale_factor) as i32;
                                        
                                        // Apply cutoff_score_max limit (final check)
                                        // This ensures cutoff_score_max is always respected
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:372-373
                                        let cutoff_before_max = cutoff;
                                        cutoff = cutoff.min(cutoff_score_max);
                                        
                                        if diag_enabled {
                                            let limiting_factor = if cutoff == cutoff_score_max && cutoff_score_max == gap_trigger {
                                                "cutoff_score_max(gap_trigger)"
                                            } else if cutoff == cutoff_score_max {
                                                "cutoff_score_max"
                                            } else if cutoff == gap_trigger {
                                                "gap_trigger"
                                            } else if cutoff == cutoff_from_evalue {
                                                "cutoff_from_evalue"
                                            } else {
                                                "unknown"
                                            };
                                            eprintln!("[DEBUG] Cutoff score for context (q_frame={}, s_frame={}): cutoff_e={}, cutoff_from_evalue={}, gap_trigger={}, scale_factor={}, cutoff_score_max={}, cutoff_before_max={}, final_cutoff_score={}, limiting_factor={}", 
                                                q_frame.frame, s_frame.frame, CUTOFF_E_TBLASTX, cutoff_from_evalue, gap_trigger, scale_factor, cutoff_score_max, cutoff_before_max, cutoff, limiting_factor);
                                        }
                                        cutoff
                                    });
                                    
                                    // NCBI BLAST checks cutoff_score AFTER extension, before E-value check
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:588
                                    if ungapped_score < cutoff_score {
                                        if diag_enabled {
                                            diag_inner.ungapped_cutoff_failed.fetch_add(1, AtomicOrdering::Relaxed);
                                            // Track score range for cutoff_failed hits
                                            use std::sync::atomic::Ordering as AtomicOrdering;
                                            // Update minimum score
                                            let mut current_min = diag_inner.ungapped_cutoff_failed_min_score.load(AtomicOrdering::Relaxed);
                                            while ungapped_score < current_min {
                                                match diag_inner.ungapped_cutoff_failed_min_score.compare_exchange_weak(
                                                    current_min, ungapped_score, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed
                                                ) {
                                                    Ok(_) => break,
                                                    Err(x) => current_min = x,
                                                }
                                            }
                                            // Update maximum score
                                            let mut current_max = diag_inner.ungapped_cutoff_failed_max_score.load(AtomicOrdering::Relaxed);
                                            while ungapped_score > current_max {
                                                match diag_inner.ungapped_cutoff_failed_max_score.compare_exchange_weak(
                                                    current_max, ungapped_score, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed
                                                ) {
                                                    Ok(_) => break,
                                                    Err(x) => current_max = x,
                                                }
                                            }
                                        }
                                        // Update mask even if cutoff_score check fails (NCBI BLAST behavior)
                                        // Use RAW coordinates for mask (se_ungapped_raw) to properly suppress the aligned region
                                        let mask_end = if right_extended {
                                            se_ungapped_raw.max(s_last_off.saturating_sub(k_size - 1))
                                        } else {
                                            s_pos  // s_pos is already RAW
                                        };
                                        mask.insert(mask_key, mask_end);
                                        // NCBI BLAST: Update last_hit and flag based on right_extended
                                        // Reference: aa_ungapped.c:596-606
                                        if right_extended {
                                            // Set flag to indicate extension happened
                                            // Next hit must be past mask_end to start new sequence
                                            diag_extended.insert(mask_key, true);
                                            last_seed.insert(mask_key, mask_end);
                                        } else {
                                            // No right extension: update last_hit to current position
                                            last_seed.insert(mask_key, s_pos);
                                        }
                                        if diag_enabled {
                                            diag_inner.base.mask_updates.fetch_add(1, AtomicOrdering::Relaxed);
                                        }
                                        continue;
                                    }
                                    
                                    let bit_score = calc_bit_score(ungapped_score, &params);
                                    let e_val = calc_evalue(bit_score, &search_space);

                                    // IMPORTANT (NCBI behavior):
                                    // At this stage, NCBI BLAST only applies the raw-score cutoff (cutoffs->cutoff_score)
                                    // and does NOT filter by the user E-value threshold yet. E-values are assigned later
                                    // (and may be modified by sum-statistics linking), and only then are HSPs reaped by
                                    // the E-value threshold.
                                    //
                                    // Therefore, we always save HSPs that pass cutoff_score here, and defer the E-value
                                    // threshold filtering to the writer/post-processing stage.
                                    let len = qe - qs;
                                    let mut match_count = 0;
                                    // Use RAW coordinates for array access (qs_raw, ss_raw)
                                    for k in 0..len {
                                        if q_aa[qs_raw + k] == s_aa[ss_raw + k] {
                                            match_count += 1;
                                        }
                                    }
                                    let identity = (match_count as f64 / len as f64) * 100.0;

                                    let (q_start_bp, q_end_bp) =
                                        convert_coords(qs, qe, q_frame.frame, q_frame.orig_len);
                                    let (s_start_bp, s_end_bp) =
                                        convert_coords(ss, se_ungapped, s_frame.frame, s_len);

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
                                            // Temporary individual E-value; may be overwritten by sum-statistics linking.
                                            e_value: e_val,
                                            bit_score,
                                        },
                                        raw_score: ungapped_score,
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
                                    
                                    // NCBI BLAST style diagonal suppression: update mask AFTER E-value check
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:593-606
                                    // IMPORTANT: NCBI BLAST updates diag_array[diag_coord].last_hit AFTER
                                    // the score check (line 588) and HSP saving (line 589-591).
                                    // This means mask is updated regardless of whether the HSP was saved or passed E-value.
                                    // The mask update happens after extension and score/E-value check, but before
                                    // the next iteration, ensuring that even low-scoring extensions prevent
                                    // re-extending in already-scanned regions.
                                    //
                                    // FIX: Use the actual alignment end (se_ungapped_raw) as the primary mask position.
                                    // Use RAW coordinates for mask (consistent with s_pos which is RAW).
                                    let mask_end = if right_extended {
                                        // Right extension happened: use the actual alignment end position (RAW)
                                        // to properly suppress the entire aligned region.
                                        // Also consider s_last_off (rightmost scanned position) for safety.
                                        se_ungapped_raw.max(s_last_off.saturating_sub(k_size - 1))
                                    } else {
                                        // No right extension: use current seed position (subject_offset)
                                        // Reference: aa_ungapped.c:604-605
                                        s_pos  // s_pos is already RAW
                                    };
                                    mask.insert(mask_key, mask_end);
                                    // NCBI BLAST: Update last_hit and flag based on right_extended
                                    // Reference: aa_ungapped.c:596-606
                                    if right_extended {
                                        // Set flag to indicate extension happened
                                        // Next hit must be past mask_end to start new sequence
                                        diag_extended.insert(mask_key, true);
                                        last_seed.insert(mask_key, mask_end);
                                    } else {
                                        // No right extension: update last_hit to current position
                                        last_seed.insert(mask_key, s_pos);
                                    }
                                    if diag_enabled {
                                        diag_inner.base.mask_updates.fetch_add(1, AtomicOrdering::Relaxed);
                                        // Track extension length
                                        let ext_len = qe.saturating_sub(qs);
                                        diag_inner.base.extension_total_length.fetch_add(ext_len, AtomicOrdering::Relaxed);
                                        // Update max length (using compare-and-swap for thread safety)
                                        let mut current_max = diag_inner.base.extension_max_length.load(AtomicOrdering::Relaxed);
                                        while ext_len > current_max {
                                            match diag_inner.base.extension_max_length.compare_exchange(
                                                current_max, ext_len, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed
                                            ) {
                                                Ok(_) => break,
                                                Err(x) => current_max = x,
                                            }
                                        }
                                    }
                                    
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
                                    qs_raw,  // Use RAW coordinates for gapped extension
                                    ss_raw,
                                    qe_raw - qs_raw,
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
                                        raw_score: gapped_score,
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
                        local_hits
                    })
                    .collect();

            // Flatten frame results and send
            let mut all_hits: Vec<ExtendedHit> = Vec::new();
            for hits in frame_results {
                all_hits.extend(hits);
            }
            if !all_hits.is_empty() {
                tx.send(all_hits).unwrap();
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
