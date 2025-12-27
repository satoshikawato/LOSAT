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
    GAP_EXTEND, GAP_OPEN, MIN_UNGAPPED_SCORE, TWO_HIT_WINDOW, X_DROP_GAPPED_FINAL,
    X_DROP_GAPPED_PRELIM,
};
use super::chaining::{chain_and_filter_hsps_protein, ExtendedHit, SequenceData, SequenceKey};
use super::diagnostics::{diagnostics_enabled, DiagnosticCounters, print_summary as print_diagnostics_summary};
use super::extension::{
    convert_coords, extend_gapped_protein, extend_hit_two_hit,
    get_score,
};
use crate::stats::karlin::{bit_score as calc_bit_score, evalue as calc_evalue};
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
        let mut all_sequences: FxHashMap<SequenceKey, SequenceData> = FxHashMap::default();

        while let Ok((hits, seq_data)) = rx.recv() {
            all_hits.extend(hits);
            for (key, data) in seq_data {
                all_sequences.entry(key).or_insert(data);
            }
        }

        // Chain nearby HSPs into longer alignments using cluster-then-extend
        let diag_ref = if diag_enabled_clone { Some(diagnostics_clone.as_ref()) } else { None };
        let filtered_hits = chain_and_filter_hsps_protein(
            all_hits,
            &all_sequences,
            db_len_aa_total,
            &params_clone,
            evalue_threshold,
            diag_ref,
        );

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
                    .map(|(_s_f_idx, s_frame)| {
                        let mut local_hits: Vec<ExtendedHit> = Vec::new();
                        let mut local_sequences: Vec<(SequenceKey, SequenceData)> = Vec::new();
                        let mut seen_pairs: std::collections::HashSet<SequenceKey> =
                            std::collections::HashSet::new();

                        // Mask to track already-extended regions on each diagonal (per frame)
                        let mut mask: FxHashMap<(u32, u8, isize), usize> = FxHashMap::default();
                        // Two-hit tracking: stores the last seed position on each diagonal for each query/frame combination
                        let mut last_seed: FxHashMap<(u32, u8, isize), usize> =
                            FxHashMap::default();
                        let s_aa = &s_frame.aa_seq;
                        let s_aa_len = s_aa.len();
                        
                        // Cache effective search space per query-subject-frame combination
                        // NCBI BLAST calculates this once per context (query frame + subject frame)
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846
                        let mut search_space_cache: FxHashMap<(u32, u8), SearchSpace> = FxHashMap::default();
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
                                let diag = s_pos as isize - q_pos as isize;
                                // Simplified mask_key since we're now per-frame (s_f_idx is implicit)
                                let mask_key = (q_idx, q_f_idx, diag);

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
                                let prev_s_pos_opt = last_seed.get(&mask_key).copied();
                                let two_hit_info = if let Some(prev_s_pos) = prev_s_pos_opt {
                                    if s_pos.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW {
                                        Some(prev_s_pos)
                                    } else {
                                        None
                                    }
                                } else {
                                    None
                                };

                                // Update the last seed position for this diagonal
                                last_seed.insert(mask_key, s_pos);

                                // NCBI BLAST uses strict two-hit requirement for TBLASTX
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:214-232
                                // When multiple_hits is true (window_size > 0), s_BlastAaWordFinder_TwoHit is used
                                // which requires two hits within the window before extension
                                // No one-hit extension is allowed, even for high-scoring seeds
                                if two_hit_info.is_none() {
                                    if diag_enabled {
                                        diag_inner.base.seeds_two_hit_failed.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // Count seeds that passed to extension
                                if diag_enabled {
                                    diag_inner.base.seeds_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                    diag_inner.base.ungapped_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                // NCBI BLAST always uses two-hit extension for TBLASTX
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:576-583
                                let k_size = 3usize;
                                // Two-hit extension: extend from current seed (R) to the left,
                                // only extend right if left extension reaches the first hit (L)
                                // NCBI BLAST passes s_left_off = last_hit + wordsize (end of first word)
                                let (qs, qe, ss, se_ungapped, ungapped_score, s_last_off) = {
                                    let (qs, qe, ss, se, score, _right_extended, s_last) = extend_hit_two_hit(
                                        q_aa,
                                        s_aa,
                                        two_hit_info.unwrap() + k_size,  // End of first hit (L + wordsize), NCBI BLAST style
                                        s_pos,                          // Second hit position (R)
                                        q_pos as usize,
                                    );
                                    (qs, qe, ss, se, score, s_last)
                                };

                                // NCBI BLAST style diagonal suppression: update mask after EVERY extension,
                                // using s_last_off (rightmost position scanned), not se_ungapped (best-scoring end).
                                // This prevents re-extending inside already-scanned regions and reduces fragmentation.
                                // NCBI BLAST uses: s_last_off - (wordsize - 1) for the mask value
                                // IMPORTANT: Update mask BEFORE score check, so even low-scoring extensions
                                // prevent re-extending in already-scanned regions (matches NCBI BLAST behavior)
                                let mask_end = s_last_off.saturating_sub(k_size - 1);
                                mask.insert(mask_key, mask_end);

                                // Skip if ungapped score is too low
                                if ungapped_score < MIN_UNGAPPED_SCORE {
                                    if diag_enabled {
                                        diag_inner.base.ungapped_low_score.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

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
                                    let context_key = (q_idx, q_f_idx);
                                    let search_space = *search_space_cache.entry(context_key).or_insert_with(|| {
                                        // Use NCBI-style length adjustment for single-sequence comparison
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-824
                                        SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, &params)
                                    });
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
                                let context_key = (q_idx, q_f_idx);
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
