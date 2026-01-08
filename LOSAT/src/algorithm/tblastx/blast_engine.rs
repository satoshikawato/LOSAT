//! TBLASTX search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:438-1497
//!
//! This module contains the main `run()` and `run_with_neighbor_map()` functions
//! that coordinate the TBLASTX search process:
//! - Query and subject sequence reading
//! - Lookup table construction
//! - Parallel subject scanning with two-hit extension
//! - HSP reevaluation, linking, and culling
//! - Output generation

use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::atomic::{AtomicI32, AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::mpsc::channel;
use std::time::{Duration, Instant};

use crate::common::{write_output_ncbi_order, Hit};
use crate::config::ScoringMatrix;
use crate::stats::{lookup_protein_params_ungapped, KarlinParams};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegMasker;

use super::args::TblastxArgs;
use super::chaining::UngappedHit;
use super::constants::{GAP_TRIGGER_BIT_SCORE, X_DROP_UNGAPPED_BITS};
use super::ncbi_cutoffs::{
    compute_eff_lengths_subject_mode_tblastx, cutoff_score_for_update_tblastx,
    cutoff_score_max_for_tblastx, gap_trigger_raw_score, x_drop_raw_score, BLAST_GAP_DECAY_RATE,
};
use super::diagnostics::{
    diagnostics_enabled, print_summary as print_diagnostics_summary, DiagnosticCounters,
};
use super::extension::{convert_coords, extend_hit_two_hit};
use super::lookup::{build_ncbi_lookup, encode_kmer, NeighborLookup};
use super::reevaluate::{
    get_num_identities_and_positives_ungapped, hsp_test,
    reevaluate_ungapped_hit_ncbi_translated,
};
use super::sum_stats_linking::{
    apply_sum_stats_even_gap_linking, compute_avg_query_length_ncbi, 
    find_smallest_lambda_params, LinkingParams,
};
use super::hsp_culling;
use super::translation::{generate_frames, QueryFrame};
use super::tracing::{
    trace_final_hit_if_match, trace_hsp_target, trace_match_target, trace_ungapped_hit_if_match,
};
use crate::stats::karlin::bit_score as calc_bit_score;

// Import OffsetPair from scan submodule
use super::scan::OffsetPair;

// Import DiagStruct from blast_extend module (NCBI blast_extend.c equivalent)
use super::blast_extend::DiagStruct;

// Import InitHSP and related functions from blast_gapalign module (NCBI blast_gapalign.c equivalent)
use super::blast_gapalign::{get_ungapped_hsp_list, trace_init_hsp_if_match, InitHSP};

// Import subject scanning functions from blast_aascan module (NCBI blast_aascan.c equivalent)
use super::blast_aascan::s_blast_aa_scan_subject;

struct WorkerState {
    tx: std::sync::mpsc::Sender<Vec<Hit>>,
    offset_pairs: Vec<OffsetPair>,
    diag_array: Vec<DiagStruct>,
    diag_offset: i32,
}

#[inline]
fn atomic_max_usize(dst: &AtomicUsize, val: usize) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
fn atomic_min_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val < cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
fn atomic_max_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

/// NCBI Blast_HSPListReevaluateUngapped equivalent
/// Reference: blast_hits.c:2609-2737, blast_engine.c:1492-1497
/// 
/// NCBI code flow:
/// 1. For each HSP in the list
/// 2. Get context and compute query_start (context start position)
/// 3. Call Blast_HSPReevaluateWithAmbiguitiesUngapped
/// 4. Call Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
/// 5. Delete HSPs that fail tests
/// 6. Sort by score (scores may have changed)
/// 
/// NCBI reference (verbatim, blast_hits.c:2694-2706):
/// ```c
/// context = hsp->context;
/// query_start = query_blk->sequence + query_info->contexts[context].query_offset;
/// if (kTranslateSubject)
///     subject_start = Blast_HSPGetTargetTranslation(target_t, hsp, NULL);
/// if (kNucleotideSubject) {
///     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesUngapped(hsp, query_start,
///         subject_start, word_params, sbp, kTranslateSubject);
/// }
/// ```
/// 
/// This function performs batch reevaluation on all HSPs in the list,
/// updating scores and trimming boundaries. HSPs that fail reevaluation
/// or subsequent tests are removed from the list.
fn reevaluate_ungapped_hsp_list(
    mut ungapped_hits: Vec<UngappedHit>,
    contexts: &[super::lookup::QueryContext],
    s_frames: &[super::translation::QueryFrame],
    cutoff_scores: &[Vec<i32>],
    args: &super::args::TblastxArgs,
    timing_enabled: bool,
    reeval_ns: &AtomicU64,
    reeval_calls: &AtomicU64,
    is_long_sequence: bool,
    stats_hsp_filtered_by_reeval: &mut usize,
) -> Vec<UngappedHit> {
    let mut kept_hits = Vec::new();
    
    for mut hit in ungapped_hits {
        let ctx = &contexts[hit.ctx_idx];
        let s_frame = &s_frames[hit.s_f_idx];
        let cutoff = cutoff_scores[hit.ctx_idx][hit.s_f_idx];
        
        // NCBI: query_start = query_blk->sequence + query_info->contexts[context].query_offset;
        // In LOSAT, ctx.aa_seq is the frame sequence (with sentinel at index 0)
        // query_start is the start of the context sequence (sentinel included)
        let query = &ctx.aa_seq;
        let subject = &s_frame.aa_seq;
        
        // NCBI: query = query_start + hsp->query.offset
        // hsp->query.offset is context-relative coordinate (sentinel included)
        // In LOSAT, q_aa_start is sentinel-exclusive, so we add 1 to get sentinel-inclusive
        let qs = (hit.q_aa_start + 1) as usize;  // +1 for sentinel
        let ss = (hit.s_aa_start + 1) as usize;  // +1 for sentinel
        let len_u = hit.q_aa_end.saturating_sub(hit.q_aa_start);
        
        if len_u == 0 {
            continue;
        }
        
        // NCBI: Blast_HSPReevaluateWithAmbiguitiesUngapped (blast_hits.c:2702-2705)
        // Reference: blast_hits.c:675-733
        // NCBI reference (verbatim, blast_hits.c:688):
        //   Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
        //   return s_UpdateReevaluatedHSPUngapped(hsp, cutoff_score, score, ...);
        // If score < cutoff_score, the HSP is deleted (s_UpdateReevaluatedHSPUngapped returns FALSE)
        let t0 = if timing_enabled { Some(std::time::Instant::now()) } else { None };
        let (new_qs, new_ss, new_len, new_score) = if let Some(result) =
            reevaluate_ungapped_hit_ncbi_translated(query, subject, qs, ss, len_u, cutoff)
        {
            result
        } else {
            // Deleted by NCBI reevaluation (score < cutoff_score)
            // NCBI reference: blast_hits.c:730-732 (s_UpdateReevaluatedHSPUngapped returns FALSE)
            if is_long_sequence {
                *stats_hsp_filtered_by_reeval += 1;
            }
            if let Some(t) = t0 {
                let elapsed = t.elapsed();
                reeval_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                reeval_calls.fetch_add(1, AtomicOrdering::Relaxed);
            }
            continue;
        };
        if let Some(t) = t0 {
            let elapsed = t.elapsed();
            reeval_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
            reeval_calls.fetch_add(1, AtomicOrdering::Relaxed);
        }
        
        // NCBI: Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
        // Reference: blast_hits.c:2708-2720
        let query_nomask = ctx.aa_seq_nomask.as_ref().unwrap_or(&ctx.aa_seq);
        let num_ident = if let Some((num_ident, _align_length)) =
            get_num_identities_and_positives_ungapped(
                query_nomask,
                subject,
                new_qs,
                new_ss,
                new_len,
            ) {
            num_ident
        } else {
            continue;
        };
        
        // NCBI: Blast_HSPTest (blast_hits.c:2719)
        let percent_identity = args.percent_identity;
        let min_hit_length = args.min_hit_length;
        if hsp_test(num_ident, new_len, percent_identity, min_hit_length) {
            continue;
        }
        
        // Update hit with reevaluated coordinates and score
        // Convert from sentinel-inclusive to sentinel-exclusive for reporting
        let new_qe = new_qs + new_len;
        let new_se = new_ss + new_len;
        let (qs_l, qe_l) = (new_qs.saturating_sub(1), new_qe.saturating_sub(1));
        let (ss_l, se_l) = (new_ss.saturating_sub(1), new_se.saturating_sub(1));
        
        hit.q_aa_start = qs_l;
        hit.q_aa_end = qe_l;
        hit.s_aa_start = ss_l;
        hit.s_aa_end = se_l;
        hit.raw_score = new_score;
        hit.num_ident = num_ident;
        
        trace_ungapped_hit_if_match("after_reevaluate", &hit);
        kept_hits.push(hit);
    }
    
    // NCBI: Sort the HSP array by score (scores may have changed!)
    // Reference: blast_hits.c:2734
    kept_hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));
    
    kept_hits
}

pub fn run(args: TblastxArgs) -> Result<()> {
    // Use neighbor_map mode for faster scanning with pre-computed neighbor relationships
    if args.neighbor_map {
        return run_with_neighbor_map(args);
    }
    
    // Optional timing breakdown (disabled by default to preserve output/parity logs)
    let timing_enabled = std::env::var_os("LOSAT_TIMING").is_some();
    let t_total = Instant::now();
    let mut t_read_queries = Duration::ZERO;
    let mut t_build_lookup = Duration::ZERO;
    let mut t_read_subjects = Duration::ZERO;

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
            anyhow::bail!("Invalid --only-qframe");
        }
    }
    if let Some(f) = only_sframe {
        if !valid_frame(f) {
            anyhow::bail!("Invalid --only-sframe");
        }
    }

    // [C] window = diag->window;
    let window = args.window_size as i32;
    // [C] wordsize = lookup->word_length;
    let wordsize: i32 = 3;

    // NCBI BLAST computes x_dropoff per-context using kbp[context]->Lambda:
    //   p->cutoffs[context].x_dropoff_init =
    //       (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // Reference: blast_parameters.c:219-221
    //
    // For tblastx, NCBI uses kbp_ideal->Lambda for all contexts (blast_stat.c:2796-2797):
    //   if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
    //      Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
    // So all contexts get the same x_dropoff, but we maintain per-context structure for parity.
    //
    // x_dropoff_per_context is populated after build_ncbi_lookup() creates the contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    let diag_enabled = diagnostics_enabled();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    let t_phase_read_queries = Instant::now();
    if args.verbose {
        eprintln!("Reading queries...");
    }
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();
    // NCBI tblastx low-complexity filtering uses SEG on translated protein sequences.
    // No nucleotide-level DUST masking is applied.
    //
    // NCBI reference (verbatim):
    //   else if (*ptr == 'L' || *ptr == 'T')
    //   { /* do low-complexity filtering; dust for blastn, otherwise seg.*/
    //       if (program_number == eBlastTypeBlastn
    //           || program_number == eBlastTypeMapping)
    //           SDustOptionsNew(&dustOptions);
    //       else
    //           SSegOptionsNew(&segOptions);
    //       ptr++;
    //   }
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:572-580

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

    if args.seg {
        let seg = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        for frames in &mut query_frames {
            for frame in frames {
                if frame.aa_seq.len() >= 3 {
                    for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }

                    // NCBI BLAST query masking semantics (SEG, etc.): keep an unmasked copy
                    // (`sequence_nomask`), then overwrite masked residues in the working
                    // query sequence buffer with X.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1 kProtMask = 21;     /* X in NCBISTDAA */
                    //   query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
                    //   query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
                    //   buffer[index] = kMaskingLetter;
                    //
                    // LOSAT uses NCBISTDAA encoding where X = 21.
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        const X_MASK_NCBISTDAA: u8 = 21; // NCBI: kProtMask = 21
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1); // keep last sentinel untouched
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
    }

    t_read_queries = t_phase_read_queries.elapsed();

    let t_phase_build_lookup = Instant::now();
    if args.verbose {
        eprintln!("Building lookup table...");
    }
    // Note: karlin_params argument is now unused - computed per context in build_ncbi_lookup()
    // We still pass it for x_dropoff calculation (which uses ideal params for all contexts in tblastx)
    let (lookup, contexts) = build_ncbi_lookup(
        &query_frames,
        args.threshold,
        args.include_stop_seeds,
        args.ncbi_stop_stop_score,
        false, // lazy_neighbors disabled - use neighbor_map mode instead
        &ungapped_params_for_xdrop, // Used for x_dropoff calculation only
    );
    t_build_lookup = t_phase_build_lookup.elapsed();

    // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
    // Compute per-context x_dropoff using kbp[context]->Lambda.
    // Reference: blast_parameters.c:219-221
    //
    // For tblastx, all contexts use kbp_ideal (BLOSUM62 ungapped Lambda = 0.3176),
    // so x_dropoff = ceil(7.0 * LN2 / 0.3176) = 16 for all contexts.
    // We maintain per-context structure for NCBI parity.
    //
    // NCBI BLAST dynamic x_dropoff (blast_parameters.c:380-383):
    //   if (curr_cutoffs->x_dropoff_init == 0)
    //      curr_cutoffs->x_dropoff = new_cutoff;  // x_dropoff = cutoff_score
    //   else
    //      curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
    //
    // For TBLASTX, x_dropoff_init = 16 (not 0), so this dynamic update is not triggered.
    // However, we store x_dropoff_init here and apply the logic during subject processing
    // where cutoff_score is available.
    let x_dropoff_init = x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ungapped_params_for_xdrop, 1.0);
    let x_dropoff_per_context: Vec<i32> = vec![x_dropoff_init; contexts.len()];

    let t_phase_read_subjects = Instant::now();
    if args.verbose {
        eprintln!("Reading subjects...");
    }
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }
    t_read_subjects = t_phase_read_subjects.elapsed();

    // NCBI BLAST Karlin params for TBLASTX (ungapped-only algorithm):
    // 
    // TBLASTX is explicitly ungapped-only (blast_options.c line 869-873):
    //   "Gapped search is not allowed for tblastx"
    //
    // For bit score and E-value calculation, NCBI uses sbp->kbp (ungapped):
    //   blast_hits.c line 1833: kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    //   blast_hits.c line 1918: same pattern in Blast_HSPListGetBitScores
    //
    // NCBI reference (blast_setup.c:768):
    //   kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
    // For tblastx, gapped_calculation = FALSE, so kbp_ptr = sbp->kbp (ungapped params)
    // Therefore, ALL calculations (eff_searchsp, cutoff, bit score, E-value) use UNGAPPED params
    //
    // BLOSUM62 ungapped: lambda=0.3176, K=0.134 (used for tblastx)
    // BLOSUM62 gapped:   lambda=0.267,  K=0.041 (NOT used for tblastx)
    let ungapped_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    // Note: gapped_params is kept for API compatibility but NOT used for tblastx
    let gapped_params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    // Use UNGAPPED params for all calculations (NCBI parity for tblastx)
    let params = ungapped_params.clone();

    // Compute NCBI-style average query length for linking cutoffs
    // Reference: blast_parameters.c:CalculateLinkHSPCutoffs (lines 1023-1026)
    // NCBI uses average over ALL contexts including zero-length (frame restriction via strand)
    let query_nucl_lengths: Vec<usize> = queries_raw.iter().map(|r| r.seq().len()).collect();
    let avg_query_length = compute_avg_query_length_ncbi(&query_nucl_lengths);

    if args.verbose {
        eprintln!(
            "Searching {} queries vs {} subjects... (avg_query_length={})",
            queries_raw.len(),
            subjects_raw.len(),
            avg_query_length
        );
    }

    let t_search_start = Instant::now();
    let scan_ns = AtomicU64::new(0);
    let scan_calls = AtomicU64::new(0);
    let ungapped_ns = AtomicU64::new(0);
    let ungapped_calls = AtomicU64::new(0);
    let reeval_ns = AtomicU64::new(0);
    let reeval_calls = AtomicU64::new(0);
    let linking_ns = AtomicU64::new(0);
    let linking_calls = AtomicU64::new(0);
    let identity_ns = AtomicU64::new(0);
    let identity_calls = AtomicU64::new(0);

    let bar = if args.verbose {
        let bar = ProgressBar::new(subjects_raw.len() as u64);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
                .unwrap(),
        );
        bar
    } else {
        ProgressBar::hidden()
    };

    let (tx, rx) = channel::<Vec<Hit>>();
    let out_path = args.out.clone();
    let evalue_threshold = args.evalue;

    let writer = std::thread::spawn(move || -> Result<()> {
        let mut all: Vec<Hit> = Vec::new();
        while let Ok(h) = rx.recv() {
            all.extend(h);
        }
        all.retain(|h| h.e_value <= evalue_threshold);
        // NCBI-style output ordering: query (input order) → subject (best_evalue/score/oid) → HSP (score/coords)
        // Reference: BLAST_LinkHsps() + s_EvalueCompareHSPLists() + ScoreCompareHSPs()
        write_output_ncbi_order(all, out_path.as_ref())?;
        Ok(())
    });

    // Diagonal array sizing MUST match NCBI's `s_BlastDiagTableNew`:
    // it depends only on (query_length + window_size), not on subject length.
    //
    // NCBI reference (verbatim):
    //   diag_array_length = 1;
    //   while (diag_array_length < (qlen+window_size))
    //       diag_array_length = diag_array_length << 1;
    //   diag_table->diag_array_length = diag_array_length;
    //   diag_table->diag_mask = diag_array_length-1;
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    //
    // `query_length` here is the concatenated translated query buffer length
    // (frames share a boundary sentinel). In LOSAT this equals the end offset
    // of the last query context buffer.
    // NCBI BLAST diag array sizing:
    // diag_array_length = next_power_of_2(qlen + window_size)
    // For tblastx, qlen = total concatenated query buffer (all 6 frames)
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    let query_length: i32 = contexts
        .last()
        .map(|c| c.frame_base + (c.aa_seq.len() as i32 - 1))
        .unwrap_or(0);
    let mut diag_array_size: i32 = 1;
    while diag_array_size < (query_length + window) {
        diag_array_size <<= 1;
    }
    let diag_mask: i32 = diag_array_size - 1;

    // [C] array_size for offset_pairs
    // NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain
    // Reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h + lookup_wrap.c
    const OFFSET_ARRAY_SIZE: i32 = 4096;
    let offset_array_size: i32 = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);

    let lookup_ref = &lookup;
    let contexts_ref = &contexts;
    let query_ids_ref = &query_ids;
    let ungapped_params_ref = &ungapped_params;
    let _gapped_params_ref = &gapped_params;  // Unused - tblastx uses ungapped params

    subjects_raw.par_iter().enumerate().for_each_init(
        || WorkerState {
            tx: tx.clone(),
            offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
            diag_array: vec![DiagStruct::default(); diag_array_size as usize],
            // NCBI: diag_table->offset = window_size;
            // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
            diag_offset: window,
        },
        |st, (s_idx, s_rec)| {
            // NCBI: Creates ewp (diagonal table) ONCE per SUBJECT via BlastExtendWordNew
            // (blast_engine.c:1002). Each subject sequence gets a fresh diagonal array.
            // Reference: blast_extend.c:109-180 (BlastExtendWordNew) allocates with calloc,
            // which zeros all entries. The diag_offset is initialized to window.
            //
            // Reset diagonal state for each subject to match NCBI behavior:
            for d in st.diag_array.iter_mut() {
                *d = DiagStruct::default();
            }
            st.diag_offset = window;

            let mut s_frames = generate_frames(s_rec.seq(), &db_code);
            if let Some(f) = only_sframe {
                s_frames.retain(|x| x.frame == f);
            }

            let s_id = s_rec
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            let s_len = s_rec.seq().len();

            // [C] BlastOffsetPair *offset_pairs
            let offset_pairs = &mut st.offset_pairs;

            // [C] DiagStruct *diag_array = diag->hit_level_array;
            let diag_array = &mut st.diag_array;

            // [C] diag_offset = diag->offset;  (reset to window per-subject)
            let mut diag_offset: i32 = st.diag_offset;

            // Precompute per-context cutoff scores using NCBI BLAST algorithm.
            // Reference: blast_parameters.c BlastInitialWordParametersUpdate
            //
            // NCBI cutoff calculation for tblastx ungapped path:
            // 1. gap_trigger from ungapped params (kbp_std)
            // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
            // 3. Per-subject update: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
            // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
            //
            // For -subject mode, subject_len_nucl is used (not per-frame AA length).
            let subject_len_nucl = s_len as i64;
            let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; contexts_ref.len()];
            // NCBI word_params->cutoff_score_min = min of cutoffs across all contexts
            // Reference: blast_parameters.c BlastInitialWordParametersUpdate line 401-403
            let mut cutoff_score_min = i32::MAX;
            
            // =======================================================================
            // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
            // =======================================================================
            // NCBI stores these in query_info->contexts[ctx].length_adjustment and
            // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths (blast_setup.c:700-850).
            // We precompute them here and pass to sum_stats_linking for NCBI parity.
            // Reference: blast_setup.c:846-847
            //   query_info->contexts[index].eff_searchsp = effective_search_space;
            //   query_info->contexts[index].length_adjustment = length_adjustment;
            let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
            let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
            
            // Precompute gap_trigger once (same for all contexts with same params)
            // Reference: blast_parameters.c:343-344
            let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ungapped_params_ref);
            
            for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
                // NCBI: query_length = query_info->contexts[context].query_length
                let query_len_aa = ctx.aa_len as i64;
                
                // =======================================================================
                // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                // length_adjustment and eff_searchsp from a single source of truth.
                // This mirrors BLAST_CalcEffLengths which computes and stores both values.
                // Reference: blast_setup.c:821-847
                // =======================================================================
                let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    ungapped_params_ref,  // tblastx uses ungapped params (kbp_gap is NULL)
                );
                let eff_searchsp = eff_lengths.eff_searchsp;
                length_adj_per_context.push(eff_lengths.length_adjustment);
                eff_searchsp_per_context.push(eff_searchsp);
                
                // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
                // This uses the effective search space WITH length adjustment
                // Reference: blast_parameters.c:942-946
                let cutoff_score_max = cutoff_score_max_for_tblastx(
                    eff_searchsp,
                    evalue_threshold,  // User's E-value (typically 10.0)
                    ungapped_params_ref,
                );
                
                // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
                // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
                // Reference: blast_parameters.c:348-374
                let cutoff = cutoff_score_for_update_tblastx(
                    query_len_aa,
                    subject_len_nucl,  // NUCLEOTIDE length, NOT divided by 3!
                    gap_trigger,
                    cutoff_score_max,
                    BLAST_GAP_DECAY_RATE,  // 0.5
                    ungapped_params_ref,
                    1.0,  // scale_factor (standard BLOSUM62)
                );
                
                // DEBUG: Print cutoff values for first context
                static PRINTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
                if diag_enabled
                    && ctx_idx == 0
                    && !PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed)
                {
                    eprintln!("[DEBUG CUTOFF] query_len_aa={}, subject_len_nucl={}", query_len_aa, subject_len_nucl);
                    eprintln!("[DEBUG CUTOFF] eff_searchsp={}", eff_searchsp);
                    eprintln!("[DEBUG CUTOFF] length_adjustment={}", eff_lengths.length_adjustment);
                    eprintln!("[DEBUG CUTOFF] cutoff_score_max={}", cutoff_score_max);
                    eprintln!("[DEBUG CUTOFF] gap_trigger={}", gap_trigger);
                    eprintln!("[DEBUG CUTOFF] final cutoff={}", cutoff);
                }
                
                // Track minimum cutoff for linking
                cutoff_score_min = cutoff_score_min.min(cutoff);
                
                // All subject frames use the same cutoff (NCBI behavior)
                for sf_idx in 0..s_frames.len() {
                    cutoff_scores[ctx_idx][sf_idx] = cutoff;
                }
            }
            // If no contexts, use 0 as fallback
            if cutoff_score_min == i32::MAX {
                cutoff_score_min = 0;
            }

            // NCBI: Each subject frame gets its own init_hitlist that is reset between frames.
            // Reference: blast_engine.c:491 BlastInitHitListReset(init_hitlist)
            // After each frame, init_hsps are converted and merged into combined_ungapped_hits.

            // Statistics for HSP saving analysis (long sequences only)
            let is_long_sequence = subject_len_nucl > 600_000;
            let mut stats_hsp_saved = 0usize;
            let mut stats_hsp_filtered_by_cutoff = 0usize;
            let mut stats_hsp_filtered_by_reeval = 0usize;
            let mut stats_hsp_filtered_by_hsp_test = 0usize;
            let mut stats_score_distribution: Vec<i32> = Vec::new();

            // NCBI: Combined HSP list across all subject frames
            // Reference: blast_engine.c:438 BlastHSPList* combined_hsp_list
            let mut combined_ungapped_hits: Vec<UngappedHit> = Vec::new();

            for (s_f_idx, s_frame) in s_frames.iter().enumerate() {
                // NCBI: BlastInitHitListReset(init_hitlist) - reset per-frame
                // Reference: blast_engine.c:491
                let mut init_hsps: Vec<InitHSP> = Vec::new();

                // NCBI: Diagonal state is NOT reset between subject frames.
                // NCBI shares ewp (diag_table) across all 6 subject frame iterations.
                // Reference: blast_engine.c:805-855
                // The per-subject reset (above, line ~1646) matches NCBI's Blast_ExtendWordNew calloc.
                // Blast_ExtendWordExit at the end of this loop increments diag_offset.

                let subject = &s_frame.aa_seq;
                let s_aa_len = s_frame.aa_len;
                if subject.len() < 5 {
                    continue;
                }

                // [C] scan_range[1] = subject->seq_ranges[0].left;
                // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                let mut scan_range: [i32; 3] = [0, 1, (subject.len() - 4) as i32];

                // [C] while (scan_range[1] <= scan_range[2])
                while scan_range[1] <= scan_range[2] {
                    let prev_scan_left = scan_range[1];
                    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                    let t0 = if timing_enabled { Some(Instant::now()) } else { None };
                    let hits = s_blast_aa_scan_subject(
                        lookup_ref,
                        subject,
                        offset_pairs,
                        offset_array_size,
                        &mut scan_range,
                    );
                    if let Some(t0) = t0 {
                        scan_ns.fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
                        scan_calls.fetch_add(1, AtomicOrdering::Relaxed);
                    }

                    if diag_enabled && hits > 0 {
                        diagnostics
                            .base
                            .kmer_matches
                            .fetch_add(hits as usize, AtomicOrdering::Relaxed);

                        // DEBUG: Check for duplicate offset pairs in scan output
                        if is_long_sequence {
                            let mut seen: HashSet<(i32, i32)> = HashSet::with_capacity(hits as usize);
                            let mut duplicate_count = 0usize;
                            for i in 0..hits as usize {
                                let pair = unsafe { &*offset_pairs.as_ptr().add(i) };
                                if !seen.insert((pair.q_off, pair.s_off)) {
                                    duplicate_count += 1;
                                }
                            }
                            if duplicate_count > 0 {
                                eprintln!("[DEBUG SCAN_DUPES] s_f_idx={} scan_range=[{},{}] hits={} duplicates={} ({:.2}%)",
                                    s_f_idx, prev_scan_left, scan_range[1], hits, duplicate_count,
                                    (duplicate_count as f64 / hits as f64) * 100.0);
                            }
                        }
                    }

                    if hits == 0 && scan_range[1] == prev_scan_left {
                        // Safety guard: with correct NCBI-sized offset arrays, this should not happen.
                        // If it does, breaking avoids an infinite loop.
                        break;
                    }

                    // [C] for (i = 0; i < hits; ++i)
                    // OPTIMIZATION: Use raw pointers to eliminate bounds checking in hot loop
                    let diag_ptr = diag_array.as_mut_ptr();
                    let offset_pairs_ptr = offset_pairs.as_ptr();
                    
                    for i in 0..hits as usize {
                        // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                        let pair = unsafe { &*offset_pairs_ptr.add(i) };
                        let query_offset = pair.q_off;
                        let subject_offset = pair.s_off;

                        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                        let diag_coord = ((query_offset - subject_offset) & diag_mask) as usize;
                        
                        // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
                        let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

                        // [C] if (diag_array[diag_coord].flag)
                        if diag_entry.flag != 0 {
                            // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                            if subject_offset + diag_offset < diag_entry.last_hit {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_masked
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                continue;
                            }
                            // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            // [C] diag_array[diag_coord].flag = 0;
                            diag_entry.last_hit = subject_offset + diag_offset;
                            diag_entry.flag = 0;
                        }
                        // [C] else
                        else {
                            // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                            let last_hit = diag_entry.last_hit - diag_offset;
                            // [C] diff = subject_offset - last_hit;
                            let diff = subject_offset - last_hit;

                            // [C] if (diff >= window)
                            if diff >= window {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_too_far
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }

                            // [C] if (diff < wordsize)
                            if diff < wordsize {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_overlap
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                continue;
                            }

                            // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                            let ctx_idx = lookup_ref.get_context_idx(query_offset);
                            let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
                            let q_raw = (query_offset - ctx.frame_base) as usize;
                            // NCBI uses masked sequence for extension (same as neighbor-map mode)
                            let query = &ctx.aa_seq;

                            // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                            if query_offset - diff < ctx.frame_base {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_ctx_boundary_fail
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }

                            if diag_enabled {
                                diagnostics
                                    .base
                                    .seeds_second_hit_window
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                diagnostics
                                    .base
                                    .seeds_passed
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                            }

                            // [C] cutoffs = word_params->cutoffs + curr_context;
                            let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx).get_unchecked(s_f_idx) };
                            // [C] cutoffs->x_dropoff (per-context x_dropoff)
                            // Reference: aa_ungapped.c:579
                            let x_dropoff = unsafe { *x_dropoff_per_context.get_unchecked(ctx_idx) };

                            // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
                            //                                   last_hit + wordsize, subject_offset, query_offset, ...)
                            // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                            let t0 = if timing_enabled { Some(Instant::now()) } else { None };
                            let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
                                extend_hit_two_hit(
                                    query,
                                    subject,
                                    (last_hit + wordsize) as usize, // s_left_off (end of first hit word)
                                    subject_offset as usize, // s_right_off (second hit start)
                                    q_raw as usize,          // q_right_off (second hit start, local)
                                    x_dropoff,               // x_dropoff (per-context, NCBI parity)
                                );
                            if let Some(t0) = t0 {
                                ungapped_ns.fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
                                ungapped_calls.fetch_add(1, AtomicOrdering::Relaxed);
                            }

                            let hsp_q: i32 = hsp_q_u as i32;
                            let hsp_s: i32 = hsp_s_u as i32;
                            let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
                            let s_last_off: i32 = s_last_off_u as i32;

                            if diag_enabled {
                                diagnostics
                                    .base
                                    .ungapped_extensions
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                if right_extend {
                                    diagnostics
                                        .base
                                        .ungapped_two_hit_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                } else {
                                    diagnostics
                                        .base
                                        .ungapped_one_hit_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                if hsp_len > 0 {
                                    diagnostics
                                        .base
                                        .extension_total_length
                                        .fetch_add(hsp_len as usize, AtomicOrdering::Relaxed);
                                    atomic_max_usize(
                                        &diagnostics.base.extension_max_length,
                                        hsp_len as usize,
                                    );
                                }
                            }

                            // NCBI: Update diagonal state based on right_extend
                            // Reference: aa_ungapped.c:596-606
                            // if (right_extend) {
                            //     diag_array[diag_coord].flag = 1;
                            //     diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                            // } else {
                            //     diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            // }
                            if right_extend {
                                // "If an extension to the right happened, reset the last hit
                                //  so that future hits to this diagonal must start over."
                                diag_entry.flag = 1;
                                diag_entry.last_hit =
                                    s_last_off - (wordsize - 1) + diag_offset;
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .mask_updates
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                            } else {
                                // "Otherwise, make the present hit into the previous hit for this diagonal"
                                diag_entry.last_hit = subject_offset + diag_offset;
                            }

                            // [C] if (score >= cutoffs->cutoff_score)
                            // NCBI reference: aa_ungapped.c:575-591 (Extension後のcutoffチェック)
                            if is_long_sequence {
                                if score >= cutoff {
                                    stats_score_distribution.push(score);
                                    stats_hsp_saved += 1;  // Count saved HSPs
                                } else {
                                    stats_hsp_filtered_by_cutoff += 1;
                                }
                            }
                            if score >= cutoff {
                                if diag_enabled {
                                    diagnostics
                                        .ungapped_only_hits
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                // Extra debug for a traced HSP: print seed/extension inputs and cutoffs.
                                if let Some(target) = trace_hsp_target() {
                                    // Compute outfmt coords for this candidate init-hsp (same logic as trace_init_hsp_if_match).
                                    let q_aa_start = (hsp_q_u as usize).saturating_sub(1);
                                    let q_aa_end = (hsp_qe_u as usize).saturating_sub(1);
                                    let s_aa_start = (hsp_s_u as usize).saturating_sub(1);
                                    let s_aa_end = (_hsp_se_u as usize).saturating_sub(1);
                                    let (q_start_dna, q_end_dna) =
                                        convert_coords(q_aa_start, q_aa_end, ctx.frame, ctx.orig_len);
                                    let (s_start_dna, s_end_dna) =
                                        convert_coords(s_aa_start, s_aa_end, s_frame.frame, s_len);
                                    if trace_match_target(target, q_start_dna, q_end_dna, s_start_dna, s_end_dna) {
                                        eprintln!(
                                            "[TRACE_HSP] seed/extend ctx_idx={} s_f_idx={} q_frame={} s_frame={} score={} cutoff={} x_dropoff={} last_hit={} subject_offset={} diff={} q_raw={} query_offset={} diag_coord={} right_extend={} s_last_off={}",
                                            ctx_idx,
                                            s_f_idx,
                                            ctx.frame,
                                            s_frame.frame,
                                            score,
                                            cutoff,
                                            x_dropoff,
                                            last_hit,
                                            subject_offset,
                                            diff,
                                            q_raw,
                                            query_offset,
                                            diag_coord,
                                            right_extend,
                                            s_last_off,
                                        );
                                    }
                                }
                                // NCBI: BlastSaveInitHsp equivalent
                                // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                                // Store HSP with absolute coordinates (before coordinate conversion)
                                //
                                // hsp_q is frame-relative coordinate (with sentinel, array index)
                                // frame_base is concatenated buffer start position (with sentinel, array index)
                                // Calculate absolute coordinate in concatenated buffer
                                // NCBI: ungapped_data->q_start = q_start (absolute, sentinel included)
                                let hsp_q_absolute = ctx.frame_base + (hsp_q as i32);
                                let hsp_qe_absolute = ctx.frame_base + ((hsp_q + hsp_len) as i32);
                                
                                let init = InitHSP {
                                    q_start_absolute: hsp_q_absolute,
                                    q_end_absolute: hsp_qe_absolute,
                                    s_start: hsp_s,
                                    s_end: hsp_s + hsp_len,
                                    score,
                                    ctx_idx,
                                    s_f_idx,
                                    q_idx: ctx.q_idx,
                                    s_idx: s_idx as u32,
                                    q_frame: ctx.frame,
                                    s_frame: s_frame.frame,
                                    q_orig_len: ctx.orig_len,
                                    s_orig_len: s_len,
                                };
                                trace_init_hsp_if_match("init_hsp_saved", &init, contexts_ref);
                                init_hsps.push(init);
                            } else if diag_enabled {
                                diagnostics
                                    .ungapped_cutoff_failed
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                atomic_min_i32(&diagnostics.ungapped_cutoff_failed_min_score, score);
                                atomic_max_i32(&diagnostics.ungapped_cutoff_failed_max_score, score);
                            }
                        }
                    }
                }

                // [C] Blast_ExtendWordExit(ewp, subject->length);
                //
                // NCBI reference (verbatim):
                //   if (ewp->diag_table->offset >= INT4_MAX / 4) {
                //       ewp->diag_table->offset = ewp->diag_table->window;
                //       s_BlastDiagClear(ewp->diag_table);
                //   } else {
                //       ewp->diag_table->offset += subject_length + ewp->diag_table->window;
                //   }
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:167-173
                if diag_offset >= i32::MAX / 4 {
                    diag_offset = window;
                    // s_BlastDiagClear(): clear all diagonal state when offset risks overflow.
                    // NCBI sets last_hit = -window (blast_extend.c:103)
                    for d in diag_array.iter_mut() {
                        *d = DiagStruct::clear(window);
                    }
                } else {
                    diag_offset += s_aa_len as i32 + window;
                }

                // NCBI: BLAST_GetUngappedHSPList equivalent - per-frame conversion
                // Reference: blast_engine.c:561-562
                // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
                //         hit_params->options, &hsp_list);
                if !init_hsps.is_empty() {
                    let frame_ungapped_hits = get_ungapped_hsp_list(
                        init_hsps,
                        contexts_ref,
                        &s_frames,
                    );

                    // NCBI: Blast_HSPListsMerge - merge per-frame results
                    // Reference: blast_engine.c:581-584
                    combined_ungapped_hits.extend(frame_ungapped_hits);
                }
            } // End of subject frame loop

            // NCBI: Blast_HSPListReevaluateUngapped equivalent
            // Reference: blast_engine.c:1492-1497, blast_hits.c:2609-2737
            // Perform batch reevaluation on all HSPs after merging all frames
            let ungapped_hits = reevaluate_ungapped_hsp_list(
                combined_ungapped_hits,
                contexts_ref,
                &s_frames,
                &cutoff_scores,
                &args,
                timing_enabled,
                &reeval_ns,
                &reeval_calls,
                is_long_sequence,
                &mut stats_hsp_filtered_by_reeval,
            );

            if !ungapped_hits.is_empty() {
                    // DEBUG: Print HSP statistics for long sequences
                    if is_long_sequence {
                        eprintln!("[DEBUG HSP_STATS] After reevaluation: {} HSPs", ungapped_hits.len());
                        eprintln!("[DEBUG HSP_STATS] Saved by cutoff: {}", stats_hsp_saved);
                        eprintln!("[DEBUG HSP_STATS] Filtered by cutoff: {}", stats_hsp_filtered_by_cutoff);
                        eprintln!("[DEBUG HSP_STATS] Filtered by reeval: {}", stats_hsp_filtered_by_reeval);
                        if !stats_score_distribution.is_empty() {
                            let min_score = stats_score_distribution.iter().min().unwrap();
                            let max_score = stats_score_distribution.iter().max().unwrap();
                            let avg_score = stats_score_distribution.iter().sum::<i32>() as f64 / stats_score_distribution.len() as f64;
                            eprintln!("[DEBUG HSP_STATS] Score range: {} - {} (avg: {:.2})", 
                                      min_score, max_score, avg_score);
                        }
                    }
                    if diag_enabled {
                        diagnostics
                            .base
                            .hsps_before_chain
                            .fetch_add(ungapped_hits.len(), AtomicOrdering::Relaxed);
                    }
                    // NCBI CalculateLinkHSPCutoffs parameters
                    let linking_params = LinkingParams {
                        avg_query_length,
                        subject_len_nucl,
                        cutoff_score_min,
                        scale_factor: 1.0, // Standard BLOSUM62
                    gap_decay_rate: 0.5, // BLAST_GAP_DECAY_RATE
                };
                // Build NCBI-style subject frame base offsets for sum-statistics linking.
                // In NCBI, HSP coords live in a concatenated translation buffer with sentinels.
                // LOSAT uses per-frame sequences; for linking we emulate absolute offsets by
                // concatenating frames in the same order as `generate_frames()`.
                let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(s_frames.len());
                let mut base: i32 = 0;
                for f in &s_frames {
                    subject_frame_bases.push(base);
                    // NCBI concatenation shares the trailing NULLB sentinel between frames:
                    //   offset += length + 1;
                    // where `length` is the number of residues (excluding sentinels).
                    // Source: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
                    base += f.aa_seq.len() as i32 - 1;
                }

                // NCBI parity: s_BlastFindSmallestLambda selects Karlin params with smallest lambda
                // Reference: blast_parameters.c:92-112, 1012-1013
                // For tblastx, all contexts use kbp_ideal (same lambda), so this is equivalent
                // to using ungapped_params directly. We maintain this structure for parity.
                let context_params: Vec<KarlinParams> = contexts_ref.iter().map(|_| params.clone()).collect();
                let linking_params_for_cutoff = find_smallest_lambda_params(&context_params)
                    .unwrap_or_else(|| params.clone());

                // Output HSP saving statistics for long sequences
                if diag_enabled
                    && is_long_sequence
                    && (stats_hsp_saved > 0
                        || stats_hsp_filtered_by_cutoff > 0
                        || stats_hsp_filtered_by_reeval > 0
                        || stats_hsp_filtered_by_hsp_test > 0)
                {
                    let total_attempted = stats_hsp_saved + stats_hsp_filtered_by_cutoff + stats_hsp_filtered_by_reeval + stats_hsp_filtered_by_hsp_test;
                    eprintln!("[DEBUG HSP_SAVING] subject_len_nucl={}, cutoff={}", subject_len_nucl, cutoff_score_min);
                    eprintln!("[DEBUG HSP_SAVING] total_attempted={}, saved={}, filtered_by_cutoff={}, filtered_by_reeval={}, filtered_by_hsp_test={}", 
                        total_attempted, stats_hsp_saved, stats_hsp_filtered_by_cutoff, stats_hsp_filtered_by_reeval, stats_hsp_filtered_by_hsp_test);
                    if !stats_score_distribution.is_empty() {
                        stats_score_distribution.sort();
                        let min_score = stats_score_distribution[0];
                        let max_score = stats_score_distribution[stats_score_distribution.len() - 1];
                        let median_score = stats_score_distribution[stats_score_distribution.len() / 2];
                        let low_score_count = stats_score_distribution.iter().filter(|&&s| s < 30).count();
                        eprintln!("[DEBUG HSP_SAVING] score_distribution: min={}, max={}, median={}, low_score(<30)={}/{} ({:.2}%)", 
                            min_score, max_score, median_score, low_score_count, stats_score_distribution.len(),
                            if stats_score_distribution.len() > 0 { (low_score_count as f64 / stats_score_distribution.len() as f64) * 100.0 } else { 0.0 });
                    }
                }
                
                // NCBI Parity: Pass pre-computed length_adjustment and eff_searchsp per context
                // These values are stored in query_info->contexts[ctx] in NCBI and referenced
                // by link_hsps.c for BLAST_SmallGapSumE/BLAST_LargeGapSumE calculations.
                let t_linking = if timing_enabled { Some(Instant::now()) } else { None };
                let linked = apply_sum_stats_even_gap_linking(
                    ungapped_hits,
                    &linking_params_for_cutoff,
                    &linking_params,
                    contexts_ref,
                    &subject_frame_bases,
                    &length_adj_per_context,
                    &eff_searchsp_per_context,
                );
                if trace_hsp_target().is_some() {
                    for h in &linked {
                        trace_ungapped_hit_if_match("after_linking", h);
                    }
                }
                if let Some(t) = t_linking {
                    let elapsed = t.elapsed();
                    linking_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                    linking_calls.fetch_add(1, AtomicOrdering::Relaxed);
                }
                if diag_enabled {
                    diagnostics
                        .base
                        .hsps_after_chain
                        .fetch_add(linked.len(), AtomicOrdering::Relaxed);
                }
                
                // NCBI HSP culling (conditional on culling_limit > 0)
                // Reference: hspfilter_culling.c - applied after linking, before output
                // Default for tblastx: culling_limit = 0 (disabled)
                let linked_before_cull = linked.len();
                let linked = if args.culling_limit > 0 {
                    hsp_culling::apply_culling(
                        linked,
                        contexts_ref,
                        args.culling_limit,
                    )
                } else {
                    // Default: no culling (matches NCBI tblastx default)
                    linked
                };
                if diag_enabled && args.culling_limit > 0 {
                    diagnostics
                        .base
                        .hsps_culled_dominated
                        .fetch_add(linked_before_cull.saturating_sub(linked.len()), AtomicOrdering::Relaxed);
                }
                
                // DEBUG: Collect statistics before filtering
                let total_linked = linked.len();
                let mut stats_single_hsps = 0usize;
                let mut stats_chain_heads = 0usize;
                let mut stats_chain_members = 0usize;
                for h in &linked {
                    if h.linked_set && !h.start_of_chain {
                        stats_chain_members += 1;
                    } else if !h.linked_set {
                        stats_single_hsps += 1;
                    } else if h.start_of_chain {
                        stats_chain_heads += 1;
                    }
                }
                
                let mut final_hits: Vec<Hit> = Vec::new();
                let mut filtered_by_evalue = 0usize;
                for h in linked {
                    // NCBI reference (verbatim, link_hsps.c:1018-1020):
                    //   /* If this is not a single piece or the start of a chain, then Skip it. */
                    //   if (H->linked_set == TRUE && H->start_of_chain == FALSE)
                    //       continue;
                    // NOTE: This "continue" in NCBI is NOT a filter - NCBI then walks the link
                    // pointer from each chain head to include all chain members (lines 1047-1059).
                    // LOSAT's linking already produces a flat list of ALL HSPs (including chain
                    // members) so we output everything without this skip logic.

                    // NCBI reference: E-value filtering is applied during output conversion
                    // The exact timing may differ, but the threshold check is standard
                    if h.e_value > evalue_threshold {
                        filtered_by_evalue += 1;
                        if diag_enabled {
                            diagnostics
                                .ungapped_evalue_failed
                                .fetch_add(1, AtomicOrdering::Relaxed);
                        }
                        continue;
                    }
                    
                    if diag_enabled {
                        diagnostics
                            .ungapped_evalue_passed
                            .fetch_add(1, AtomicOrdering::Relaxed);
                    }

                    let ctx = &contexts_ref[h.ctx_idx];
                    let s_frame = &s_frames[h.s_f_idx];

                    // Compute identity/mismatch on demand (final hits only)
                    // NCBI computes identities using the *unmasked* query sequence buffer
                    // (`sequence_nomask`), even when the working query sequence was masked
                    // (e.g., SEG) to influence seeding/extension.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1* query_nomask = query_blk->sequence_nomask +
                    //       query_info->contexts[context].query_offset;
                    let len = h.q_aa_end.saturating_sub(h.q_aa_start);
                    let q0 = h.q_aa_start + 1; // +1 for sentinel
                    let s0 = h.s_aa_start + 1; // +1 for sentinel
                    let q_seq_nomask: &[u8] = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
                    // Use SIMD-optimized identity calculation (already implemented in reevaluate.rs)
                    let t_identity = if timing_enabled { Some(Instant::now()) } else { None };
                    let matches = get_num_identities_and_positives_ungapped(
                        q_seq_nomask,
                        &s_frame.aa_seq,
                        q0,
                        s0,
                        len,
                    ).map(|(num_ident, _)| num_ident).unwrap_or(0);
                    if let Some(t) = t_identity {
                        let elapsed = t.elapsed();
                        identity_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                        identity_calls.fetch_add(1, AtomicOrdering::Relaxed);
                    }
                    let mismatch = len.saturating_sub(matches);
                    let identity = if len > 0 {
                        (matches as f64 / len as f64) * 100.0
                    } else {
                        0.0
                    };

                    let bit = calc_bit_score(h.raw_score, &params);
                    let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, ctx.frame, ctx.orig_len);
                    let (s_start, s_end) =
                        convert_coords(h.s_aa_start, h.s_aa_end, s_frame.frame, s_len);

                    let out_hit = Hit {
                        query_id: query_ids_ref[ctx.q_idx as usize].clone(),
                        subject_id: s_id.clone(),
                        identity,
                        length: len,
                        mismatch,
                        gapopen: 0,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        e_value: h.e_value,
                        bit_score: bit,
                        q_idx: ctx.q_idx,
                        s_idx: h.s_idx,
                        raw_score: h.raw_score,
                    };
                    trace_final_hit_if_match("output_hit", &out_hit);
                    final_hits.push(out_hit);
                }
                
                // DEBUG: Print output filtering statistics
                // NCBI reference: link_hsps.c:1018-1020 - continue is NOT an output filter
                // Chain members are included in output via link pointer traversal
                // Always print for debugging hit count discrepancy
                eprintln!("[DEBUG OUTPUT_FILTER] Total linked HSPs: {}", total_linked);
                eprintln!("[DEBUG OUTPUT_FILTER] Single HSPs (linked_set=false): {}", stats_single_hsps);
                eprintln!("[DEBUG OUTPUT_FILTER] Chain heads (linked_set=true, start_of_chain=true): {}", stats_chain_heads);
                eprintln!("[DEBUG OUTPUT_FILTER] Chain members (linked_set=true, start_of_chain=false): {} (included in output)", stats_chain_members);
                eprintln!("[DEBUG OUTPUT_FILTER] Filtered by E-value (threshold={}): {}", evalue_threshold, filtered_by_evalue);
                eprintln!("[DEBUG OUTPUT_FILTER] Expected after E-value filter: {} (all HSPs - E-value filtered)", total_linked - filtered_by_evalue);
                eprintln!("[DEBUG OUTPUT_FILTER] Final hits after filtering: {}", final_hits.len());

                if !final_hits.is_empty() {
                    if diag_enabled {
                        diagnostics
                            .base
                            .hsps_after_overlap_filter
                            .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                        diagnostics
                            .output_from_ungapped
                            .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                    }
                    st.tx.send(final_hits).unwrap();
                }
            }
            bar.inc(1);
            st.diag_offset = diag_offset;
        },
    );

    // Close channel so the writer can exit.
    drop(tx);

    bar.finish();
    writer.join().unwrap()?;
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }

    if timing_enabled {
        let t_search = t_search_start.elapsed();
        let scan_s = scan_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let scan_n = scan_calls.load(AtomicOrdering::Relaxed);
        let ungapped_s = ungapped_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let ungapped_n = ungapped_calls.load(AtomicOrdering::Relaxed);
        let reeval_s = reeval_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let reeval_n = reeval_calls.load(AtomicOrdering::Relaxed);
        let linking_s = linking_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let linking_n = linking_calls.load(AtomicOrdering::Relaxed);
        let identity_s = identity_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let identity_n = identity_calls.load(AtomicOrdering::Relaxed);

        eprintln!("[TIMING] read_queries: {:.3}s", t_read_queries.as_secs_f64());
        eprintln!("[TIMING] build_lookup: {:.3}s", t_build_lookup.as_secs_f64());
        eprintln!("[TIMING] read_subjects: {:.3}s", t_read_subjects.as_secs_f64());
        eprintln!("[TIMING] scan_subject: {:.3}s (calls={})", scan_s, scan_n);
        eprintln!("[TIMING] ungapped_extend: {:.3}s (calls={})", ungapped_s, ungapped_n);
        eprintln!("[TIMING] reevaluate: {:.3}s (calls={})", reeval_s, reeval_n);
        eprintln!("[TIMING] sum_stats_linking: {:.3}s (calls={})", linking_s, linking_n);
        eprintln!("[TIMING] identity_calc: {:.3}s (calls={})", identity_s, identity_n);
        eprintln!("[TIMING] search_total: {:.3}s", t_search.as_secs_f64());
        eprintln!("[TIMING] total: {:.3}s", t_total.elapsed().as_secs_f64());
    }
    Ok(())
}

/// Run TBLASTX with pre-computed neighbor map using subject-side indexing.
/// This approach:
/// 1. Index query k-mers: query_lookup[kmer] = [(q_idx, f_idx, pos), ...]
/// 2. Pre-compute neighbor relationships once
/// 3. For each subject k-mer, find all matching query positions via expanded lookup
/// 4. Apply NCBI-style sum_stats_even_gap_linking for HSP merging
///
/// Two-hit/diag logic is ported verbatim from NCBI aa_ungapped.c:s_BlastAaWordFinder_TwoHit
fn run_with_neighbor_map(args: TblastxArgs) -> Result<()> {
    use super::lookup::LOOKUP_ALPHABET_SIZE;

    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    // NCBI: window = diag->window
    let window: i32 = args.window_size as i32;
    // NCBI: wordsize = lookup->word_length (always 3 for protein/tblastx)
    let wordsize: i32 = 3;
    let threshold = args.threshold;

    // NCBI BLAST computes x_dropoff per-context using kbp[context]->Lambda:
    //   p->cutoffs[context].x_dropoff_init =
    //       (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // Reference: blast_parameters.c:219-221
    //
    // For tblastx, NCBI uses kbp_ideal->Lambda for all contexts (blast_stat.c:2796-2797).
    // x_dropoff_per_context is populated per-subject after counting total_contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    if args.verbose {
        eprintln!("Reading queries (neighbor map mode)...");
    }
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();

    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();

    // Apply SEG masking if enabled
    if args.seg {
        let seg_masker = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        const X_MASK_NCBISTDAA: u8 = 21;
        for frames in query_frames.iter_mut() {
            for frame in frames.iter_mut() {
                if frame.aa_seq.len() >= 3 {
                    for m in seg_masker.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1);
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
    }

    if args.verbose {
        eprintln!("Building neighbor lookup...");
    }
    let neighbor_lookup = NeighborLookup::build(&query_frames, threshold, &ungapped_params_for_xdrop);
    
    // Use compressed neighbor index: no expanded_lookup pre-computation
    // Instead, resolve neighbors on-the-fly during scan
    if args.verbose {
        eprintln!("Using compressed neighbor index (no expanded_lookup)...");
    }

    if args.verbose {
        eprintln!("Reading subjects...");
    }
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    // NCBI BLAST Karlin params for TBLASTX (ungapped-only algorithm):
    // 
    // TBLASTX is explicitly ungapped-only (blast_options.c line 869-873):
    //   "Gapped search is not allowed for tblastx"
    //
    // For bit score and E-value calculation, NCBI uses sbp->kbp (ungapped):
    //   blast_hits.c line 1833: kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    //   blast_hits.c line 1918: same pattern in Blast_HSPListGetBitScores
    //
    // For cutoff score search space calculation, NCBI uses gapped params:
    //   blast_parameters.c uses kbp_gap for eff_searchsp in cutoff calculation
    //
    // BLOSUM62 ungapped: lambda=0.3176, K=0.134
    // BLOSUM62 gapped:   lambda=0.267,  K=0.041
    let ungapped_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    // Gapped params for cutoff score search space calculation
    let gapped_params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    // Use UNGAPPED params for bit score and E-value (NCBI parity for tblastx)
    let params = ungapped_params.clone();
    
    // NCBI reference (blast_setup.c:768):
    //   kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
    // For tblastx, gapped_calculation = FALSE, so kbp_ptr = sbp->kbp (ungapped params)
    // Therefore, ALL calculations (eff_searchsp, cutoff, bit score, E-value) use UNGAPPED params
    //
    // BLOSUM62 ungapped: lambda=0.3176, K=0.134 (used for tblastx)
    // BLOSUM62 gapped:   lambda=0.267,  K=0.041 (NOT used for tblastx)
    // Note: gapped_params is kept for API compatibility but NOT used for tblastx
    
    // Optional timing breakdown (disabled by default to preserve output/parity logs)
    let timing_enabled = std::env::var_os("LOSAT_TIMING").is_some();
    let reeval_ns = AtomicU64::new(0);
    let reeval_calls = AtomicU64::new(0);
    let linking_ns = AtomicU64::new(0);
    let linking_calls = AtomicU64::new(0);
    let identity_ns = AtomicU64::new(0);
    let identity_calls = AtomicU64::new(0);

    // Compute NCBI-style average query length for linking cutoffs
    // Reference: blast_parameters.c:CalculateLinkHSPCutoffs (lines 1023-1026)
    let query_nucl_lengths: Vec<usize> = queries_raw.iter().map(|r| r.seq().len()).collect();
    let avg_query_length = compute_avg_query_length_ncbi(&query_nucl_lengths);

    if args.verbose {
        eprintln!(
            "Searching {} queries vs {} subjects (neighbor map mode, avg_query_length={})...",
            queries_raw.len(),
            subjects_raw.len(),
            avg_query_length
        );
    }

    let bar = if args.verbose {
        let bar = ProgressBar::new(subjects_raw.len() as u64 * 6);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
                .unwrap(),
        );
        bar
    } else {
        ProgressBar::hidden()
    };

    let evalue_threshold = args.evalue;
    let query_frames_ref = &query_frames;
    let query_ids_ref = &query_ids;
    let neighbor_map_ref = &neighbor_lookup.neighbor_map.map;
    let query_lookup_ref = &neighbor_lookup.query_lookup;
    let query_contexts_ref = &neighbor_lookup.contexts;
    let params_ref = &params;
    let ungapped_params_ref = &ungapped_params;
    let _gapped_params_ref = &gapped_params;  // Unused - tblastx uses ungapped params

    // Collect UngappedHit for sum_stats_linking
    // Key: (q_idx, s_idx)
    let all_ungapped: std::sync::Mutex<Vec<UngappedHit>> = std::sync::Mutex::new(Vec::new());

    // Process each subject
    for (s_idx, s_rec) in subjects_raw.iter().enumerate() {
        let s_id = s_rec
            .id()
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();
        let s_len = s_rec.seq().len();
        let s_frames = generate_frames(s_rec.seq(), &db_code);
        
        // Precompute per-context cutoff scores using NCBI BLAST algorithm.
        // Reference: blast_parameters.c BlastInitialWordParametersUpdate
        //
        // NCBI cutoff calculation for tblastx ungapped path:
        // 1. gap_trigger from ungapped params (kbp_std)
        // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
        // 3. Per-subject update: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
        // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
        //
        // For -subject mode, subject_len_nucl is used (not per-frame AA length).
        let subject_len_nucl = s_len as i64;
        let total_contexts: usize = query_frames.iter().map(|f| f.len()).sum();
        let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; total_contexts];
        
        // =======================================================================
        // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
        // =======================================================================
        // NCBI stores these in query_info->contexts[ctx].length_adjustment and
        // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths (blast_setup.c:700-850).
        // We precompute them here and pass to sum_stats_linking for NCBI parity.
        // Reference: blast_setup.c:846-847
        //   query_info->contexts[index].eff_searchsp = effective_search_space;
        //   query_info->contexts[index].length_adjustment = length_adjustment;
        let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
        let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
        
        // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
        // Compute per-context x_dropoff using kbp[context]->Lambda.
        // Reference: blast_parameters.c:219-221
        //
        // For tblastx, all contexts use kbp_ideal (BLOSUM62 ungapped Lambda = 0.3176),
        // so x_dropoff = ceil(7.0 * LN2 / 0.3176) = 16 for all contexts.
        // We maintain per-context structure for NCBI parity.
        let x_dropoff_per_context: Vec<i32> = (0..total_contexts)
            .map(|_| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ungapped_params_for_xdrop, 1.0))
            .collect();
        
        // Precompute gap_trigger once (same for all contexts with same params)
        // Reference: blast_parameters.c:343-344
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ungapped_params_ref);
        
        let mut ctx_idx = 0;
        for q_frames in query_frames.iter() {
            for q_frame in q_frames.iter() {
                // NCBI: query_length = query_info->contexts[context].query_length
                let query_len_aa = q_frame.aa_len as i64;
                
                // =======================================================================
                // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                // length_adjustment and eff_searchsp from a single source of truth.
                // This mirrors BLAST_CalcEffLengths which computes and stores both values.
                // Reference: blast_setup.c:821-847
                // =======================================================================
                let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    ungapped_params_ref,  // tblastx uses ungapped params (kbp_gap is NULL)
                );
                let eff_searchsp = eff_lengths.eff_searchsp;
                length_adj_per_context.push(eff_lengths.length_adjustment);
                eff_searchsp_per_context.push(eff_searchsp);
                
                // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
                // This uses the effective search space WITH length adjustment
                // Reference: blast_parameters.c:942-946
                let cutoff_score_max = cutoff_score_max_for_tblastx(
                    eff_searchsp,
                    evalue_threshold,  // User's E-value (typically 10.0)
                    ungapped_params_ref,
                );
                
                // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
                // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
                // Reference: blast_parameters.c:348-374
                let cutoff = cutoff_score_for_update_tblastx(
                    query_len_aa,
                    subject_len_nucl,  // NUCLEOTIDE length, NOT divided by 3!
                    gap_trigger,
                    cutoff_score_max,
                    BLAST_GAP_DECAY_RATE,  // 0.5
                    ungapped_params_ref,
                    1.0,  // scale_factor (standard BLOSUM62)
                );
                
                // All subject frames use the same cutoff (NCBI behavior)
                for sf_idx in 0..s_frames.len() {
                    cutoff_scores[ctx_idx][sf_idx] = cutoff;
                }
                ctx_idx += 1;
            }
        }
        let cutoff_scores_ref = &cutoff_scores;
        let x_dropoff_per_context_ref = &x_dropoff_per_context;
        let length_adj_per_context_ref = &length_adj_per_context;
        let eff_searchsp_per_context_ref = &eff_searchsp_per_context;

        // Process each subject frame in parallel
        let frame_hits: Vec<UngappedHit> = s_frames
            .par_iter()
            .enumerate()
            .flat_map(|(s_f_idx, s_frame)| {
                let s_aa = &s_frame.aa_seq;
                if s_aa.len() < 5 {
                    bar.inc(1);
                    return Vec::new();
                }

                // NCBI: BlastSaveInitHsp equivalent - store initial HSPs with absolute coordinates
                // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let s_aa_len = s_aa.len();
                
                // Count total query contexts for Vec-based diagonal tracking
                let total_q_contexts: usize = query_frames_ref.iter().map(|f| f.len()).sum();
                let max_q_aa_len = query_frames_ref.iter()
                    .flat_map(|frames| frames.iter().map(|f| f.aa_seq.len()))
                    .max()
                    .unwrap_or(1);
                
                // Compute diag_array dimensions
                let diag_array_len = s_aa_len + max_q_aa_len + 1;
                let diag_offset = max_q_aa_len as i32;
                
                // Build context index mapping: (q_idx, q_f_idx) -> flat index
                let mut ctx_base: Vec<usize> = Vec::with_capacity(query_frames_ref.len() + 1);
                ctx_base.push(0);
                for frames in query_frames_ref.iter() {
                    ctx_base.push(ctx_base.last().unwrap() + frames.len());
                }
                
                // NCBI-style DiagStruct arrays: [ctx_flat_idx][diag] = DiagStruct { last_hit, flag }
                // NCBI: last_hit stores (subject_offset + diag_offset) to avoid subtraction in hot path
                // NCBI: flag=1 means we just extended, flag=0 means ready for two-hit check
                let mut diag_array: Vec<Vec<DiagStruct>> = vec![vec![DiagStruct::default(); diag_array_len]; total_q_contexts];
                
                // NCBI aa_ungapped.c:502-610 - main scan loop
                for s_pos in 1..=(s_aa_len.saturating_sub(4)) {
                    let c0 = unsafe { *s_aa.get_unchecked(s_pos) } as usize;
                    let c1 = unsafe { *s_aa.get_unchecked(s_pos + 1) } as usize;
                    let c2 = unsafe { *s_aa.get_unchecked(s_pos + 2) } as usize;

                    if c0 >= LOOKUP_ALPHABET_SIZE || c1 >= LOOKUP_ALPHABET_SIZE || c2 >= LOOKUP_ALPHABET_SIZE {
                        continue;
                    }

                    let s_kmer = encode_kmer(c0, c1, c2);
                    
                    // Compressed neighbor index: resolve neighbors on-the-fly
                    let neighbor_kmers = &neighbor_map_ref[s_kmer];
                    if neighbor_kmers.is_empty() {
                        continue;
                    }
                    
                    // subject_offset in NCBI terms (0-based position in subject AA sequence)
                    let subject_offset: i32 = (s_pos - 1) as i32;
                    
                    // For each neighbor query k-mer, get query positions
                    for &neighbor_kmer in neighbor_kmers {
                        let query_hits = &query_lookup_ref[neighbor_kmer as usize];
                        for &(q_idx, q_f_idx, q_pos) in query_hits {
                            let ctx_flat = ctx_base[q_idx as usize] + q_f_idx as usize;
                            
                            // query_offset in NCBI terms (0-based position in query AA sequence)
                            let query_offset: i32 = q_pos as i32;
                            
                            // NCBI: diag_coord = (query_offset - subject_offset) & diag_mask
                            // We use a Vec so we compute a direct index instead of masking
                            let diag_coord = (query_offset - subject_offset + diag_offset) as usize;
                            
                            if diag_coord >= diag_array_len {
                                continue;
                            }
                            
                            let diag_entry = &mut diag_array[ctx_flat][diag_coord];
                            
                            // NCBI aa_ungapped.c:519-531 - flag==1 block
                            // "If the reset bit is set, an extension just happened."
                            if diag_entry.flag != 0 {
                                // "If we've already extended past this hit, skip it."
                                if subject_offset + diag_offset < diag_entry.last_hit {
                                    continue;
                                }
                                // "Otherwise, start a new hit." - reset flag and last_hit, then FALL THROUGH
                                // NCBI: After reset, control falls through to the else block (flag==0) below
                                // to continue with extension logic. DO NOT continue here!
                                diag_entry.last_hit = subject_offset + diag_offset;
                                diag_entry.flag = 0;
                                // NO continue here - fall through to extension logic below
                            }
                            
                            // NCBI aa_ungapped.c:533-606 - flag==0 block
                            // "If the reset bit is cleared, try to start an extension."
                            
                            // NCBI: last_hit = diag_array[diag_coord].last_hit - diag_offset
                            let last_hit = diag_entry.last_hit - diag_offset;
                            // NCBI: diff = subject_offset - last_hit
                            let diff = subject_offset - last_hit;
                            
                            // NCBI aa_ungapped.c:538-544: "if (diff >= window)"
                            if diff >= window {
                                // "We are beyond the window for this diagonal; start a new hit"
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:549-551: "if (diff < wordsize)"
                            if diff < wordsize {
                                // "If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up"
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:560-573: Check if last hit is in current query context
                            // "if (query_offset - diff < query_info->contexts[curr_context].query_offset)"
                            // In neighbor-map mode, each (q_idx, q_f_idx) is its own context with frame_base=0
                            // So we check: query_offset - diff < 0
                            if query_offset - diff < 0 {
                                // "there was no last hit for this diagonal; start a new hit"
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // PASSED all two-hit checks - now extend
                            let q_frame = &query_frames_ref[q_idx as usize][q_f_idx as usize];
                            let q_aa = &q_frame.aa_seq;
                            
                            // Convert to raw positions (+1 for sentinel)
                            let q_raw = (query_offset + 1) as usize;
                            let s_raw = (subject_offset + 1) as usize;
                            
                            if q_raw + 2 >= q_aa.len() || s_raw + 2 >= s_aa.len() {
                                continue;
                            }
                            
                            // [C] cutoffs->x_dropoff (per-context x_dropoff)
                            // Reference: aa_ungapped.c:579
                            let x_dropoff = x_dropoff_per_context_ref[ctx_flat];
                            
                            // NCBI aa_ungapped.c:576-583: s_BlastAaExtendTwoHit
                            // s_left_off = last_hit + wordsize (end of first hit word)
                            // s_right_off = subject_offset (second hit start)
                            // q_right_off = query_offset (second hit start)
                            // last_hit = subject_offset which is 0-indexed, so add +1 for array index
                            // (consistent with s_raw and q_raw which also have +1)
                            let s_left_off = (last_hit + wordsize + 1) as usize;
                            let (hsp_q, hsp_qe, hsp_s, _hsp_se, score, right_extend, s_last_off) =
                                extend_hit_two_hit(
                                    q_aa,
                                    s_aa,
                                    s_left_off,      // s_left_off (end of first hit word, raw coords)
                                    s_raw,           // s_right_off (second hit start, raw coords)
                                    q_raw,           // q_right_off (second hit start, raw coords)
                                    x_dropoff,       // x_dropoff (per-context, NCBI parity)
                                );
                            
                            let hsp_len = hsp_qe.saturating_sub(hsp_q);
                            
                            // NCBI aa_ungapped.c:596-606: Update diag_array after extension
                            if right_extend {
                                // "If an extension to the right happened, reset the last hit so that
                                //  future hits to this diagonal must start over."
                                diag_entry.flag = 1;
                                diag_entry.last_hit = (s_last_off as i32) - (wordsize - 1) + diag_offset;
                            } else {
                                // "Otherwise, make the present hit into the previous hit for this diagonal"
                                diag_entry.last_hit = subject_offset + diag_offset;
                            }
                            
                            // NCBI aa_ungapped.c:588-591: Check score threshold
                            // cutoffs = word_params->cutoffs + curr_context
                            // if (score >= cutoffs->cutoff_score) ...
                            let cutoff = cutoff_scores_ref[ctx_flat][s_f_idx];
                            if score < cutoff {
                                continue;
                            }
                            
                            // NCBI: BlastSaveInitHsp equivalent
                            // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                            // Store HSP with absolute coordinates (before coordinate conversion)
                            //
                            // In neighbor-map mode, each query frame is independent with frame_base=0
                            // hsp_q is frame-relative coordinate (with sentinel)
                            // Calculate absolute coordinate in concatenated buffer
                            let frame_base = 0i32;  // Each query frame is independent in neighbor-map mode
                            let hsp_q_absolute = frame_base + (hsp_q as i32) - 1;  // -1 for sentinel
                            let hsp_qe_absolute = frame_base + ((hsp_q + hsp_len) as i32) - 1;  // -1 for sentinel
                            
                            init_hsps.push(InitHSP {
                                q_start_absolute: hsp_q_absolute,
                                q_end_absolute: hsp_qe_absolute,
                                s_start: hsp_s as i32,
                                s_end: (hsp_s + hsp_len) as i32,
                                score,
                                ctx_idx: ctx_flat,
                                s_f_idx,
                                q_idx,
                                s_idx: s_idx as u32,
                                q_frame: q_frame.frame,
                                s_frame: s_frame.frame,
                                q_orig_len: q_frame.orig_len,
                                s_orig_len: s_len,
                            });
                        } // end for query_hits
                    } // end for neighbor_kmers
                } // end for s_pos

                // NCBI: BLAST_GetUngappedHSPList equivalent
                // Reference: blast_gapalign.c:4719-4775
                // Convert InitHSPs (with absolute coordinates) to UngappedHits (with context-relative coordinates)
                // and perform batch reevaluation
                //
                // In neighbor-map mode, each query frame is independent with frame_base=0
                // Create temporary contexts for get_ungapped_hsp_list
                let mut temp_contexts = Vec::new();
                for (q_idx, frames) in query_frames_ref.iter().enumerate() {
                    for (f_idx, frame) in frames.iter().enumerate() {
                        temp_contexts.push(super::lookup::QueryContext {
                            q_idx: q_idx as u32,
                            f_idx: f_idx as u8,
                            frame: frame.frame,
                            aa_seq: frame.aa_seq.clone(),
                            aa_seq_nomask: frame.aa_seq_nomask.clone(),
                            aa_len: frame.aa_len,
                            orig_len: frame.orig_len,
                            frame_base: 0,  // Each query frame is independent
                            karlin_params: ungapped_params_ref.clone(),
                        });
                    }
                }
                // NCBI: BLAST_GetUngappedHSPList equivalent
                // Convert InitHSPs (with absolute coordinates) to UngappedHits (with context-relative coordinates)
                let local_hits = get_ungapped_hsp_list(
                    init_hsps,
                    &temp_contexts,
                    &s_frames,
                );
                
                // NCBI: Blast_HSPListReevaluateUngapped equivalent
                // Perform batch reevaluation on all HSPs
                // Note: Statistics collection is disabled in neighbor-map mode (parallel processing)
                let mut dummy_stats = 0usize;
                let local_hits = reevaluate_ungapped_hsp_list(
                    local_hits,
                    &temp_contexts,
                    &s_frames,
                    cutoff_scores_ref,
                    &args,
                    timing_enabled,
                    &reeval_ns,
                    &reeval_calls,
                    false,  // is_long_sequence: statistics not collected in neighbor-map mode
                    &mut dummy_stats,
                );

                bar.inc(1);
                local_hits
            })
            .collect();

        // Add to global collection
        all_ungapped.lock().unwrap().extend(frame_hits);
    }

    bar.finish();

    let all_ungapped = all_ungapped.into_inner().unwrap();
    if args.verbose {
        eprintln!("=== Stage Counters ===");
        eprintln!("[1] Raw ungapped hits (after extension): {}", all_ungapped.len());
    }

    // NOTE: NCBI does NOT call Blast_HSPListPurgeHSPsWithCommonEndpoints in the
    // ungapped tblastx path. The purge is only in the gapped path:
    //   blast_engine.c line 545: inside `if (score_options->gapped_calculation)`
    // Since tblastx is ungapped, we skip purge for NCBI parity.
    //
    // Additionally, if purge were needed, it must be per-subject (BlastHSPList),
    // not on the entire all_ungapped collection which mixes multiple subjects.

    // Apply sum-stats linking per subject (NCBI CalculateLinkHSPCutoffs is subject-specific)
    // Reference: blast_parameters.c:CalculateLinkHSPCutoffs
    // Step 1: Group hits by s_idx
    let mut hits_by_subject: rustc_hash::FxHashMap<u32, Vec<UngappedHit>> = rustc_hash::FxHashMap::default();
    for hit in all_ungapped {
        hits_by_subject.entry(hit.s_idx).or_default().push(hit);
    }
    
    // Step 2: Apply linking per subject with proper cutoffs
    let linked_hits: Vec<UngappedHit> = hits_by_subject
        .into_par_iter()
        .flat_map(|(s_idx, subject_hits)| {
            // Get subject length for this subject
            let subject_len_nucl = subjects_raw[s_idx as usize].seq().len() as i64;
            
            // =======================================================================
            // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
            // =======================================================================
            // NCBI stores these in query_info->contexts[ctx].length_adjustment and
            // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths (blast_setup.c:700-850).
            // We precompute them here for each subject and pass to sum_stats_linking.
            let total_contexts: usize = query_frames_ref.iter().map(|f| f.len()).sum();
            let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
            let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
            
            // Compute cutoff_score_min as minimum across all query contexts
            // using the NCBI BlastInitialWordParametersUpdate logic
            // Reference: blast_parameters.c:348-374, 401-403
            let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ungapped_params_ref);
            let mut cutoff_score_min = i32::MAX;
            for frames in query_frames_ref.iter() {
                for frame in frames.iter() {
                    let query_len_aa = frame.aa_len as i64;
                    
                    // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                    // length_adjustment and eff_searchsp from a single source of truth.
                    // Reference: blast_setup.c:821-847
                    let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                        query_len_aa,
                        subject_len_nucl,
                        ungapped_params_ref,  // tblastx uses ungapped params (kbp_gap is NULL)
                    );
                    let eff_searchsp = eff_lengths.eff_searchsp;
                    length_adj_per_context.push(eff_lengths.length_adjustment);
                    eff_searchsp_per_context.push(eff_searchsp);
                    
                    // Step 1: cutoff_score_max from BlastHitSavingParametersNew
                    let cutoff_score_max = cutoff_score_max_for_tblastx(
                        eff_searchsp,
                        evalue_threshold,
                        ungapped_params_ref,
                    );
                    
                    // Step 2: Per-subject cutoff from BlastInitialWordParametersUpdate
                    let cutoff = cutoff_score_for_update_tblastx(
                        query_len_aa,
                        subject_len_nucl,
                        gap_trigger,
                        cutoff_score_max,
                        BLAST_GAP_DECAY_RATE,
                        ungapped_params_ref,
                        1.0,
                    );
                    cutoff_score_min = cutoff_score_min.min(cutoff);
                }
            }
            if cutoff_score_min == i32::MAX {
                cutoff_score_min = 0;
            }
            
            // Create NCBI-style linking params for this subject
            let linking_params = LinkingParams {
                avg_query_length,
                subject_len_nucl,
                cutoff_score_min,
                scale_factor: 1.0,
                gap_decay_rate: 0.5,
            };

            // Subject frame bases in the concatenated translation buffer order.
            // We match `generate_frames()` ordering (forward 1,2,3 then reverse -1,-2,-3)
            // and each frame has its own leading+trailing sentinel in LOSAT's current model.
            let subject_len = subject_len_nucl as usize;
            let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(6);
            let mut base: i32 = 0;
            for shift in 0..3usize {
                if shift + 3 <= subject_len {
                    let aa_len = ((subject_len - shift) / 3) as i32;
                    subject_frame_bases.push(base);
                    // NCBI concatenation shares the boundary NULLB between frames:
                    // increment is (aa_len + 1), not (aa_len + 2).
                    base += aa_len + 1;
                }
            }
            for shift in 0..3usize {
                if shift + 3 <= subject_len {
                    let aa_len = ((subject_len - shift) / 3) as i32;
                    subject_frame_bases.push(base);
                    base += aa_len + 1;
                }
            }

            // NCBI parity: s_BlastFindSmallestLambda selects Karlin params with smallest lambda
            // Reference: blast_parameters.c:92-112, 1012-1013
            // For tblastx, all contexts use kbp_ideal (same lambda), so this is equivalent
            // to using params_ref directly. We maintain this structure for parity.
            let context_params: Vec<KarlinParams> = query_contexts_ref.iter().map(|_| params_ref.clone()).collect();
            let linking_params_for_cutoff = find_smallest_lambda_params(&context_params)
                .unwrap_or_else(|| params_ref.clone());

            // NCBI Parity: Pass pre-computed length_adjustment and eff_searchsp per context
            // These values are stored in query_info->contexts[ctx] in NCBI and referenced
            // by link_hsps.c for BLAST_SmallGapSumE/BLAST_LargeGapSumE calculations.
            let t_linking = if timing_enabled { Some(Instant::now()) } else { None };
            let linked = apply_sum_stats_even_gap_linking(
                subject_hits,
                &linking_params_for_cutoff,
                &linking_params,
                query_contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
            );
            if let Some(t) = t_linking {
                let elapsed = t.elapsed();
                linking_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                linking_calls.fetch_add(1, AtomicOrdering::Relaxed);
            }
            linked
        })
        .collect();

    if trace_hsp_target().is_some() {
        for h in &linked_hits {
            trace_ungapped_hit_if_match("after_linking", h);
        }
    }
    if args.verbose {
        eprintln!("[2] After sum_stats_linking: {} hits", linked_hits.len());
    }
    
    if args.verbose {
        // Count E-value distribution before final filter
        let mut e_0 = 0usize;      // E <= 0 (0.0)
        let mut e_tiny = 0usize;   // 0 < E <= 1e-50
        let mut e_small = 0usize;  // 1e-50 < E <= 1e-10
        let mut e_med = 0usize;    // 1e-10 < E <= 10
        let mut e_large = 0usize;  // E > 10
        let mut e_inf = 0usize;    // E == INF
        for h in &linked_hits {
            if h.e_value == f64::INFINITY {
                e_inf += 1;
            } else if h.e_value <= 0.0 {
                e_0 += 1;
            } else if h.e_value <= 1e-50 {
                e_tiny += 1;
            } else if h.e_value <= 1e-10 {
                e_small += 1;
            } else if h.e_value <= 10.0 {
                e_med += 1;
            } else {
                e_large += 1;
            }
        }
        eprintln!("[3] E-value distribution (before final filter):");
        eprintln!(
            "    E=0: {}, E<=1e-50: {}, E<=1e-10: {}, E<=10: {}, E>10: {}, E=INF: {}",
            e_0, e_tiny, e_small, e_med, e_large, e_inf
        );
    }

    // Pre-generate subject frames to avoid redundant translation
    if args.verbose {
        eprintln!("Pre-generating subject frames for identity calculation...");
    }
    let subject_frames_cache: Vec<Vec<QueryFrame>> = subjects_raw
        .iter()
        .map(|s| generate_frames(s.seq(), &db_code))
        .collect();
    
    // Convert UngappedHit to Hit
    let mut final_hits: Vec<Hit> = Vec::with_capacity(linked_hits.len());
    
    for h in linked_hits {
        if h.e_value > evalue_threshold {
            continue;
        }
        
        let q_id = &query_ids_ref[h.q_idx as usize];
        let s_id = subjects_raw[h.s_idx as usize]
            .id()
            .split_whitespace()
            .next()
            .unwrap_or("unknown");
        
        // `ctx_idx` is the global query context index (NCBI `hsp->context` equivalent).
        // Use the pre-built `QueryContext` buffer to access the (possibly unmasked) query AA sequence.
        let q_ctx = &query_contexts_ref[h.ctx_idx];
        let q_aa = q_ctx.aa_seq_nomask.as_deref().unwrap_or(&q_ctx.aa_seq);
        
        // Get subject frame from cache (no redundant translation)
        let s_frame_obj = &subject_frames_cache[h.s_idx as usize][h.s_f_idx];
        let s_aa = &s_frame_obj.aa_seq;
        
        let len = h.q_aa_end.saturating_sub(h.q_aa_start);
        
        // Calculate identity using raw positions (+1 for sentinel)
        let q0 = h.q_aa_start + 1;
        let s0 = h.s_aa_start + 1;
        // Use SIMD-optimized identity calculation (already implemented in reevaluate.rs)
        let t_identity = if timing_enabled { Some(Instant::now()) } else { None };
        let matches = get_num_identities_and_positives_ungapped(
            q_aa,
            s_aa,
            q0,
            s0,
            len,
        ).map(|(num_ident, _)| num_ident).unwrap_or(0);
        if let Some(t) = t_identity {
            let elapsed = t.elapsed();
            identity_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
            identity_calls.fetch_add(1, AtomicOrdering::Relaxed);
        }
        let identity = if len > 0 {
            (matches as f64 / len as f64) * 100.0
        } else {
            0.0
        };

        let bit_score = calc_bit_score(h.raw_score, &params);
        let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, h.q_frame, h.q_orig_len);
        let (s_start, s_end) = convert_coords(h.s_aa_start, h.s_aa_end, h.s_frame, h.s_orig_len);

        let out_hit = Hit {
            query_id: q_id.clone(),
            subject_id: s_id.to_string(),
            identity,
            length: len,
            mismatch: len.saturating_sub(matches),
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: h.e_value,
            bit_score,
            q_idx: h.q_idx,
            s_idx: h.s_idx,
            raw_score: h.raw_score,
        };
        trace_final_hit_if_match("output_hit", &out_hit);
        final_hits.push(out_hit);
    }

    if args.verbose {
        eprintln!(
            "[4] Final hits after E-value filter (threshold={}): {}",
            evalue_threshold,
            final_hits.len()
        );
    }
    
    // Report top hit alignment length
    if !final_hits.is_empty() {
        let max_len = final_hits.iter().map(|h| h.length).max().unwrap_or(0);
        if args.verbose {
            eprintln!("[5] Top alignment length: {} AA", max_len);
        }
    }
    if args.verbose {
        eprintln!("=== End Stage Counters ===");
    }
    
    if timing_enabled {
        let linking_s = linking_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let linking_n = linking_calls.load(AtomicOrdering::Relaxed);
        let identity_s = identity_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let identity_n = identity_calls.load(AtomicOrdering::Relaxed);
        eprintln!("[TIMING] sum_stats_linking: {:.3}s (calls={})", linking_s, linking_n);
        eprintln!("[TIMING] identity_calc: {:.3}s (calls={})", identity_s, identity_n);
    }
    
    // NCBI-style output ordering: query (input order) → subject (best_evalue/score/oid) → HSP (score/coords)
    // Reference: BLAST_LinkHsps() + s_EvalueCompareHSPLists() + ScoreCompareHSPs()
    write_output_ncbi_order(final_hits, args.out.as_ref())?;
    Ok(())
}
