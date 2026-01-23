//! Neighbor map mode implementation
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use super::*;
use std::sync::Arc;

/// Run TBLASTX with pre-computed neighbor map using subject-side indexing.
/// This approach:
/// 1. Index query k-mers: query_lookup[kmer] = [(q_idx, f_idx, pos), ...]
/// 2. Pre-compute neighbor relationships once
/// 3. For each subject k-mer, find all matching query positions via expanded lookup
/// 4. Apply NCBI-style sum_stats_even_gap_linking for HSP merging
///
/// Two-hit/diag logic is ported verbatim from NCBI aa_ungapped.c:s_BlastAaWordFinder_TwoHit
pub(crate) fn run_with_neighbor_map(args: TblastxArgs) -> Result<()> {
    use crate::algorithm::tblastx::lookup::LOOKUP_ALPHABET_SIZE;

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
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes kbp_std per context and applies check_ideal:
    //   if (check_ideal && kbp->Lambda >= kbp_ideal->Lambda) Blast_KarlinBlkCopy(kbp, kbp_ideal);
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
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
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    let query_ids: Vec<Arc<str>> = queries_raw
        .iter()
        .map(|r| Arc::<str>::from(r.id().split_whitespace().next().unwrap_or("unknown")))
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

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    let subject_ids: Vec<Arc<str>> = subjects_raw
        .iter()
        .map(|r| Arc::<str>::from(r.id().split_whitespace().next().unwrap_or("unknown")))
        .collect();

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
    //   blast_parameters.c selects kbp_gap if available, else kbp for tblastx.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:860-865
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
    
    // NCBI reference (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:768):
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
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
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
    let neighbor_map_ref = &neighbor_lookup.neighbor_map.map;
    let query_lookup_ref = &neighbor_lookup.query_lookup;
    let query_contexts_ref = &neighbor_lookup.contexts;
    let params_ref = &params;
    let _gapped_params_ref = &gapped_params;  // Unused - tblastx uses ungapped params

    // Collect UngappedHit for sum_stats_linking
    // Key: (q_idx, s_idx)
    let all_ungapped: std::sync::Mutex<Vec<UngappedHit>> = std::sync::Mutex::new(Vec::new());

    // Process each subject
    for (s_idx, s_rec) in subjects_raw.iter().enumerate() {
        let s_len = s_rec.seq().len();
        let s_frames = generate_frames(s_rec.seq(), &db_code);
        
        // Precompute per-context cutoff scores using NCBI BLAST algorithm.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:280-419
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
        // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths
        // (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:700-850).
        // We precompute them here and pass to sum_stats_linking for NCBI parity.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846-847
        //   query_info->contexts[index].eff_searchsp = effective_search_space;
        //   query_info->contexts[index].length_adjustment = length_adjustment;
        let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
        let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
        
        // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
        // Compute per-context x_dropoff using kbp[context]->Lambda.
        // NCBI uses kbp_std[context]->Lambda for x_dropoff_init per context.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:219-221
        let x_dropoff_per_context: Vec<i32> = neighbor_lookup
            .contexts
            .iter()
            .map(|ctx| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ctx.karlin_params, 1.0))
            .collect();
        
        let mut ctx_idx = 0;
        for q_frames in query_frames.iter() {
            for q_frame in q_frames.iter() {
                // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
                let ctx_params = &neighbor_lookup.contexts[ctx_idx].karlin_params;

                // NCBI: query_length = query_info->contexts[context].query_length
                let query_len_aa = q_frame.aa_len as i64;

                // NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-345
                let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params);
                
                // =======================================================================
                // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                // length_adjustment and eff_searchsp from a single source of truth.
                // This mirrors BLAST_CalcEffLengths which computes and stores both values.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847
                // =======================================================================
                let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    ctx_params,  // tblastx uses per-context ungapped params (kbp_gap is NULL)
                );
                let eff_searchsp = eff_lengths.eff_searchsp;
                length_adj_per_context.push(eff_lengths.length_adjustment);
                eff_searchsp_per_context.push(eff_searchsp);
                
                // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
                // This uses the effective search space WITH length adjustment
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:942-946
                let cutoff_score_max = cutoff_score_max_for_tblastx(
                    eff_searchsp,
                    evalue_threshold,  // User's E-value (typically 10.0)
                    ctx_params,
                );
                
                // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
                // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374
                let cutoff = cutoff_score_for_update_tblastx(
                    query_len_aa,
                    subject_len_nucl,  // NUCLEOTIDE length, NOT divided by 3!
                    gap_trigger,
                    cutoff_score_max,
                    BLAST_GAP_DECAY_RATE,  // 0.5
                    ctx_params,
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
                let s_aa_full = &s_frame.aa_seq;
                // NCBI subject->sequence points past the leading NULLB sentinel, so offsets
                // are 0-based from the first residue; see blast_engine.c:811-812 and
                // blast_util.c:112-116.
                let s_aa = &s_aa_full[1..s_aa_full.len() - 1];
                if s_aa.len() < wordsize as usize {
                    bar.inc(1);
                    return Vec::new();
                }

                // NCBI: BlastSaveInitHsp equivalent - store initial HSPs with absolute coordinates
                // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let s_aa_len = s_aa.len();
                
                // Count total query contexts for Vec-based diagonal tracking
                let total_q_contexts: usize = query_frames_ref.iter().map(|f| f.len()).sum();
                // NCBI query/subject lengths exclude sentinels.
                // Reference: blast_query_info.c:311-315, blast_engine.c:811-812.
                let max_q_aa_len = query_frames_ref.iter()
                    .flat_map(|frames| frames.iter().map(|f| f.aa_len))
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
                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:509-519.
                for s_pos in 0..=s_aa_len.saturating_sub(wordsize as usize) {
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
                    let subject_offset: i32 = s_pos as i32;
                    
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
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:536-553
                            // "If the reset bit is set, an extension just happened."
                            if diag_entry.flag() != 0 {
                                // "If we've already extended past this hit, skip it."
                                if subject_offset + diag_offset < diag_entry.last_hit() {
                                    continue;
                                }
                                // "Otherwise, start a new hit." - reset flag and last_hit, then FALL THROUGH
                                // NCBI: After reset, control falls through to the else block (flag==0) below
                                // to continue with extension logic. DO NOT continue here!
                                diag_entry.set_last_hit(subject_offset + diag_offset);
                                diag_entry.set_flag(0);
                                // NO continue here - fall through to extension logic below
                            }
                            
                            // NCBI aa_ungapped.c:533-606 - flag==0 block
                            // "If the reset bit is cleared, try to start an extension."
                            
                            // NCBI: last_hit = diag_array[diag_coord].last_hit - diag_offset
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:558-560
                            let last_hit = diag_entry.last_hit() - diag_offset;
                            // NCBI: diff = subject_offset - last_hit
                            let diff = subject_offset - last_hit;
                            
                            // NCBI aa_ungapped.c:538-544: "if (diff >= window)"
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:562-569
                            if diff >= window {
                                // "We are beyond the window for this diagonal; start a new hit"
                                diag_entry.set_last_hit(subject_offset + diag_offset);
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:549-551: "if (diff < wordsize)"
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:573-580
                            if diff < wordsize {
                                // "If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up"
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:560-573: Check if last hit is in current query context
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:592-606
                            // "if (query_offset - diff < query_info->contexts[curr_context].query_offset)"
                            // In neighbor-map mode, each (q_idx, q_f_idx) is its own context with frame_base=0
                            // So we check: query_offset - diff < 0
                            if query_offset - diff < 0 {
                                // "there was no last hit for this diagonal; start a new hit"
                                diag_entry.set_last_hit(subject_offset + diag_offset);
                                continue;
                            }
                            
                            // PASSED all two-hit checks - now extend
                            let q_frame = &query_frames_ref[q_idx as usize][q_f_idx as usize];
                            // NCBI query->sequence points past the leading NULLB sentinel, so
                            // query offsets are 0-based into that buffer.
                            // Reference: blast_query_info.c:311-315, blast_util.c:112-116.
                            let q_aa_full = &q_frame.aa_seq;
                            let q_aa = &q_aa_full[1..q_aa_full.len() - 1];
                            
                            let q_raw = query_offset as usize;
                            let s_raw = subject_offset as usize;
                            
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
                            let s_left_off = (last_hit + wordsize) as usize;
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
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:636-648
                            if right_extend {
                                // "If an extension to the right happened, reset the last hit so that
                                //  future hits to this diagonal must start over."
                                diag_entry.set_flag(1);
                                diag_entry.set_last_hit(
                                    (s_last_off as i32) - (wordsize - 1) + diag_offset,
                                );
                            } else {
                                // "Otherwise, make the present hit into the previous hit for this diagonal"
                                diag_entry.set_last_hit(subject_offset + diag_offset);
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
                            // In neighbor-map mode, each query frame is independent with frame_base=0.
                            // hsp_q is frame-relative coordinate in query->sequence (0-based).
                            // Reference: blast_gapalign.c:4756-4768, blast_query_info.c:311-315.
                            let frame_base = 0i32;  // Each query frame is independent in neighbor-map mode
                            let hsp_q_absolute = frame_base + (hsp_q as i32);
                            let hsp_qe_absolute = frame_base + ((hsp_q + hsp_len) as i32);
                            
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
                // NCBI uses per-context kbp_std computed from query composition.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
                let mut temp_contexts: Vec<QueryContext> = Vec::with_capacity(query_contexts_ref.len());
                for ctx in query_contexts_ref.iter() {
                    let mut ctx_copy = ctx.clone();
                    // Neighbor-map mode treats each query frame independently.
                    ctx_copy.frame_base = 0;
                    temp_contexts.push(ctx_copy);
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
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
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
            // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths
            // (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:700-850).
            // We precompute them here for each subject and pass to sum_stats_linking.
            let total_contexts: usize = query_contexts_ref.len();
            let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
            let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
            
            // Compute cutoff_score_min as minimum across all query contexts
            // using the NCBI BlastInitialWordParametersUpdate logic
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374, 401-403
            let mut cutoff_score_min = i32::MAX;
            for ctx in query_contexts_ref.iter() {
                // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
                let ctx_params = &ctx.karlin_params;
                let query_len_aa = ctx.aa_len as i64;

                // NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-345
                let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params);

                // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                // length_adjustment and eff_searchsp from a single source of truth.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847
                let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    ctx_params,  // tblastx uses per-context ungapped params (kbp_gap is NULL)
                );
                let eff_searchsp = eff_lengths.eff_searchsp;
                length_adj_per_context.push(eff_lengths.length_adjustment);
                eff_searchsp_per_context.push(eff_searchsp);

                // Step 1: cutoff_score_max from BlastHitSavingParametersNew
                let cutoff_score_max = cutoff_score_max_for_tblastx(
                    eff_searchsp,
                    evalue_threshold,
                    ctx_params,
                );

                // Step 2: Per-subject cutoff from BlastInitialWordParametersUpdate
                let cutoff = cutoff_score_for_update_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    gap_trigger,
                    cutoff_score_max,
                    BLAST_GAP_DECAY_RATE,
                    ctx_params,
                    1.0,
                );
                cutoff_score_min = cutoff_score_min.min(cutoff);
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

            // NCBI parity: s_BlastFindSmallestLambda selects smallest lambda across contexts.
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:92-112
            // Per-context kbp_std is computed from query composition (with check_ideal).
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
            let context_params: Vec<KarlinParams> = query_contexts_ref
                .iter()
                .map(|ctx| ctx.karlin_params)
                .collect();
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
        
        // `ctx_idx` is the global query context index (NCBI `hsp->context` equivalent).
        // Use the pre-built `QueryContext` buffer to access the (possibly unmasked) query AA sequence.
        let q_ctx = &query_contexts_ref[h.ctx_idx];
        // NCBI uses query_nomask = sequence_nomask + query_offset (0-based).
        // Reference: blast_filter.c:1381-1382, blast_util.c:112-116.
        let q_aa_full = q_ctx.aa_seq_nomask.as_deref().unwrap_or(&q_ctx.aa_seq);
        let q_aa = &q_aa_full[1..q_aa_full.len() - 1];
        
        // Get subject frame from cache (no redundant translation)
        let s_frame_obj = &subject_frames_cache[h.s_idx as usize][h.s_f_idx];
        let s_aa_full = &s_frame_obj.aa_seq;
        let s_aa = &s_aa_full[1..s_aa_full.len() - 1];
        
        let len = h.q_aa_end.saturating_sub(h.q_aa_start);
        
        let q0 = h.q_aa_start;
        let s0 = h.s_aa_start;
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

        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        let out_hit = Hit {
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
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame: h.q_frame as i32,
            query_length: 0,
            q_idx: h.q_idx,
            s_idx: h.s_idx,
            raw_score: h.raw_score,
            gap_info: None,
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
    write_output_ncbi_order(final_hits, args.out.as_ref(), &query_ids, &subject_ids)?;
    Ok(())
}
