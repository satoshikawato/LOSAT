//! Main run() function for backbone mode
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use super::*;
use std::sync::Arc;

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
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes kbp_std per context and applies check_ideal:
    //   if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
    //      Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
    // We still maintain per-context structure for parity.
    //
    // x_dropoff_per_context is populated after build_ncbi_lookup() creates the contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    let diag_enabled = diagnostics_enabled();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    // Debug: optional scan output dump around a specific subject offset.
    // The offset_pairs correspond to s_BlastAaScanSubject copy loop.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115
    let scan_debug_center = std::env::var("LOSAT_DEBUG_SCAN_SOFF")
        .ok()
        .and_then(|v| v.parse::<i32>().ok());
    let scan_debug_window = std::env::var("LOSAT_DEBUG_SCAN_WINDOW")
        .ok()
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    let scan_debug_range = scan_debug_center.map(|center| (center - scan_debug_window, center + scan_debug_window));
    if let Some((lo, hi)) = scan_debug_range {
        eprintln!(
            "[DEBUG SCAN_OFF] enabled s_off_range=[{},{}] (center={} window={})",
            lo,
            hi,
            scan_debug_center.unwrap_or(0),
            scan_debug_window
        );
    }

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
    // NCBI lookup build always precomputes neighbors (no lazy mode).
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-543
    let (lookup, contexts) = build_ncbi_lookup(
        &query_frames,
        args.threshold,
        args.include_stop_seeds,
        args.ncbi_stop_stop_score,
        &ungapped_params_for_xdrop, // Used for x_dropoff calculation only
    );
    t_build_lookup = t_phase_build_lookup.elapsed();

    // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
    // Compute per-context x_dropoff using kbp[context]->Lambda.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes kbp_std per context and applies check_ideal:
    //   if (check_ideal && kbp->Lambda >= kbp_ideal->Lambda) Blast_KarlinBlkCopy(kbp, kbp_ideal);
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
    //
    // NCBI BLAST dynamic x_dropoff (ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:380-383):
    //   if (curr_cutoffs->x_dropoff_init == 0)
    //      curr_cutoffs->x_dropoff = new_cutoff;  // x_dropoff = cutoff_score
    //   else
    //      curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
    //
    // For TBLASTX, x_dropoff_init is non-zero, so this dynamic update is not triggered.
    // However, we store x_dropoff_init here and apply the logic during subject processing
    // where cutoff_score is available.
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(|ctx| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ctx.karlin_params, 1.0))
        .collect();

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
    // NCBI reference (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:768):
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
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
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
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    // } BlastHSPList;
    // ```
    let query_ids_out = query_ids.clone();
    let subject_ids_out = subject_ids.clone();

    let writer = std::thread::spawn(move || -> Result<()> {
        let mut all: Vec<Hit> = Vec::new();
        while let Ok(h) = rx.recv() {
            all.extend(h);
        }
        all.retain(|h| h.e_value <= evalue_threshold);
        // NCBI-style output ordering: query (input order) → subject (best_evalue/score/oid) → HSP (score/coords)
        // Reference: BLAST_LinkHsps() + s_EvalueCompareHSPLists() + ScoreCompareHSPs()
        write_output_ncbi_order(all, out_path.as_ref(), &query_ids_out, &subject_ids_out)?;
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
    // `query_length` here must match BLAST_SequenceBlk->length, computed as
    // last_context.query_offset + last_context.query_length (excludes the final trailing NULLB).
    // References: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-315, 378-381
    // NCBI BLAST diag array sizing:
    // diag_array_length = next_power_of_2(qlen + window_size)
    // For tblastx, qlen = total concatenated query buffer (all 6 frames)
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    let query_length: i32 = contexts
        .last()
        .map(|c| c.frame_base + c.aa_len as i32)
        .unwrap_or(0);
    let mut diag_array_size: i32 = 1;
    while diag_array_size < (query_length + window) {
        diag_array_size <<= 1;
    }
    let diag_mask: i32 = diag_array_size - 1;
    if trace_hsp_target().is_some() {
        // NCBI: diag_array_length = next_power_of_2(qlen + window_size).
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
        eprintln!(
            "[TRACE_HSP] diag_table query_length={} window={} diag_array_size={} diag_mask={}",
            query_length,
            window,
            diag_array_size,
            diag_mask
        );
    }

    // [C] array_size for offset_pairs
    // NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain
    // Reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h + lookup_wrap.c
    const OFFSET_ARRAY_SIZE: i32 = 4096;
    let offset_array_size: i32 = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);

    let lookup_ref = &lookup;
    let contexts_ref = &contexts;
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

            let s_len = s_rec.seq().len();

            // [C] BlastOffsetPair *offset_pairs
            let offset_pairs = &mut st.offset_pairs;

            // [C] DiagStruct *diag_array = diag->hit_level_array;
            let diag_array = &mut st.diag_array;

            // [C] diag_offset = diag->offset;  (reset to window per-subject)
            let mut diag_offset: i32 = st.diag_offset;

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
            let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; contexts_ref.len()];
            // NCBI word_params->cutoff_score_min = min of cutoffs across all contexts
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:401-403
            let mut cutoff_score_min = i32::MAX;
            
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
            let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
            let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
            
            for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
                // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
                let ctx_params = &ctx.karlin_params;

                // NCBI: query_length = query_info->contexts[context].query_length
                let query_len_aa = ctx.aa_len as i64;

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

                let subject_full = &s_frame.aa_seq;
                // NCBI subject->sequence points past the leading NULLB sentinel, so offsets
                // are 0-based from the first residue; see blast_engine.c:811-812 and
                // blast_util.c:112-116.
                let subject = &subject_full[1..subject_full.len() - 1];
                let s_aa_len = s_frame.aa_len;
                if subject.len() < wordsize as usize {
                    // NCBI still advances the diagonal table offset even when no hits can be found.
                    // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:445-446,
                    // ncbi-blast/c++/src/algo/blast/core/blast_extend.c:167-173
                    if diag_offset >= i32::MAX / 4 {
                        diag_offset = window;
                        // s_BlastDiagClear(): clear all diagonal state when offset risks overflow.
                        // NCBI sets last_hit = -window (blast_extend.c:101-103).
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:101-103
                        for d in diag_array.iter_mut() {
                            *d = DiagStruct::clear(window);
                        }
                    } else {
                        diag_offset += s_aa_len as i32 + window;
                    }
                    continue;
                }

                // NCBI subject seq_ranges are used by s_DetermineScanningOffsets (masksubj.inl).
                // With no subject masking, the range is [0, subject->length].
                // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:509-511,
                // ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
                let seq_ranges: [(i32, i32); 1] = [(0, s_aa_len as i32)];
                // [C] scan_range[0] = 0;
                // [C] scan_range[1] = subject->seq_ranges[0].left;
                // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];

                // [C] while (scan_range[1] <= scan_range[2])
                while scan_range[1] <= scan_range[2] {
                    let prev_scan_left = scan_range[1];
                    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                    let t0 = if timing_enabled { Some(Instant::now()) } else { None };
                    let hits = s_blast_aa_scan_subject(
                        lookup_ref,
                        subject,
                        &seq_ranges,
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
                            // NCBI BlastOffsetPair uses Uint4 offsets.
                            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                            let mut seen: HashSet<(u32, u32)> = HashSet::with_capacity(hits as usize);
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

                        // Debug: dump scan output around a target subject offset.
                        // The scan output is produced by s_BlastAaScanSubject.
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-123
                        if let Some((lo, hi)) = scan_debug_range {
                            // NCBI offsets are Uint4; cast for debug range comparison only.
                            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                            let subject_offset_i32 = subject_offset as i32;
                            if subject_offset_i32 >= lo && subject_offset_i32 <= hi {
                                eprintln!(
                                    "[DEBUG SCAN_OFF] s_f_idx={} s_off={} q_off={} scan_range=[{},{}] diag_offset={}",
                                    s_f_idx,
                                    subject_offset,
                                    query_offset,
                                    prev_scan_left,
                                    scan_range[1],
                                    diag_offset
                                );
                            }
                        }

                        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                        // NCBI uses Uint4 for offsets; apply unsigned wrapping.
                        // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                        //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:534
                        let diag_coord =
                            (query_offset.wrapping_sub(subject_offset) & (diag_mask as u32)) as usize;
                        
                        // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
                        let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

                        // [C] if (diag_array[diag_coord].flag)
                        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:536-553
                        if diag_entry.flag() != 0 {
                            // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                            let subject_plus_offset =
                                subject_offset.wrapping_add(diag_offset as u32);
                            if subject_plus_offset < diag_entry.last_hit() as u32 {
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
                            diag_entry.set_last_hit(subject_plus_offset as i32);
                            diag_entry.set_flag(0);
                            // Track flag reset (hit after previous extension zone)
                            if diag_enabled {
                                diagnostics
                                    .base
                                    .seeds_flag_reset
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                            }
                        }
                        // [C] else
                        else {
                            // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                            let last_hit = diag_entry.last_hit() - diag_offset;
                            // [C] diff = subject_offset - last_hit;
                            // NCBI uses Uint4 for subject_offset; compute with unsigned wrap.
                            // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                            //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:559-560
                            let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

                            // [C] if (diff >= window)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:562-569
                            if diff >= window {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_too_far
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                diag_entry.set_last_hit(
                                    subject_offset.wrapping_add(diag_offset as u32) as i32,
                                );
                                continue;
                            }

                            // [C] if (diff < wordsize)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:573-580
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
                            // NCBI passes Uint4 query_offset into BSearchContextInfo (Int4).
                            // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                            //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:590
                            let ctx_idx = lookup_ref.get_context_idx(query_offset as i32);
                            let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
                            let q_raw =
                                query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
                            // NCBI uses masked sequence for extension; query->sequence is
                            // sequence_start + 1, so offsets are 0-based in that buffer.
                            // Reference: blast_query_info.c:311-315, blast_util.c:112-116.
                            let query_full = &ctx.aa_seq;
                            let query = &query_full[1..query_full.len() - 1];

                            // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:592-606
                            let q_minus_diff = query_offset.wrapping_sub(diff as u32);
                            if q_minus_diff < ctx.frame_base as u32 {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_ctx_boundary_fail
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                diag_entry.set_last_hit(
                                    subject_offset.wrapping_add(diag_offset as u32) as i32,
                                );
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
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:636-648
                            // if (right_extend) {
                            //     diag_array[diag_coord].flag = 1;
                            //     diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                            // } else {
                            //     diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            // }
                            if right_extend {
                                // "If an extension to the right happened, reset the last hit
                                //  so that future hits to this diagonal must start over."
                                diag_entry.set_flag(1);
                                diag_entry.set_last_hit(
                                    s_last_off - (wordsize - 1) + diag_offset,
                                );
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .mask_updates
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                            } else {
                                // "Otherwise, make the present hit into the previous hit for this diagonal"
                                diag_entry.set_last_hit(
                                    subject_offset.wrapping_add(diag_offset as u32) as i32,
                                );
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
                                    // NCBI offsets are 0-based in query/subject->sequence buffers.
                                    // Reference: blast_gapalign.c:4756-4768, blast_aascan.c:110-113.
                                    let q_aa_start = hsp_q_u as usize;
                                    let q_aa_end = hsp_qe_u as usize;
                                    let s_aa_start = hsp_s_u as usize;
                                    let s_aa_end = _hsp_se_u as usize;
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
                                        // NCBI two-hit gating checks (diff/window/wordsize/context).
                                        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:531-606
                                        // NCBI uses Uint4 offsets; apply unsigned wrap here as well.
                                        // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                        let q_minus_diff = query_offset.wrapping_sub(diff as u32);
                                        eprintln!(
                                            "[TRACE_HSP] two_hit_pass diff={} window={} wordsize={} diff>=window={} diff<wordsize={} q_minus_diff={} ctx_frame_base={} q_minus_diff<base={} diag_offset={} diag_mask={} diag_array_size={}",
                                            diff,
                                            window,
                                            wordsize,
                                            diff >= window,
                                            diff < wordsize,
                                            q_minus_diff,
                                            ctx.frame_base,
                                            q_minus_diff < ctx.frame_base as u32,
                                            diag_offset,
                                            diag_mask,
                                            diag_array_size
                                        );
                                    }
                                }
                                // NCBI: BlastSaveInitHsp equivalent
                                // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                                // Store HSP with absolute coordinates (before coordinate conversion)
                                //
                                // hsp_q is frame-relative coordinate in query->sequence (0-based),
                                // frame_base is the context query_offset in the concatenated buffer.
                                // NCBI: ungapped_data->q_start is absolute query offset.
                                // Reference: blast_gapalign.c:4756-4768, blast_query_info.c:311-315.
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

                // NCBI parity: s_BlastFindSmallestLambda selects smallest lambda across contexts.
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:92-112
                // Per-context kbp_std is computed from query composition (with check_ideal).
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
                let context_params: Vec<KarlinParams> = contexts_ref
                    .iter()
                    .map(|ctx| ctx.karlin_params)
                    .collect();
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
                    // NCBI uses query_nomask = query_blk->sequence_nomask + query_offset,
                    // with sequence_nomask pointing past the leading NULLB.
                    // Reference: blast_filter.c:1381-1382, blast_util.c:112-116.
                    let q0 = h.q_aa_start;
                    let s0 = h.s_aa_start;
                    let q_seq_nomask_full: &[u8] = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
                    let q_seq_nomask = &q_seq_nomask_full[1..q_seq_nomask_full.len() - 1];
                    let s_seq = &s_frame.aa_seq[1..s_frame.aa_seq.len() - 1];
                    // Use SIMD-optimized identity calculation (already implemented in reevaluate.rs)
                    let t_identity = if timing_enabled { Some(Instant::now()) } else { None };
                    let matches = get_num_identities_and_positives_ungapped(
                        q_seq_nomask,
                        s_seq,
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
                        mismatch,
                        gapopen: 0,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        e_value: h.e_value,
                        bit_score: bit,
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
                        // ```c
                        // if (hsp->query.frame != hsp->subject.frame) {
                        //    *q_end = query_length - hsp->query.offset;
                        //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
                        // }
                        // ```
                        query_frame: ctx.frame as i32,
                        query_length: 0,
                        q_idx: ctx.q_idx,
                        s_idx: h.s_idx,
                        raw_score: h.raw_score,
                        gap_info: None,
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
