//! The set-up of a TBLASTX query batch that `search_query_batch` does before its subjects: the
//! search's flags and word values (`batch_flags`), the frames of the batch's queries (six, or
//! the three of the strand that `-strand` names) with their masks (`query_frames`), and the
//! batch's lookup table with its contexts and X-dropoffs (`batch_lookup`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:158-277
//! (`BlastInitialWordParametersNew`); c++/src/algo/blast/api/blast_setup_cxx.cpp:484-658
//! (`SetupQueries_OMF`) and c++/src/algo/blast/core/blast_filter.c:1261-1304
//! (`BlastSetUp_GetFilteringLocations`); c++/src/algo/blast/core/lookup_wrap.c:61-195
//! (`LookupTableWrapInit_MT`), which calls c++/src/algo/blast/core/blast_aalookup.c:227-254
//! (`BlastAaLookupTableNew`).
//!
//! This module is a child of `run_impl` for the modules it keeps private (`x_query_frames`,
//! `x_pair_par`). The functions are statements of `search_query_batch` moved here in their order;
//! what they read is a struct of references or `Copy` values or a parameter, and what they make is
//! what `search_query_batch` then reads under the same names.

use super::{
    diagnostics_enabled, generate_frames, lookup_protein_params_ungapped, trace_hsp_target,
    x_drop_raw_score, x_pair_par, x_query_frames, DiagnosticCounters, GeneticCode, KarlinParams,
    QueryContext, QueryFrame, SegMasker, Subjects, TblastxArgs, TbxScore,
};
use crate::algorithm::blastn::args::QueryStrand;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::utils::threading::SearchPool;
use anyhow::Result;
use std::time::{Duration, Instant};

/// The search's flags and word values that `search_query_batch` reads under the same names.
pub(super) struct BatchFlags {
    pub(super) extension_debug_enabled: bool,
    pub(super) timing_enabled: bool,
    /// The clock of the batch's `[TIMING] total` line.
    pub(super) t_total: Instant,
    pub(super) num_threads: usize,
    pub(super) requested_parallel: bool,
    pub(super) use_serial_scan_chunks: bool,
    pub(super) query_code: GeneticCode,
    pub(super) db_code: GeneticCode,
    pub(super) window: i32,
    pub(super) wordsize: i32,
    /// The matrix's ungapped Karlin block, which the X-dropoffs and another matrix's contexts
    /// read (`batch_lookup`).
    pub(super) ungapped_params_for_xdrop: KarlinParams,
    pub(super) diag_enabled: bool,
    pub(super) debug_cutoffs_all: bool,
    pub(super) debug_output_filter: bool,
    pub(super) debug_hsp_saving: bool,
    pub(super) diagnostics: std::sync::Arc<DiagnosticCounters>,
    pub(super) scan_debug_range: Option<(i32, i32)>,
}

/// The search's debug and timing flags, its threads, the genetic codes, the window and word size,
/// and the ungapped Karlin block of the X-dropoffs.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3202-3243
/// (`CMTArgs::x_ExtractAlgorithmOptions`) and :1042-1056
/// (`CGeneticCodeArgs::ExtractAlgorithmOptions`);
/// c++/src/algo/blast/core/blast_parameters.c:158-277 (`BlastInitialWordParametersNew`, the
/// X-dropoffs at :219-221).
pub(super) fn batch_flags<S: TbxScore>(
    sc: S,
    args: &TblastxArgs,
    parallel_pool: &SearchPool<'_>,
) -> BatchFlags {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-583
    // LOSAT-only existing diagnostics: capture configuration once per search;
    // the NCBI extension inputs and all biological decisions stay unchanged.
    let extension_debug_enabled = crate::utils::xdiag::var("LOSAT_DEBUG_EXTENSION").is_ok();
    // Optional timing breakdown (disabled by default to preserve output/parity logs)
    let timing_enabled = std::env::var_os("LOSAT_TIMING").is_some();
    let t_total = Instant::now();

    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1297
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    let num_threads = parallel_pool.threads();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    //
    // A negative `-window_size` (an imported strategy's) searches serially, in NCBI's order of
    // subjects and chunks on one diagonal table, which is NCBI's `-num_threads 1` search: its
    // table keeps its state from one call of the word finder to the next, so what NCBI's threads
    // find would depend on how they share the subjects (`word_finder::WideOffset`).
    let requested_parallel = parallel_pool.enabled() && args.window() >= 0;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003
    //
    // The scan-chunk experiment keeps the canonical two-hit state serial inside
    // each NCBI WordFinder search unit, so subject-level parallelism is disabled
    // while this gate is active.
    let use_serial_scan_chunks = std::env::var_os("LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS").is_some();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1048-1055
    // Setting each code makes its table (`AddGeneticCode`, which calls `FindGeneticCode`,
    // c++/src/algo/blast/api/blast_options_cxx.cpp:1478-1483 `SetQueryGeneticCode`, :2411-2416
    // `SetDbGeneticCode`, c++/src/algo/blast/api/blast_aux.cpp:653-660). The queries are
    // translated with theirs in `SetupQueries_OMF`
    // (c++/src/algo/blast/api/blast_setup_cxx.cpp:554-562), the `-db` subjects in
    // `BLAST_PreliminarySearchEngine` (c++/src/algo/blast/core/blast_engine.c:1458-1465), and
    // `-subject` ones also in `SetupSubjects_OMF` (blast_setup_cxx.cpp:801-809).
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    // LOSAT intentionally treats `--db-gencode` as the subject genetic code for
    // local `-s/--subject` searches. Unlike BLAST+, search/scoring and reporting
    // both honor the explicit subject code.

    // Two-hit window size (NCBI `diag->window`).
    let window = args.window();
    // Seed word length (NCBI `lookup->word_length`).
    // (`-word_size`; BLOSUM62's instance of the search is for 3, `TbxScore::word`)
    let wordsize: i32 = sc.word() as i32;

    // NCBI BLAST computes x_dropoff per-context using kbp[context]->Lambda:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes kbp_std per context and applies check_ideal:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2796-2797
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
    // We still maintain per-context structure for parity.
    //
    // x_dropoff_per_context is populated after build_ncbi_lookup() creates the contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(sc.matrix());

    let diag_enabled = diagnostics_enabled();
    let debug_cutoffs_all = crate::utils::xdiag::var_os("LOSAT_DEBUG_CUTOFFS_ALL").is_some();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3265-3274
    let debug_output_filter = std::env::var_os("LOSAT_DEBUG_OUTPUT_FILTER").is_some();
    let debug_hsp_saving = crate::utils::xdiag::var_os("LOSAT_DEBUG_HSP_SAVING").is_some();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    // Debug: optional scan output dump around a specific subject offset.
    // The offset_pairs correspond to s_BlastAaScanSubject copy loop.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-115
    let scan_debug_center = crate::utils::xdiag::var("LOSAT_DEBUG_SCAN_SOFF")
        .ok()
        .and_then(|v| v.parse::<i32>().ok());
    let scan_debug_window = crate::utils::xdiag::var("LOSAT_DEBUG_SCAN_WINDOW")
        .ok()
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    let scan_debug_range =
        scan_debug_center.map(|center| (center - scan_debug_window, center + scan_debug_window));
    if let Some((lo, hi)) = scan_debug_range {
        eprintln!(
            "[DEBUG SCAN_OFF] enabled s_off_range=[{},{}] (center={} window={})",
            lo,
            hi,
            scan_debug_center.unwrap_or(0),
            scan_debug_window
        );
    }
    BatchFlags {
        extension_debug_enabled,
        timing_enabled,
        t_total,
        num_threads,
        requested_parallel,
        use_serial_scan_chunks,
        query_code,
        db_code,
        window,
        wordsize,
        ungapped_params_for_xdrop,
        diag_enabled,
        debug_cutoffs_all,
        debug_output_filter,
        debug_hsp_saving,
        diagnostics,
        scan_debug_range,
    }
}

/// The six frames of the batch's queries, which `search_query_batch` reads under the same names.
pub(super) struct QueryFrames {
    pub(super) query_frames: Vec<Vec<QueryFrame>>,
    pub(super) t_read_queries: Duration,
}

/// The six frames of each query of the batch, translated with the query's genetic code and
/// masked with SEG and, with `-lcase_masking`, the lower-case runs.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:484-658
/// (`SetupQueries_OMF`, which translates each context of a query with `BLAST_GetTranslation`
/// at :576-588) and c++/src/algo/blast/core/blast_filter.c:1261-1304
/// (`BlastSetUp_GetFilteringLocations`).
pub(super) fn query_frames(
    args: &TblastxArgs,
    queries_raw: &[FastaRecord],
    query_code: &GeneticCode,
    parallel_pool: &SearchPool<'_>,
) -> Result<QueryFrames> {
    let t_phase_read_queries = Instant::now();
    // NCBI tblastx low-complexity filtering uses SEG on translated protein sequences.
    // No nucleotide-level DUST masking is applied.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:572-580

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:576-588 (the six
    // contexts of a tblastx query, each translated with the query's genetic code)
    // Diagnostic (LOSAT_TIMING): the queries' six-frame translation for the common line
    // (utils/xtiming.rs). No NCBI counterpart.
    let t_translate = crate::utils::xtiming::clock();
    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), query_code))
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:200-215 (with
    // `-strand plus` contexts 3-5 of every query, with `minus` contexts 0-2, have length 0;
    // `s_QueryInfo_SetContext` at :65-87 gives such a context no letters in the buffer and no
    // separator, and marks it invalid)
    // The frames of the strand that is not searched are dropped here, once per batch: the batch's
    // buffer, lookup table, Karlin blocks, cutoffs and search spaces then hold the searched
    // frames only, as NCBI's (a frame of no letters is already absent, `generate_frames`).
    let (plus, minus) = QueryStrand::from_option(&args.strand).searched();
    if !(plus && minus) {
        for frames in &mut query_frames {
            frames.retain(|frame| if frame.frame > 0 { plus } else { minus });
        }
    }
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Translate, t_translate);

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:404-406
    if args.lcase_masking {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1239-1256 (see
        // `mask_lcase`): the lower-case runs join the SEG intervals of each frame.
        let seg = args
            .seg_spec()?
            .params()
            .map(|params| SegMasker::new(params.window, params.locut, params.hicut));
        x_query_frames::mask_lcase(
            &mut query_frames,
            queries_raw,
            seg.as_ref(),
            args.soft_masking(),
            parallel_pool,
        );
    } else if let Some(params) = args.seg_spec()?.params() {
        // Diagnostic (LOSAT_TIMING): the queries' SEG for the common line (utils/xtiming.rs). No
        // NCBI counterpart.
        let t_seg = crate::utils::xtiming::clock();
        let seg = SegMasker::new(params.window, params.locut, params.hicut);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
        // Dispatch point (LOSAT_X_TBXQFRAMES, strict set): with two or more threads the frames
        // are masked on the search pool, one task per frame; otherwise, and with the switch off,
        // the frames are masked one after the other as before (tblastx/x_query_frames.rs, which
        // holds the per-frame bodies `mask_frame` and `mask_frame_soft`).
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
        // `-soft_masking true` keeps the intervals for the lookup table and the frames' letters.
        if args.soft_masking() {
            x_query_frames::mask_soft(&mut query_frames, &seg, parallel_pool);
        } else {
            x_query_frames::mask(&mut query_frames, &seg, parallel_pool);
        }
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
    }

    let t_read_queries = t_phase_read_queries.elapsed();
    Ok(QueryFrames {
        query_frames,
        t_read_queries,
    })
}

/// What `batch_lookup` reads.
pub(super) struct LookupInputs<'a, S> {
    pub(super) args: &'a TblastxArgs,
    pub(super) subjects: Subjects<'a>,
    pub(super) query_frames: &'a [Vec<QueryFrame>],
    pub(super) window: i32,
    pub(super) scan_debug_range: Option<(i32, i32)>,
    pub(super) extension_debug_enabled: bool,
    pub(super) diag_enabled: bool,
    pub(super) debug_hsp_saving: bool,
    pub(super) requested_parallel: bool,
    pub(super) use_serial_scan_chunks: bool,
    pub(super) parallel_pool: &'a SearchPool<'a>,
    pub(super) sc: S,
    pub(super) ungapped_params_for_xdrop: KarlinParams,
}

/// The batch's lookup table and contexts, which `search_query_batch` reads under the same names.
pub(super) struct BatchLookup {
    /// `-lcase_masking` with a lower-case letter in a subject.
    pub(super) x_soft_subjects: bool,
    /// PAIRLUT: the groups built at dispatch point 0, for `prepare` (dispatch point 1).
    pub(super) x_pair_lut_batch: Option<x_pair_par::XLutBatch>,
    pub(super) x_pair_lut_is_stub: bool,
    pub(super) lookup: BlastAaLookupTable,
    pub(super) contexts: Vec<QueryContext>,
    pub(super) t_build_lookup: Duration,
    /// One X-dropoff per context of `contexts`.
    pub(super) x_dropoff_per_context: Vec<i32>,
}

/// The batch's lookup table (or PAIRLUT's stub with its unit tables), the contexts of the queries'
/// frames with their Karlin blocks, and the X-dropoff of each context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:61-195
/// (`LookupTableWrapInit_MT`), c++/src/algo/blast/core/blast_aalookup.c:227-254
/// (`BlastAaLookupTableNew`) and c++/src/algo/blast/core/blast_parameters.c:158-277
/// (`BlastInitialWordParametersNew`, the X-dropoffs at :219-221).
pub(super) fn batch_lookup<S: TbxScore>(inputs: &LookupInputs<'_, S>) -> BatchLookup {
    let &LookupInputs {
        args,
        subjects,
        query_frames,
        window,
        scan_debug_range,
        extension_debug_enabled,
        diag_enabled,
        debug_hsp_saving,
        requested_parallel,
        use_serial_scan_chunks,
        parallel_pool,
        sc,
        ungapped_params_for_xdrop,
    } = inputs;
    let t_phase_build_lookup = Instant::now();
    // Note: karlin_params argument is now unused - computed per context in build_ncbi_lookup()
    // We still pass it for x_dropoff calculation (which uses ideal params for all contexts in tblastx)
    // NCBI lookup build always precomputes neighbors (no lazy mode).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-543
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456 (exact match indexing before neighbor expansion)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:245
    // The x86-64 conversion: `INT_MIN` for a threshold beyond `Int4` (every word).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:813-826
    // `-lcase_masking` with a lower-case letter in a subject: only the reference frame loop scans
    // ranges (`x_soft_subject_calls`), so the batch takes it as a debugging run does. A `-db`
    // subject has no case.
    // A masked `-db` search (`-db_soft_mask`, `-db_hard_mask`; series NP, M4b PR 2) takes it too:
    // every subject's frames have their own ranges or chunks (`x_db_mask_calls`).
    let x_soft_subjects = (args.lcase_masking
        && matches!(subjects, Subjects::Fasta(records)
            if records.iter().any(|record| record.seq().iter().any(u8::is_ascii_lowercase))))
        || matches!(subjects, Subjects::Stream(stream) if stream.mask().is_some());
    // EXPERIMENT (PAIRLUT, graduated on 2026-10-22), dispatch point 0 (per batch, at the
    // lookup build): when PAIRPAR will run the seed stage of every subject of this batch
    // as units, the units' tables (one per query-context group) are built directly from the
    // query, and `lookup` becomes a stub with the full table's header and no entries
    // (x_pair_par.rs `lut_build`; lookup/x_pair_lut.rs for why the tables are the same). The
    // inputs are the ones dispatch point 1 (`seed_tables::scan_tables`, after this call) gives
    // `prepare`, by the same expressions; `prepare` decides again and stops the run if it
    // disagrees. When the unit path will not run, the reference builder below runs.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
    // PAIRLUT builds the group tables with the search's matrix (SC7b; `XLutEarly::sc`), and the
    // batch's full table of TBXLUT holds the scores and words of the search's matrix and word size.
    let mut x_pair_lut_batch = {
        x_pair_par::lut_build(x_pair_par::XLutEarly {
            query_frames,
            threshold: crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
            window,
            subjects,
            reference_loop: trace_hsp_target().is_some()
                || scan_debug_range.is_some()
                || extension_debug_enabled
                || diag_enabled
                || debug_hsp_saving
                || x_soft_subjects
                // The one-hit word finder (`-window_size 0` or less) runs in the reference frame
                // loop only (`scan_chunk_one_hit`).
                || window <= 0,
            requested_parallel,
            use_serial_scan_chunks,
            pool: parallel_pool,
            sc,
        })
    };
    let x_pair_lut_stub = x_pair_lut_batch
        .as_mut()
        .and_then(x_pair_par::XLutBatch::take_stub);
    // True when `lookup` is the stub: the reference frame loop must then run zero times
    // (asserted at dispatch point 2).
    let x_pair_lut_is_stub = x_pair_lut_stub.is_some();
    // Dispatch point (LOSAT_X_TBXLUT, on unless `LOSAT_X_TBXLUT=0`): the full table of the batch
    // from the chain builder (`build_ncbi_lookup_for_profile`) or filled in place from per-cell
    // counts; the same table.
    let (lookup, contexts) = match x_pair_lut_stub {
        Some(stub) if S::BLOSUM62 => {
            crate::algorithm::tblastx::lookup::x_tbxlut_reach_stub();
            stub
        }
        Some((lookup, mut contexts)) => {
            crate::algorithm::tblastx::lookup::x_tbxlut_reach_stub();
            // Another matrix's contexts (SC7b): the alpha and beta of its ungapped row, as in the
            // reference branch below (`blast_setup.c:812-822`).
            for ctx in contexts.iter_mut().filter(|ctx| ctx.is_valid) {
                ctx.karlin_params.alpha = ungapped_params_for_xdrop.alpha;
                ctx.karlin_params.beta = ungapped_params_for_xdrop.beta;
            }
            (lookup, contexts)
        }
        None => {
            // The table of the search's matrix and word size (`sc`), built without chains (another
            // matrix or word size than BLOSUM62's 3 letters: TBLASTN's builders; item X2e of line RG).
            let (lookup, mut contexts) = crate::algorithm::tblastx::lookup::x_tbx_build_ncbi_lookup(
                query_frames,
                crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
                &ungapped_params_for_xdrop, // Used for x_dropoff calculation only
                true,
                sc.matrix(),
                sc.word(),
            );
            if !S::BLOSUM62 {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:812-822 (the alpha and
                // beta of the length adjustment: the matrix's ungapped row,
                // `BLAST_GetAlphaBeta`, for the Karlin block of every context)
                // (`compute_karlin_params_ungapped` gives every block BLOSUM62's.)
                for ctx in contexts.iter_mut().filter(|ctx| ctx.is_valid) {
                    ctx.karlin_params.alpha = ungapped_params_for_xdrop.alpha;
                    ctx.karlin_params.beta = ungapped_params_for_xdrop.beta;
                }
            }
            (lookup, contexts)
        }
    };
    let t_build_lookup = t_phase_build_lookup.elapsed();
    // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs). No
    // NCBI counterpart.
    crate::utils::xtiming::add_duration(crate::utils::xtiming::Phase::BuildLookup, t_build_lookup);

    // Initial per-context ungapped X-dropoff (NCBI `word_params->cutoffs[context].x_dropoff_init`).
    // Compute per-context x_dropoff using `kbp[context]->Lambda`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes `kbp_std` per context and applies `check_ideal`,
    // copying `kbp_ideal` into `kbp` via `Blast_KarlinBlkCopy` whenever `kbp->Lambda` is at least `kbp_ideal->Lambda`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
    //
    // NCBI BLAST dynamic x_dropoff (c++/src/algo/blast/core/blast_parameters.c:380-383):
    // when initialized to zero, `curr_cutoffs->x_dropoff` is set dynamically to the cutoff score.
    //
    // For TBLASTX, x_dropoff_init is non-zero, so this dynamic update is not triggered.
    // However, we store x_dropoff_init here and apply the logic during subject processing
    // where cutoff_score is available.
    // The X-drop in bits: `-xdrop_ungap` or TBLASTX's default (`TblastxArgs::x_drop_ungapped_bits`).
    let x_drop_ungapped_bits = args.x_drop_ungapped_bits();
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(|ctx| x_drop_raw_score(x_drop_ungapped_bits, &ctx.karlin_params, 1.0))
        .collect();
    BatchLookup {
        x_soft_subjects,
        x_pair_lut_batch,
        x_pair_lut_is_stub,
        lookup,
        contexts,
        t_build_lookup,
        x_dropoff_per_context,
    }
}
