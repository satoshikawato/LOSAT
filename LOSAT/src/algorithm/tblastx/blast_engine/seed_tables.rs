//! The seed-stage tables of a TBLASTX query batch, made before its subjects are searched: the
//! sizes of the diagonal table and of the offset arrays, LOSAT_X_TBXSCANBF's choice, the batch's
//! unit path (PAIRPAR, dispatch point 1) and its batch on the board (`scan_tables`), then the
//! holders of the seed release, the size of the workers' diagonal tables and the extensions'
//! sentinel walks (`seed_batch`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:970-1069
//! (`s_BlastSetUpAuxStructures`), which sizes the offset arrays with
//! c++/src/algo/blast/core/lookup_wrap.c:255-289 (`GetOffsetArraySize`) and makes the diagonal
//! table with c++/src/algo/blast/core/blast_extend.c:110-159 (`BlastExtendWordNew`) and :41-67
//! (`s_BlastDiagTableNew`).
//!
//! This module is a child of `run_impl` for the modules it keeps private (`x_pair_par`,
//! `x_seed_release`, `x_subject_next`). Its functions are statements of `search_query_batch`
//! moved here in their order. `search_query_batch` makes the guard that takes the batch off the
//! board (`XEndBatch`) between the two calls: it borrows the first one's result. What they read is
//! a struct of references or `Copy` values and, by value, the tables they take over; what they
//! make is what `search_query_batch` then reads under the same names.

use super::{
    trace_hsp_target, x_pair_par, x_seed_release, x_subject_next, x_tbx_sent_contexts,
    QueryContext, Subjects, TblastxArgs, TbxScore,
};
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::algorithm::tblastx::x_work_board::XBoard;
use crate::utils::threading::SearchPool;
use anyhow::Result;
use std::sync::Arc;

/// What `scan_tables` reads.
pub(super) struct ScanInputs<'a, S> {
    pub(super) contexts: &'a [QueryContext],
    pub(super) lookup: &'a BlastAaLookupTable,
    pub(super) window: i32,
    pub(super) subjects: Subjects<'a>,
    pub(super) scan_debug_range: Option<(i32, i32)>,
    pub(super) extension_debug_enabled: bool,
    pub(super) diag_enabled: bool,
    pub(super) debug_hsp_saving: bool,
    pub(super) x_soft_subjects: bool,
    pub(super) use_serial_scan_chunks: bool,
    pub(super) use_parallel_chunks: bool,
    pub(super) use_parallel: bool,
    pub(super) parallel_pool: &'a SearchPool<'a>,
    pub(super) args: &'a TblastxArgs,
    pub(super) sc: S,
}

/// The sizes of the batch's diagonal table and offset arrays, its unit path and its batch on the
/// board, which `search_query_batch` reads under the same names.
pub(super) struct ScanTables<'b, S: TbxScore> {
    pub(super) diag_array_size: i32,
    pub(super) diag_mask: i32,
    pub(super) offset_array_size: i32,
    pub(super) offset_pairs_len: usize,
    pub(super) x_scan_bf: bool,
    /// PAIRPAR: the batch's unit path, `None` when the reference frame loop takes every subject.
    pub(super) x_pair_par_batch: Option<x_pair_par::XPairPar<S>>,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): the search's board and the batch on it, `None` when the
    /// batch does not use the board.
    pub(super) x_board_batch: Option<(
        &'b XBoard,
        Arc<crate::algorithm::tblastx::x_work_board::XBatch>,
    )>,
}

/// The sizes of the batch's diagonal table and offset arrays, LOSAT_X_TBXSCANBF's choice, the
/// unit path of PAIRPAR (dispatch point 1, which takes `x_pair_lut_batch`, the groups PAIRLUT
/// built), and the batch on the search's board (`x_board`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:41-67
/// (`s_BlastDiagTableNew`) and c++/src/algo/blast/core/lookup_wrap.c:255-289
/// (`GetOffsetArraySize`), which c++/src/algo/blast/core/blast_engine.c:970-1069
/// (`s_BlastSetUpAuxStructures`) calls at :992 and, through `BlastExtendWordNew`, at :1002-1004.
pub(super) fn scan_tables<'b, S: TbxScore>(
    inputs: &ScanInputs<'_, S>,
    x_pair_lut_batch: Option<x_pair_par::XLutBatch>,
    x_board: Option<&'b XBoard>,
) -> Result<ScanTables<'b, S>> {
    let &ScanInputs {
        contexts,
        lookup,
        window,
        subjects,
        scan_debug_range,
        extension_debug_enabled,
        diag_enabled,
        debug_hsp_saving,
        x_soft_subjects,
        use_serial_scan_chunks,
        use_parallel_chunks,
        use_parallel,
        parallel_pool,
        args,
        sc,
    } = inputs;
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
    crate::blastinput::app::check_diag_table_window(query_length, window, "TBLASTX")?;
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
            query_length, window, diag_array_size, diag_mask
        );
    }

    // [C] array_size for offset_pairs
    // NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain
    // Reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h + lookup_wrap.c
    const OFFSET_ARRAY_SIZE: i32 = 4096;
    let offset_array_size: i32 = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);
    // EXPERIMENT (LOSAT_X_TBXSCANBF): with the switch, the offset arrays have X_SCAN_BF_SLACK
    // entries after `offset_array_size`.
    let offset_pairs_len =
        crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size);
    // EXPERIMENT (LOSAT_X_TBXSCANBF): whether the reference frame loop scans this batch's table with
    // the branch-free loop (decided once per table).
    let x_scan_bf = crate::algorithm::tblastx::x_scan_bf::x_scan_bf_choice(lookup);

    // EXPERIMENT (PAIRPAR), dispatch point 1 (per batch): the gates, the exactness guards (the
    // first one only without LOSAT_X_TBXBAND) and the cost rule, then one lookup table per query
    // context (each backbone chain filtered to the context, order kept). None (a gate, a guard
    // or the cost rule) keeps the reference loop for every subject of the batch.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-63
    // ```c
    // while (diag_array_length < (qlen+window_size))
    // {
    //         diag_array_length = diag_array_length << 1;
    // }
    // ...
    // diag_table->diag_mask = diag_array_length-1;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-343
    // ```c
    // if (lookup->thin_backbone[i] ) {
    //     Int4 * dest = NULL;
    //     /* set the corresponding bit in the pv_array */
    //     PV_SET(pv, i, PV_ARRAY_BTS);
    //     bbc[i].num_used = lookup->thin_backbone[i][1];
    // ```
    // The guards compare the subject lengths (and, without LOSAT_X_TBXBAND, the context lengths)
    // with the slack of this table (size above minus query_length); the per-context tables use the
    // layout above.
    let x_pair_par_batch = {
        x_pair_par::XPairPar::prepare(x_pair_par::XBatch {
            lookup,
            contexts,
            query_length,
            diag_array_size,
            subjects,
            debugging: trace_hsp_target().is_some()
                || scan_debug_range.is_some()
                || extension_debug_enabled
                || diag_enabled
                || debug_hsp_saving
                || x_soft_subjects
                // The one-hit word finder (`-window_size 0`) runs in the reference frame loop
                // only (`scan_chunk_one_hit`).
                || window == 0,
            chunk_switches: use_serial_scan_chunks || use_parallel_chunks,
            // PAIRLUT: the groups built at dispatch point 0.
            lut: x_pair_lut_batch,
            // LOSAT_X_PAIRPARGATE: the thread count and threshold dispatch point 0 decided with.
            threads: parallel_pool.threads(),
            threshold: crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
            sc,
        })
    };
    // The reference frame loop scans the batch's table (LOSAT_X_TBXSCANBF's choice above).
    if x_pair_par_batch.is_none() {
        crate::algorithm::tblastx::x_scan_bf::x_reach_scan_bf(x_scan_bf);
    }
    // EXPERIMENT (LOSAT_X_TBXBOARD): every batch of a search with a board (2+ threads) uses it: its
    // units (the unit path), link groups, frame merges and chain sorts go on the search's board
    // (tblastx/x_work_board.rs). Several subjects run in TBXNEXT's loop (the other subject loop
    // waits in the pool); one subject runs on this thread and its units bring their own helpers.
    let x_board_batch = x_board
        .filter(|_| {
            use_parallel
                && !use_serial_scan_chunks
                && (subjects.count() == 1 || x_subject_next::on())
        })
        .map(|board| (board, board.batch(subjects.count())));
    Ok(ScanTables {
        diag_array_size,
        diag_mask,
        offset_array_size,
        offset_pairs_len,
        x_scan_bf,
        x_pair_par_batch,
        x_board_batch,
    })
}

/// What `seed_batch` reads.
pub(super) struct SeedInputs<'a, S> {
    pub(super) x_board_batch: &'a Option<(
        &'a XBoard,
        Arc<crate::algorithm::tblastx::x_work_board::XBatch>,
    )>,
    pub(super) subjects: Subjects<'a>,
    pub(super) diag_array_size: i32,
    pub(super) contexts: &'a [QueryContext],
    pub(super) sc: S,
    pub(super) extension_debug_enabled: bool,
}

/// The holders of the seed release, which own the batch's tables from here on, and what the
/// subject loop and `process_subject` read of the batch, which `search_query_batch` reads under
/// the same names.
pub(super) struct SeedBatch<S: TbxScore> {
    /// EXPERIMENT (LOSAT_X_TBXBOARD): whether the batch puts units on the board.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    pub(super) x_board_unit_path: bool,
    pub(super) x_seed_batch: x_seed_release::XSeedBatch<S>,
    pub(super) x_worker_diag_cells: usize,
    pub(super) tbx_sent: bool,
}

/// The board's slots of the batch's subjects, closed when the batch has no unit path; the holders
/// of the seed release, which take `lookup` and the unit path; the size the workers' diagonal
/// tables are made with; and whether the extensions take the walks the sentinels stop
/// (LOSAT_X_TBXUNGAPSENT).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1004 (the diagonal table
/// that `s_BlastSetUpAuxStructures` makes with `BlastExtendWordNew`) and :1604-1605 (where
/// `BLAST_PreliminarySearchEngine` frees it with the word parameters).
pub(super) fn seed_batch<S: TbxScore>(
    inputs: &SeedInputs<'_, S>,
    lookup: BlastAaLookupTable,
    x_pair_par_batch: Option<x_pair_par::XPairPar<S>>,
) -> SeedBatch<S> {
    let &SeedInputs {
        x_board_batch,
        subjects,
        diag_array_size,
        contexts,
        sc,
        extension_debug_enabled,
    } = inputs;
    // A batch without the unit path puts no units on the board: its subjects' slots are closed
    // now, so the subject loop's threads without a subject return at once (`XBoard::help`).
    if let Some((_, batch)) = x_board_batch
        .as_ref()
        .filter(|_| x_pair_par_batch.is_none())
    {
        for s in 0..subjects.count() {
            batch.close_without_units(s);
        }
    }
    if crate::algorithm::tblastx::x_work_board::on() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXBOARD", x_board_batch.is_some());
    }
    // Whether the batch puts units on the board (`x_board_help` in `subject_loop::run_subject_loop`).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let x_board_unit_path = x_pair_par_batch.is_some();
    // Dispatch point 1 of the seed release (per batch): the batch's tables move
    // into two holders, `x_seed_batch.lookup` (the lookup table) and `x_seed_batch.units` (the
    // unit path with its group tables). From here on they are reached only through the stages
    // that each call of `process_subject` opens (dispatch point 2), and the
    // batch's last units stage frees the unit path (dispatch point 3) and its last lookup stage
    // the lookup table (dispatch point 5). `x_worker_diag_cells` is the size the workers'
    // diagonal tables are created with: `diag_array_size`, or 0 when the
    // unit path takes every subject, so that the reference frame loop, the
    // only reader of that table, runs for none.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1004
    // ```c
    // if ((status = BlastExtendWordNew(query->length, word_params,
    //                                 &aux_struct->ewp)) != 0)
    //   return status;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
    // ```c
    // word_params = BlastInitialWordParametersFree(word_params);
    // s_BlastCoreAuxStructFree(aux_struct);
    // ```
    // NCBI keeps its lookup table and diagonal table through the subject loop; when LOSAT
    // frees its copies is not a value NCBI computes (x_seed_release.rs).
    let x_seed_batch = x_seed_release::XSeedBatch::new(lookup, x_pair_par_batch, subjects.count());
    let x_worker_diag_cells = x_seed_batch.worker_table_cells(diag_array_size as usize);
    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the reference frame loops' extensions may take the walks
    // the sentinels stop (not with the debug output of `extend_hit_two_hit_fitting`); another matrix's
    // need its byte table (`TbxScore::x_sent_walks`).
    let tbx_sent = sc.x_sent_walks()
        && !extension_debug_enabled
        && x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice()));
    SeedBatch {
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        x_board_unit_path,
        x_seed_batch,
        x_worker_diag_cells,
        tbx_sent,
    }
}
