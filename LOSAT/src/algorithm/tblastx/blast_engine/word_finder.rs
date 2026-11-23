//! The word finder of TBLASTX's frame loop: `scan_chunk` scans one chunk of a subject frame with
//! `scan_chunk_one_hit` (`-window_size 0`), `scan_chunk_plain` (the port of NCBI's two-hit loop) or
//! `scan_chunk_bucketed` (SEEDBUCKET's order of the same hits). Each starts the chunk with
//! `begin_chunk_scan` and ends it with `finish_chunk_scan`; the two two-hit scans scan with
//! `scan_next_hits` and run the two-hit test and extension of each hit with `two_hit_one_hit`. They
//! read the subject's values (`WordFinderCx`: the batch's, with the subject's own `cutoff_scores`
//! where NCBI computes the statistics of each subject), the subject frame's (`FrameCx`) and the
//! chunk's (`ChunkCx`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:440-619
//! (`s_BlastAaWordFinder_TwoHit`: the scan :502-507, the hit loop :509-608, `Blast_ExtendWordExit`
//! :614), :712-812 (`s_BlastAaWordFinder_OneHit`) and :200-236 (`BlastAaWordFinder`, which chooses
//! between them and sorts the init hit list, :234).
//!
//! This module is a child of `run_impl` for the chunk types and helpers it keeps private
//! (`SubjectChunk`, `TblastxChunkScanStats`, `TblastxChunkScanResult`,
//! `advance_tblastx_diag_offset`, `adjust_tblastx_chunk_subject_offsets`, the `-lcase_masking`
//! ranges). The frame loop (`frame_chunks::search_subject_frames`) calls `scan_chunk` for each
//! chunk. The helpers are the statements the scans have in common, in their order; each scan keeps
//! what only it does (plain: the `-lcase_masking` ranges, which the one-hit scan also reads, the
//! order check and the scan dump; bucketed: the buckets and their flushes; one-hit: its loop) and
//! the order in which it calls them.

use super::{
    adjust_tblastx_chunk_subject_offsets, advance_tblastx_diag_offset, atomic_max_i32,
    atomic_max_usize, atomic_min_i32, clip_tblastx_seq_ranges_for_scan_interior, convert_coords,
    extend_hit_two_hit_with, get_ungapped_hsp_list, get_ungapped_hsp_list_int4,
    sort_init_hsps_by_score_ncbi, tblastx_scan_interiors, trace_hsp_target,
    trace_init_hsp_if_match, trace_match_target, x_scan_ranges, x_soft_chunk, x_soft_scan_letters,
    x_tbx_sent_frame, x_tbx_sent_subject, AtomicOrdering, AtomicU64, DiagStruct,
    DiagnosticCounters, HashSet, InitHSP, Instant, OffsetPair, QueryBatch, QueryContext,
    QueryFrame, SubjectChunk, TblastxChunkScanResult, TblastxChunkScanStats, TbxScore,
    TranslatedChunk, UngappedHit,
};
use crate::algorithm::common::translated_chunks::ChunkLetters;
use crate::algorithm::tblastx::blast_aascan::{
    determine_scanning_offsets, x_pairs_at_chunk, ScanOrder,
};
use crate::algorithm::tblastx::extension::OneHitAlone;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::algorithm::tblastx::subjects::XReadFailure;

/// The values the word finder reads, made once per batch by `search_query_batch` under the names
/// the frame loop read them by: the diagonal table's window and the word size, the table's mask
/// and size, the diagnostics (`diag_enabled`, `diagnostics`), LOSAT_TIMING (`timing_enabled` and
/// the clocks of the scan and the extension), LOSAT_DEBUG_EXTENSION (`extension_debug_enabled`),
/// the query contexts with their cutoffs and X-drops, the scan's LOSAT_X_TBXSCANBF (`x_scan_bf`)
/// and offset-array size, LOSAT_X_TBXUNGAPSENT (`tbx_sent`), LOSAT_DEBUG_SCAN_SOFF's range
/// (`scan_debug_range`), the search's failure (`x_read_failure`, where the one-hit word finder
/// reports an error), and the query batch as the extensions read it past a context's ends
/// (`query_batch`, which refuses into `x_read_failure`). `Sync`: every worker's `process_subject`
/// reads it (`SubjectStatistics::word_finder_cx`), and where NCBI computes the statistics of each subject copies it once per
/// subject with that subject's own `cutoff_scores` (`context_cutoffs::subject_cutoffs`); the frame
/// loop's two-hit and one-hit word finders save their extensions against the `cutoff_scores` of
/// the one they are given (`FrameCx::cx`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:473-493 (what
/// `s_BlastAaWordFinder_TwoHit` reads of the diagonal table and the lookup table) and :575 (the
/// context's cutoffs).
pub(super) struct WordFinderCx<'a> {
    pub(super) window: i32,
    pub(super) wordsize: i32,
    pub(super) diag_mask: i32,
    pub(super) diag_array_size: i32,
    pub(super) diag_enabled: bool,
    pub(super) diagnostics: &'a DiagnosticCounters,
    pub(super) timing_enabled: bool,
    pub(super) extension_debug_enabled: bool,
    pub(super) contexts_ref: &'a [QueryContext],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) x_dropoff_per_context: &'a [i32],
    pub(super) scan_ns: &'a AtomicU64,
    pub(super) scan_calls: &'a AtomicU64,
    pub(super) ungapped_ns: &'a AtomicU64,
    pub(super) ungapped_calls: &'a AtomicU64,
    pub(super) x_scan_bf: bool,
    pub(super) offset_array_size: i32,
    pub(super) tbx_sent: bool,
    pub(super) scan_debug_range: Option<(i32, i32)>,
    pub(super) x_read_failure: &'a XReadFailure,
    /// The query batch as the extensions read it past a context's ends (its refusals go to
    /// `x_read_failure`; extension/query_batch.rs).
    pub(super) query_batch: QueryBatch<'a>,
    /// `-subject_besthit`: each subject frame's list is filtered after each merge of a chunk
    /// (`frame_chunks::add_chunk_result`).
    pub(super) subject_besthit: bool,
}

/// The subject's and the subject frame's values the word finder reads, made per frame by the
/// frame loop (`frame_chunks::search_subject_frames`) under the names it reads them by: the
/// batch's lookup table as the subject's seed stage lends it (`lookup_ref`, borrowed from
/// `x_seed_stage`, so the frame loop's body holds this struct), the subject's index, length and
/// listed frames, whether it collects the HSP-saving statistics (`is_long_sequence`,
/// `collect_hsp_saving_stats`), and the frame (`s_f_idx`, `s_frame`, its letters with and without
/// the sentinels, and its `-lcase_masking` chunks with their ranges, `None` without lower-case
/// letters). It also carries the subject's `WordFinderCx` (`cx`: the batch's, with the subject's
/// `cutoff_scores`), which the four chunk scans read through it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812 (the subject frame of
/// the frame loop of `s_BlastSearchEngineCore`).
pub(super) struct FrameCx<'a> {
    // The chunk scans take only the frame and read the subject's values (`WordFinderCx`: the
    // batch's, or their copy with the subject's own `cutoff_scores`) through this field, a
    // pointer loaded from the frame, as the closures they were made from read `word_finder_cx`
    // through their captures. Given as a parameter of its own (`&WordFinderCx`, which LLVM knows
    // nothing writes during the call), its fields were held in registers across the scan loops,
    // and the bucketed order spilled a register in its push loop: +0.85 % instructions on its
    // frame-loop run, -0.32 % this way (callgrind n1, SW7 tx-6).
    // The tx-5 helpers (`begin_chunk_scan`, `scan_next_hits`, `two_hit_one_hit`,
    // `finish_chunk_scan`) still take `cx` as a parameter; every caller passes `frame.cx`.
    pub(super) cx: &'a WordFinderCx<'a>,
    pub(super) lookup_ref: &'a BlastAaLookupTable,
    pub(super) s_idx: usize,
    pub(super) s_len: usize,
    pub(super) s_frames_listed: &'a [QueryFrame],
    pub(super) is_long_sequence: bool,
    pub(super) collect_hsp_saving_stats: bool,
    pub(super) s_f_idx: usize,
    pub(super) s_frame: &'a QueryFrame,
    pub(super) subject_full: &'a [u8],
    pub(super) subject_all: &'a [u8],
    pub(super) x_soft_calls: Option<&'a [(TranslatedChunk, Vec<(i32, i32)>)]>,
}

/// The chunk's values a hit reads, made once per chunk by each two-hit scan before its scan: the
/// chunk, its letters, whether its extensions take the sentinel walks (LOSAT_X_TBXUNGAPSENT), the
/// diagonal offset (fixed during the chunk's scan: only `begin_chunk_scan`, for a chunk shorter
/// than a word, and `finish_chunk_scan` advance it), and the scores of the search's matrix.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:473 (the diagonal
/// offset is read once per call of the word finder).
pub(super) struct ChunkCx<'s, S> {
    pub(super) chunk: SubjectChunk,
    pub(super) subject: &'s [u8],
    pub(super) subject_sent: bool,
    pub(super) diag_offset: i32,
    pub(super) sc: S,
    /// The frame's letters around the chunk that the scan reads (`x_soft_scan_letters`: from
    /// before the chunk for ranges stored out of order, past its end for a first range raised to
    /// one word), which the extension of a hit before the chunk reads
    /// (`extend_hit_two_hit_outside`); the chunk's letters alone without soft masks.
    pub(super) x_around: ChunkLetters<'s>,
}

/// The start of a chunk's scan, which the scans read under the same names: the chunk's
/// empty lists of init HSPs and statistics, its letters, and whether its extensions take the
/// sentinel walks.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:243-256 (the chunk's sequence
/// and length in `s_GetNextSubjectChunk`) and :491 (`BlastInitHitListReset`).
pub(super) struct ChunkStart<'s> {
    pub(super) init_hsps: Vec<InitHSP>,
    pub(super) stats: TblastxChunkScanStats,
    pub(super) subject: &'s [u8],
    pub(super) subject_sent: bool,
}

/// Starts the scan of `chunk` of the frame: its letters (`subject`) and whether its extensions
/// take the sentinel walks (where the batch's may, `tbx_sent`: `WordFinderCx::tbx_sent` for the
/// two-hit scans, `false` for the one-hit one, whose extension has none). A chunk shorter than a
/// word has no hit: its diagonal offset is advanced as after a scan, and its empty result is the
/// `Err`, which the scan returns.
///
/// Callers pass `frame.cx` as `cx`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:243-256 (the chunk of
/// `s_GetNextSubjectChunk`) and c++/src/algo/blast/core/aa_ungapped.c:614 (`Blast_ExtendWordExit`,
/// c++/src/algo/blast/core/blast_extend.c:162-185); for the one-hit word finder
/// (`s_BlastAaWordFinder_OneHit`), :762-766 (the scan range, whose `while` is not entered for a
/// subject shorter than a word) and :807 (its `Blast_ExtendWordExit`).
pub(super) fn begin_chunk_scan<'s>(
    cx: &WordFinderCx<'_>,
    frame: &FrameCx<'s>,
    chunk: SubjectChunk,
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    tbx_sent: bool,
) -> Result<ChunkStart<'s>, TblastxChunkScanResult> {
    let &WordFinderCx {
        window, wordsize, ..
    } = cx;
    let &FrameCx {
        subject_full,
        subject_all,
        ..
    } = frame;
    let init_hsps: Vec<InitHSP> = Vec::new();
    let stats = TblastxChunkScanStats::default();
    let chunk_end = chunk.offset.saturating_add(chunk.length);
    let subject = &subject_all[chunk.offset..chunk_end];
    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the chunk's extensions take the walks the
    // sentinels stop where it is the whole frame; never on wasm32.
    let subject_sent = cfg!(not(target_arch = "wasm32"))
        && x_tbx_sent_subject(
            x_tbx_sent_frame(tbx_sent, subject_full),
            subject_full,
            chunk.offset,
            chunk.length,
        );

    if subject.len() < wordsize as usize {
        // NCBI still advances the diagonal table offset even when no hits can be found.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:499-500 (two-hit: the clamp
        // of the empty scan range; one-hit: the loop is not entered, :762-766) and :614 (two-hit)
        // or :807 (one-hit) (the advance),
        // c++/src/algo/blast/core/blast_extend.c:167-173
        advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);
        return Err(TblastxChunkScanResult {
            chunk,
            hits: Vec::new(),
            stats,
        });
    }
    Ok(ChunkStart {
        init_hsps,
        stats,
        subject,
        subject_sent,
    })
}

/// One scan call of a chunk's scan: the hits of `letters` from `scan_range` into `offset_pairs`
/// (LOSAT_X_TBXSCANBF: the branch-free scan), with LOSAT_TIMING's clock and the diagnostics' count
/// and duplicate check. `None` when the call found no hit and did not move from `prev_scan_left`,
/// where it started: the caller's scan ends there.
///
/// Callers pass `frame.cx` as `cx`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-507 (`scansub` in the scan
/// loop of `s_BlastAaWordFinder_TwoHit`).
// `#[inline(always)]` (SW7 DECISIONS 5): out of line, the scan it calls left the bucketed order
// (`scan_chunk_bucketed`, then a closure), which inlines it: +0.36 % on its frame-loop run
// (callgrind n1, SW7 tx-5).
#[inline(always)]
pub(super) fn scan_next_hits(
    cx: &WordFinderCx<'_>,
    frame: &FrameCx<'_>,
    letters: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    scan_range: &mut [i32; 3],
    prev_scan_left: i32,
) -> Option<i32> {
    let &WordFinderCx {
        diag_enabled,
        diagnostics,
        timing_enabled,
        scan_ns,
        scan_calls,
        x_scan_bf,
        offset_array_size,
        ..
    } = cx;
    let &FrameCx {
        lookup_ref,
        is_long_sequence,
        s_f_idx,
        ..
    } = frame;
    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    let t0 = if timing_enabled {
        Some(Instant::now())
    } else {
        None
    };
    // EXPERIMENT (LOSAT_X_TBXSCANBF): the branch-free scan with the switch.
    let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
        x_scan_bf,
        lookup_ref,
        letters,
        seq_ranges,
        offset_pairs,
        offset_array_size,
        scan_range,
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
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
            let mut seen: HashSet<(u32, u32)> = HashSet::with_capacity(hits as usize);
            let mut duplicate_count = 0usize;
            for i in 0..hits as usize {
                // SAFETY: `i < hits`; a scan call returns at most
                // `offset_array_size` pairs, and `offset_pairs` was made with
                // `offset_pairs_len` (>= `offset_array_size`) entries.
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
        return None;
    }
    Some(hits)
}

/// The two-hit test and extension of one hit (`query_offset`, `subject_offset`) on its diagonal's
/// cell of `diag_array`, as NCBI's loop body: the cell decides whether the hit is skipped, starts a
/// new hit or is extended with the last one; an extension updates the cell, and its HSP is
/// returned when it meets the context's cutoff (the caller saves it). The diagnostics, the
/// HSP-saving statistics (`stats`), LOSAT_TIMING's clock and the HSP trace count on the way.
/// `scan_chunk_plain` runs it on the hits in scan order and `scan_chunk_bucketed` bucket by bucket:
/// it reads and writes only the hit's cell, the sequences and the matrix, so both orders give each
/// diagonal the same sequence of states.
///
/// Callers pass `frame.cx` as `cx`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607 (the body of the hit
/// loop of `s_BlastAaWordFinder_TwoHit`, :440-619; `s_BlastAaExtendTwoHit`, :1089)
// `#[inline(always)]` (SW7 DECISIONS 5): without it the three call sites of each `S` (plain's hit
// loop, the two flushes) call it out of line, and the frame-loop runs rose by up to +50.2 %
// (callgrind n1, SW7 tx-5: 46 instructions a hit instead of 21); with it none rises.
#[inline(always)]
pub(super) fn two_hit_one_hit<S: TbxScore>(
    cx: &WordFinderCx<'_>,
    frame: &FrameCx<'_>,
    chunk_cx: &ChunkCx<'_, S>,
    diag_array: &mut [DiagStruct],
    stats: &mut TblastxChunkScanStats,
    query_offset: u32,
    subject_offset: u32,
) -> Option<InitHSP> {
    let &WordFinderCx {
        window,
        wordsize,
        diag_mask,
        diag_array_size,
        timing_enabled,
        extension_debug_enabled,
        contexts_ref,
        cutoff_scores,
        x_dropoff_per_context,
        ungapped_ns,
        ungapped_calls,
        ..
    } = cx;
    let &FrameCx {
        lookup_ref,
        s_idx,
        s_len,
        collect_hsp_saving_stats,
        s_f_idx,
        s_frame,
        subject_full,
        ..
    } = frame;
    let &ChunkCx {
        chunk,
        subject,
        subject_sent,
        diag_offset,
        sc,
        ..
    } = chunk_cx;
    let diag_ptr = diag_array.as_mut_ptr();

    // NCBI uses Uint4 for offsets; apply unsigned wrapping.
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
    //             c++/src/algo/blast/core/aa_ungapped.c:516
    let diag_coord = (query_offset.wrapping_sub(subject_offset) & (diag_mask as u32)) as usize;

    // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
    let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-531
    if diag_entry.flag() != 0 {
        // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
        // (A signed test: a hit before the chunk, of ranges stored out of order, can be below 0.)
        let subject_plus_offset = subject_offset.wrapping_add(diag_offset as u32);
        if (subject_plus_offset as i32) < diag_entry.last_hit() {
            if cx.diag_enabled {
                cx.diagnostics
                    .base
                    .seeds_masked
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            return None;
        }
        // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
        // [C] diag_array[diag_coord].flag = 0;
        diag_entry.set_last_hit(subject_plus_offset as i32);
        diag_entry.set_flag(0);
        // Track flag reset (hit after previous extension zone)
        if cx.diag_enabled {
            cx.diagnostics
                .base
                .seeds_flag_reset
                .fetch_add(1, AtomicOrdering::Relaxed);
        }
    }
    // [C] else
    else {
        // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
        let last_hit = diag_entry.last_hit() - diag_offset;
        // NCBI uses Uint4 for subject_offset; compute with unsigned wrap.
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
        //             c++/src/algo/blast/core/aa_ungapped.c:535-536
        let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-544
        if diff >= window {
            if cx.diag_enabled {
                cx.diagnostics
                    .base
                    .seeds_second_hit_too_far
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            diag_entry.set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
            return None;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:549-551
        if diff < wordsize {
            if cx.diag_enabled {
                cx.diagnostics
                    .base
                    .seeds_second_hit_overlap
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            return None;
        }

        // NCBI passes Uint4 query_offset into BSearchContextInfo (Int4).
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
        //             c++/src/algo/blast/core/aa_ungapped.c:560
        let ctx_idx = lookup_ref.get_context_idx(query_offset as i32);
        // SAFETY: `ctx_idx < contexts_ref.len()`: `get_context_idx`
        // gives at most `max(num_contexts, 1) - 1`, the builders set
        // `num_contexts` to `contexts.len()`, and `contexts` is not empty
        // (a batch without a valid context returns before the subject
        // loop).
        let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
        let q_raw = query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
        // NCBI uses masked sequence for extension; query->sequence is
        // sequence_start + 1, so offsets are 0-based in that buffer.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:311-315, c++/src/algo/blast/core/blast_util.c:112-116.
        let query_full = &ctx.aa_seq;
        let query = &query_full[1..query_full.len() - 1];

        // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
        let q_minus_diff = query_offset.wrapping_sub(diff as u32);
        if q_minus_diff < ctx.frame_base as u32 {
            if cx.diag_enabled {
                cx.diagnostics
                    .base
                    .seeds_ctx_boundary_fail
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            diag_entry.set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
            return None;
        }

        if cx.diag_enabled {
            cx.diagnostics
                .base
                .seeds_second_hit_window
                .fetch_add(1, AtomicOrdering::Relaxed);
            cx.diagnostics
                .base
                .seeds_passed
                .fetch_add(1, AtomicOrdering::Relaxed);
        }

        // [C] cutoffs = word_params->cutoffs + curr_context;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:686-688
        // SAFETY: `ctx_idx < contexts_ref.len()` (the `ctx` comment above),
        // and `cutoff_scores` is `vec![0; contexts_ref.len()]`.
        let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx) };
        // [C] cutoffs->x_dropoff (per-context x_dropoff)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:579
        // SAFETY: `ctx_idx < contexts_ref.len()` (the `ctx` comment above),
        // and `x_dropoff_per_context` maps `contexts`, one entry each.
        let x_dropoff = unsafe { *x_dropoff_per_context.get_unchecked(ctx_idx) };

        //                                   last_hit + wordsize, subject_offset, query_offset, ...)
        // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
        let t0 = if timing_enabled {
            Some(Instant::now())
        } else {
            None
        };
        // NCBI's walks run on the query batch (`query`, `query_offset`): past the context's ends
        // they go on into it (extension/query_batch.rs; aa_ungapped.c:841, 897).
        let in_batch = frame.cx.query_batch.at(ctx_idx, frame.s_idx);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582 and
        // :1089-1158
        // A first hit before the chunk (its `last_hit + wordsize` below 0: a range of masks stored
        // out of order starts before the chunk) reads the frame's letters around the chunk as
        // NCBI's extension does (`extend_hit_two_hit_outside`); a chunk without soft masks has
        // none. A word past the chunk's last letter (a first range raised to one word at its left,
        // aa_ungapped.c:495-500) is never extended: it is the chunk's first word scanned, so it
        // starts its diagonal (the cells of the chunks and subjects before lie a window or more
        // behind it), and every later word of the chunk starts at least a word before the chunk's
        // end, before it on that diagonal.
        let s_left_off = last_hit + wordsize;
        debug_assert!(
            subject_offset as i32 + wordsize <= chunk.length as i32,
            "a word past the chunk's last letter is never a second hit"
        );
        let (hsp_q_u, hsp_qe_u, hsp_s, hsp_se, score, right_extend, s_last_off) = if s_left_off < 0
        {
            crate::algorithm::tblastx::extension::extend_hit_two_hit_outside(
                sc,
                query,
                in_batch,
                &chunk_cx.x_around,
                s_left_off,
                subject_offset as i32,
                q_raw,
                x_dropoff,
            )
        } else if !subject_sent {
            let (q_start, q_end, s_start, s_end, score, right_extend, s_last_off) =
                extend_hit_two_hit_with(
                    sc,
                    query,
                    in_batch,
                    subject,
                    (last_hit + wordsize) as usize,
                    subject_offset as usize,
                    q_raw,
                    x_dropoff,
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582 (call above):
                    // score = s_BlastAaExtendTwoHit(...);
                    // Search-local LOSAT diagnostic flag only.
                    extension_debug_enabled,
                );
            let (s_start, s_end, s_last_off) = (s_start as i32, s_end as i32, s_last_off as i32);
            (
                q_start,
                q_end,
                s_start,
                s_end,
                score,
                right_extend,
                s_last_off,
            )
        } else {
            // SAFETY: the contexts passed `x_tbx_sent_contexts` (tbx_sent)
            // and `x_tbx_sent_subject` accepted the chunk
            // (`subject_sent`).
            let (q_start, q_end, s_start, s_end, score, right_extend, s_last_off) = unsafe {
                sc.sent_two_hit(
                    query_full,
                    in_batch,
                    subject_full,
                    (last_hit + wordsize) as usize,
                    subject_offset as usize,
                    q_raw,
                    x_dropoff,
                )
            };
            let (s_start, s_end, s_last_off) = (s_start as i32, s_end as i32, s_last_off as i32);
            (
                q_start,
                q_end,
                s_start,
                s_end,
                score,
                right_extend,
                s_last_off,
            )
        };
        if let Some(t0) = t0 {
            ungapped_ns.fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
            ungapped_calls.fetch_add(1, AtomicOrdering::Relaxed);
        }

        let hsp_q: i32 = hsp_q_u as i32;
        let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;

        if cx.diag_enabled {
            cx.diagnostics
                .base
                .ungapped_extensions
                .fetch_add(1, AtomicOrdering::Relaxed);
            if right_extend {
                cx.diagnostics
                    .base
                    .ungapped_two_hit_extensions
                    .fetch_add(1, AtomicOrdering::Relaxed);
            } else {
                cx.diagnostics
                    .base
                    .ungapped_one_hit_extensions
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            if hsp_len > 0 {
                cx.diagnostics
                    .base
                    .extension_total_length
                    .fetch_add(hsp_len as usize, AtomicOrdering::Relaxed);
                atomic_max_usize(&cx.diagnostics.base.extension_max_length, hsp_len as usize);
            }
        }

        // NCBI: Update diagonal state based on right_extend
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:596-606
        if right_extend {
            diag_entry.set_flag(1);
            diag_entry.set_last_hit(s_last_off - (wordsize - 1) + diag_offset);
            if cx.diag_enabled {
                cx.diagnostics
                    .base
                    .mask_updates
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
        } else {
            diag_entry.set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-591 (Extension後のcutoffチェック)
        if collect_hsp_saving_stats {
            if score >= cutoff {
                stats.score_distribution.push(score);
                stats.hsp_saved += 1;
            } else {
                stats.hsp_filtered_by_cutoff += 1;
            }
        }
        if score >= cutoff {
            if cx.diag_enabled {
                cx.diagnostics
                    .ungapped_only_hits
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }

            // Extra debug for a traced HSP: print seed/extension inputs and cutoffs.
            if let Some(target) = trace_hsp_target() {
                // Compute outfmt coords for this candidate init-hsp (same logic as trace_init_hsp_if_match).
                // NCBI offsets are 0-based in query/subject->sequence buffers.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768, c++/src/algo/blast/core/blast_aascan.c:110-113.
                let q_aa_start = hsp_q_u;
                let q_aa_end = hsp_qe_u;
                let s_aa_start = (hsp_s + chunk.offset as i32) as usize;
                let s_aa_end = (hsp_se + chunk.offset as i32) as usize;
                let (q_start_dna, q_end_dna) =
                    convert_coords(q_aa_start, q_aa_end, ctx.frame, ctx.orig_len);
                let (s_start_dna, s_end_dna) =
                    convert_coords(s_aa_start, s_aa_end, s_frame.frame, s_len);
                if trace_match_target(target, q_start_dna, q_end_dna, s_start_dna, s_end_dna) {
                    eprintln!(
                            "[TRACE_HSP] seed/extend ctx_idx={} s_f_idx={} q_frame={} s_frame={} score={} cutoff={} x_dropoff={} last_hit={} subject_offset={} chunk_offset={} diff={} q_raw={} query_offset={} diag_coord={} right_extend={} s_last_off={}",
                            ctx_idx,
                            s_f_idx,
                            ctx.frame,
                            s_frame.frame,
                            score,
                            cutoff,
                            x_dropoff,
                            last_hit,
                            subject_offset,
                            chunk.offset,
                            diff,
                            q_raw,
                            query_offset,
                            diag_coord,
                            right_extend,
                            s_last_off,
                        );
                    // NCBI two-hit gating checks (diff/window/wordsize/context).
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:531-606
                    // NCBI uses Uint4 offsets; apply unsigned wrap here as well.
                    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:360-375 BlastSaveInitHsp
            // Store HSP with absolute coordinates (before coordinate conversion)
            //
            // hsp_q is frame-relative coordinate in query->sequence (0-based),
            // frame_base is the context query_offset in the concatenated buffer.
            // NCBI: ungapped_data->q_start is absolute query offset.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768, c++/src/algo/blast/core/blast_query_info.c:311-315.
            let hsp_q_absolute = ctx.frame_base + hsp_q;
            let hsp_qe_absolute = ctx.frame_base + (hsp_q + hsp_len);

            let init = InitHSP {
                q_start_absolute: hsp_q_absolute,
                q_end_absolute: hsp_qe_absolute,
                s_start: hsp_s,
                s_end: hsp_s + hsp_len,
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:589-592
                q_seed_absolute: query_offset as i32,
                s_seed: subject_offset as i32,
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
            return Some(init);
        } else if cx.diag_enabled {
            cx.diagnostics
                .ungapped_cutoff_failed
                .fetch_add(1, AtomicOrdering::Relaxed);
            atomic_min_i32(&cx.diagnostics.ungapped_cutoff_failed_min_score, score);
            atomic_max_i32(&cx.diagnostics.ungapped_cutoff_failed_max_score, score);
        }
    }
    None
}

/// Ends the scan of `chunk`: advances the diagonal offset past the chunk, sorts its init HSPs by
/// score and converts them to the chunk's HSP list (`hsp_list`: `get_ungapped_hsp_list` for the
/// two-hit scans, `get_ungapped_hsp_list_int4` for the one-hit one) in the frame's subject offsets.
///
/// Callers pass `frame.cx` as `cx`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:614 (`Blast_ExtendWordExit`;
/// for the one-hit word finder, :807) and :234 (`Blast_InitHitListSortByScore` in
/// `BlastAaWordFinder`, :200-236), c++/src/algo/blast/core/blast_engine.c:561-562
/// (`BLAST_GetUngappedHSPList`) and :578 (`Blast_HSPListAdjustOffsets`).
pub(super) fn finish_chunk_scan<L>(
    cx: &WordFinderCx<'_>,
    frame: &FrameCx<'_>,
    chunk: SubjectChunk,
    mut init_hsps: Vec<InitHSP>,
    stats: TblastxChunkScanStats,
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    hsp_list: L,
) -> TblastxChunkScanResult
where
    L: FnOnce(Vec<InitHSP>, &[QueryContext], &[QueryFrame]) -> Vec<UngappedHit>,
{
    let &WordFinderCx {
        window,
        contexts_ref,
        ..
    } = cx;
    let &FrameCx {
        s_frames_listed, ..
    } = frame;
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:614
    advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234-235
    // BlastAaWordFinder sorts the chunk's init hit list by
    // score_compare_match (blast_extend.c:274-310) before
    // BLAST_GetUngappedHSPList sees it.
    sort_init_hsps_by_score_ncbi(&mut init_hsps);
    let mut hits = if init_hsps.is_empty() {
        Vec::new()
    } else {
        // NCBI: BLAST_GetUngappedHSPList equivalent - per chunk conversion,
        // then Blast_HSPListAdjustOffsets and Blast_HSPListsMerge.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-584
        hsp_list(init_hsps, contexts_ref, s_frames_listed)
    };
    adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
    TblastxChunkScanResult { chunk, hits, stats }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// This function is the reference loop: the port of the C loop above. The dispatcher
// `scan_chunk` below chooses between it, `scan_chunk_bucketed` and `scan_chunk_one_hit`.
// The chunk's start and end are shared by the three scans; the scan calls and the per-hit body by
// the two two-hit scans.
/// Scans `chunk` of the frame with NCBI's two-hit word finder, the hits of each scan call in scan
/// order, and returns the chunk's HSP list. The `-lcase_masking` ranges (`FrameCx::x_soft_calls`,
/// which the one-hit scan also reads), the order check and the scan dump (LOSAT_DEBUG_SCAN_SOFF)
/// are its own.
fn scan_chunk_plain<S: TbxScore>(
    frame: &FrameCx<'_>,
    chunk: SubjectChunk,
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    scan_chunk_size: Option<usize>,
    sc: S,
) -> TblastxChunkScanResult {
    // The subject's values through the frame (`FrameCx::cx`: the batch's, or their copy with the
    // subject's own `cutoff_scores`, `WordFinderCx`).
    let cx = frame.cx;
    let &WordFinderCx {
        wordsize,
        scan_debug_range,
        tbx_sent,
        ..
    } = cx;
    let &FrameCx {
        s_f_idx,
        subject_full,
        x_soft_calls,
        ..
    } = frame;
    // The chunk's letters and its empty lists; a chunk shorter than a word ends here
    // (`begin_chunk_scan`).
    let ChunkStart {
        mut init_hsps,
        mut stats,
        subject,
        subject_sent,
    } = match begin_chunk_scan(cx, frame, chunk, diag_array, diag_offset, tbx_sent) {
        Ok(start) => start,
        Err(short) => return short,
    };

    // NCBI subject seq_ranges are used by s_DetermineScanningOffsets (c++/src/algo/blast/core/masksubj.inl).
    // With no subject masking, the range is [0, subject->length].
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-511,
    // c++/src/algo/blast/core/masksubj.inl:43-58
    let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
    // `-lcase_masking`: this chunk's ranges, whole (`x_soft_chunk`).
    let x_soft = x_soft_calls.map(|calls| x_soft_chunk(calls, &chunk));
    if let Some(None) = x_soft {
        return TblastxChunkScanResult {
            chunk,
            hits: Vec::new(),
            stats,
        };
    }
    let x_soft = x_soft.flatten();
    // The scan's letters, from `x_lead` letters before the chunk (`x_soft_scan_letters`).
    let (x_scan_letters, x_lead) =
        x_soft_scan_letters(x_soft, subject, subject_full, &chunk, wordsize as usize);
    let scan_interiors = tblastx_scan_interiors(
        chunk.length,
        if x_soft.is_some() {
            None
        } else {
            scan_chunk_size
        },
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58 (cited at `ScanOrder`)
    // The order check (debug builds, `ScanOrder`): the chunk's hits come in the order of NCBI's
    // scan of its ranges, across the scan calls and the interiors, which cut the chunk's one range
    // without soft masks.
    let order_ranges: &[(i32, i32)] = x_soft.unwrap_or(&base_seq_ranges);
    let mut scan_order = ScanOrder::default();
    // What a hit reads of this chunk (`ChunkCx`), once per chunk: the
    // diagonal offset does not change during the chunk's scan.
    let chunk_cx = ChunkCx {
        chunk,
        subject,
        subject_sent,
        diag_offset: *diag_offset,
        sc,
        x_around: ChunkLetters {
            letters: x_scan_letters,
            lead: x_lead,
            length: chunk.length,
        },
    };
    for (interior_start, interior_end) in scan_interiors.into_iter() {
        // The ranges on the scan's letters (`x_scan_ranges`, moved by `x_lead`).
        let seq_ranges = match x_soft {
            Some(ranges) => x_scan_ranges(ranges, x_lead),
            None => clip_tblastx_seq_ranges_for_scan_interior(
                &base_seq_ranges,
                interior_start,
                interior_end,
                wordsize as usize,
                subject.len(),
            ),
        };
        if seq_ranges.is_empty() {
            continue;
        }
        // [C] scan_range[0] = 0;
        // [C] scan_range[1] = subject->seq_ranges[0].left;
        // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
        let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
        // [C] if (scan_range[2] < scan_range[1]) scan_range[2] = scan_range[1];
        // (aa_ungapped.c:499-500: only a masked subject's first range can be so short.)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1098-1101,
        // c++/src/algo/blast/core/blast_stat.c:1559-1575 and
        // c++/src/algo/blast/core/masksubj.inl:50-58
        // A first range raised to one word on the last letters of a frame's last chunk gives a
        // word that holds the frame's closing NULLB (and NCBI's next byte of its buffer of the
        // six frames), whose scores are BLAST_SCORE_MIN: no lookup word holds it, so NCBI's scan
        // finds no hit there and goes on with the later ranges. That word alone is not scanned
        // here; the scan goes on from the next range as NCBI's does. (On a chunk that is not
        // the frame's last the word lies in the scan's letters, `x_soft_scan_letters`.)
        if x_soft.is_some() && scan_range[2] < scan_range[1] {
            scan_range[2] = scan_range[1];
            if scan_range[1] + wordsize > x_scan_letters.len() as i32 {
                scan_range[1] += 1;
                determine_scanning_offsets(&seq_ranges, wordsize, wordsize, &mut scan_range);
            }
        }

        // [C] while (scan_range[1] <= scan_range[2])
        while scan_range[1] <= scan_range[2] {
            let prev_scan_left = scan_range[1];
            // One scan call (`scan_next_hits`); `None` ends the scan.
            let Some(hits) = scan_next_hits(
                cx,
                frame,
                x_scan_letters,
                &seq_ranges,
                offset_pairs,
                &mut scan_range,
                prev_scan_left,
            ) else {
                break;
            };
            // The hits at the chunk's offsets (NCBI's, below 0 before it, as Uint4).
            if x_lead > 0 {
                x_pairs_at_chunk(&mut offset_pairs[..hits as usize], x_lead);
            }

            // [C] for (i = 0; i < hits; ++i)
            // OPTIMIZATION: Use raw pointers to eliminate bounds checking in hot loop
            let offset_pairs_ptr = offset_pairs.as_ptr();

            for i in 0..hits as usize {
                // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                let pair = unsafe { &*offset_pairs_ptr.add(i) };
                let query_offset = pair.q_off;
                let subject_offset = pair.s_off;
                // No NCBI counterpart: LOSAT's debug check of the hit order (`ScanOrder`).
                debug_assert!(
                    scan_order.takes(order_ranges, wordsize, subject_offset as i32),
                    "TBLASTX scan chunks must replay seeds in the scan's order: the ranges in \
                     list order, the subject offsets rising within each"
                );

                // Debug: dump scan output around a target subject offset.
                // The scan output is produced by s_BlastAaScanSubject.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:83-123
                if let Some((lo, hi)) = scan_debug_range {
                    // NCBI offsets are Uint4; cast for debug range comparison only.
                    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
                    let subject_offset_i32 = subject_offset as i32 + chunk.offset as i32;
                    if subject_offset_i32 >= lo && subject_offset_i32 <= hi {
                        eprintln!(
                            "[DEBUG SCAN_OFF] s_f_idx={} s_off={} local_s_off={} q_off={} scan_range=[{},{}] chunk_offset={} diag_offset={}",
                            s_f_idx,
                            subject_offset_i32,
                            subject_offset,
                            query_offset,
                            prev_scan_left,
                            scan_range[1],
                            chunk.offset,
                            diag_offset
                        );
                    }
                }

                // The two-hit test and extension of the hit
                // (`two_hit_one_hit`); the HSP it returns is saved.
                if let Some(init) = two_hit_one_hit(
                    cx,
                    frame,
                    &chunk_cx,
                    diag_array,
                    &mut stats,
                    query_offset,
                    subject_offset,
                ) {
                    init_hsps.push(init);
                }
            }
        }
    }

    // The chunk's end (`finish_chunk_scan`): the diagonal offset advanced
    // past the chunk, the init HSPs sorted by score and converted to the chunk's HSP
    // list.
    finish_chunk_scan(
        cx,
        frame,
        chunk,
        init_hsps,
        stats,
        diag_array,
        diag_offset,
        get_ungapped_hsp_list,
    )
}

// EXPERIMENT (SEEDBUCKET, graduated on 2026-10-22: no variable): `scan_chunk_plain`
// above is the reference loop.  This variant appends the hits of every scan call to
// per-diagonal-range buckets and runs the same per-hit body
// (`two_hit_one_hit`, which `scan_chunk_plain` runs per hit) bucket by
// bucket, then puts the saved HSPs back into scan order.  See `x_seed_bucket`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// This variant ports the same C loop with one change of order: the hits of a scan call are
// grouped by range of diag_coord and processed range by range (see x_seed_bucket.rs).
// The order of the hits within one diagonal, and so every diagonal cell state, is as in C.
/// Scans `chunk` of the frame as `scan_chunk_plain`, with SEEDBUCKET's order of the hits (bucket
/// by bucket of their diagonals, flushed when the buffer is full and at the chunk's end), and
/// returns the chunk's HSP list in the order the scan order gives.
fn scan_chunk_bucketed<S: TbxScore>(
    frame: &FrameCx<'_>,
    chunk: SubjectChunk,
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    scan_chunk_size: Option<usize>,
    sc: S,
) -> TblastxChunkScanResult {
    // The subject's values through the frame (`FrameCx::cx`: the batch's, or their copy with the
    // subject's own `cutoff_scores`, `WordFinderCx`).
    let cx = frame.cx;
    let &WordFinderCx {
        wordsize,
        diag_mask,
        diag_array_size,
        tbx_sent,
        ..
    } = cx;
    // The chunk's letters and its empty lists; a chunk shorter than a word ends here
    // (`begin_chunk_scan`).
    let ChunkStart {
        mut init_hsps,
        mut stats,
        subject,
        subject_sent,
    } = match begin_chunk_scan(cx, frame, chunk, diag_array, diag_offset, tbx_sent) {
        Ok(start) => start,
        Err(short) => return short,
    };

    // NCBI subject seq_ranges are used by s_DetermineScanningOffsets (c++/src/algo/blast/core/masksubj.inl).
    // With no subject masking, the range is [0, subject->length].
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-511,
    // c++/src/algo/blast/core/masksubj.inl:43-58
    let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
    let scan_interiors = tblastx_scan_interiors(chunk.length, scan_chunk_size);
    let mut buckets = crate::algorithm::tblastx::x_seed_bucket::SeedBuckets::new(
        diag_array_size as u32,
        diag_mask as u32,
    );
    let x_bucket_budget = crate::algorithm::tblastx::x_seed_bucket::budget();
    // Per-flush positions of the HSPs the current flush has saved.
    let mut seq_keys: Vec<u32> = Vec::new();
    // What a hit reads of this chunk (`ChunkCx`), once per chunk: the
    // diagonal offset does not change during the chunk's scan.
    let chunk_cx = ChunkCx {
        chunk,
        subject,
        subject_sent,
        diag_offset: *diag_offset,
        sc,
        // Without soft masks (`scan_chunk` takes this scan only then) every hit and its
        // extension lie in the chunk.
        x_around: ChunkLetters {
            letters: subject,
            lead: 0,
            length: chunk.length,
        },
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516,588-591
    // One flush: the buffered hits go through the per-hit body of `scan_chunk_plain`
    // (`two_hit_one_hit`) bucket by bucket, then the HSPs this flush saved
    // (the tail of `init_hsps`, each with `seq`, the hit's position within the flush)
    // are put back into scan order, the order in which the C loop saves them. The body
    // reads and writes only diag_array[diag_coord], the sequences and the matrix, so
    // applying it bucket by bucket gives each diagonal the same sequence of states. The
    // flushes follow each other in scan order, so `init_hsps` is in scan order after
    // every flush. The macro is expanded at the two flushes below (each its own
    // instance of `SeedBuckets::flush`).
    macro_rules! x_flush {
        () => {{
            let x_flush_start = init_hsps.len();
            buckets.flush(|query_offset, subject_offset, seq| {
                if let Some(init) = two_hit_one_hit(
                    cx,
                    frame,
                    &chunk_cx,
                    diag_array,
                    &mut stats,
                    query_offset,
                    subject_offset,
                ) {
                    init_hsps.push(init);
                    seq_keys.push(seq);
                }
            });
            crate::algorithm::tblastx::x_seed_bucket::restore_flush_order(
                &mut init_hsps,
                x_flush_start,
                &mut seq_keys,
            );
        }};
    }
    for (interior_start, interior_end) in scan_interiors.into_iter() {
        let seq_ranges = clip_tblastx_seq_ranges_for_scan_interior(
            &base_seq_ranges,
            interior_start,
            interior_end,
            wordsize as usize,
            subject.len(),
        );
        if seq_ranges.is_empty() {
            continue;
        }
        // [C] scan_range[0] = 0;
        // [C] scan_range[1] = subject->seq_ranges[0].left;
        // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
        let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];

        // [C] while (scan_range[1] <= scan_range[2])
        while scan_range[1] <= scan_range[2] {
            let prev_scan_left = scan_range[1];
            // One scan call (`scan_next_hits`); `None` ends the scan.
            let Some(hits) = scan_next_hits(
                cx,
                frame,
                subject,
                &seq_ranges,
                offset_pairs,
                &mut scan_range,
                prev_scan_left,
            ) else {
                break;
            };

            // The hits of this scan call are appended to the buckets and
            // processed, bucket by bucket, when the buffer is full.
            //
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
            // The C loop processes each offset pair here. This loop only appends the pairs to buckets,
            // in scan order (a buffered hit keeps its position, which is used to restore that order).
            let offset_pairs_ptr = offset_pairs.as_ptr();
            for i in 0..hits as usize {
                // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                let pair = unsafe { &*offset_pairs_ptr.add(i) };
                buckets.push(pair.q_off, pair.s_off);
            }
            if buckets.flush_due(x_bucket_budget) {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
                // The buffered hits are processed when the budget is reached (a bucket at a time).
                x_flush!();
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
    // The remaining hits are processed at the end of the subject chunk.
    x_flush!();
    // The saved HSPs are in scan order here, as in C: every flush restored its own. The
    // score sort of `finish_chunk_scan` (sort_init_hsps_by_score_ncbi,
    // with its NCBI reference) is stable, so the order of comparator-equal HSPs
    // matters.

    // The chunk's end (`finish_chunk_scan`): the diagonal offset advanced
    // past the chunk, the init HSPs sorted by score and converted to the chunk's HSP
    // list.
    finish_chunk_scan(
        cx,
        frame,
        chunk,
        init_hsps,
        stats,
        diag_array,
        diag_offset,
        get_ungapped_hsp_list,
    )
}

/// The one-hit extension of the word at batch offset `query_offset` (in context `ctx_idx`, `ctx`)
/// and `subject_offset` of the chunk's letters `subject`, on the query batch, from where its first
/// pass on the context alone stopped (`again`: a walk ended at the context's end where NCBI's could
/// go on past the NULLB; `extend_hit_one_hit_alone`, extension/query_batch.rs). Out of the scan's
/// loop, from values the loop keeps.
#[cold]
#[inline(never)]
fn one_hit_in_batch<S: TbxScore>(
    sc: S,
    frame: &FrameCx<'_>,
    ctx: &QueryContext,
    ctx_idx: usize,
    query_offset: u32,
    subject: &[u8],
    subject_offset: u32,
    x_dropoff_per_context: &[i32],
    again: crate::algorithm::tblastx::extension::Again,
) -> Option<crate::algorithm::tblastx::extension::OneHit> {
    let q_raw = query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
    crate::algorithm::tblastx::extension::extend_hit_one_hit_resume(
        sc,
        &ctx.aa_seq[1..ctx.aa_seq.len() - 1],
        frame.cx.query_batch.at(ctx_idx, frame.s_idx),
        subject,
        q_raw,
        subject_offset as usize,
        x_dropoff_per_context[ctx_idx],
        again,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-807
// `s_BlastAaWordFinder_OneHit` (`-window_size 0`: `BlastAaWordFinder` takes it when the
// diagonal table has no `multiple_hits`, aa_ungapped.c:214-224) on one chunk of a frame:
// the scan of `scan_chunk_plain` without the clamp of a short first range (that loop
// then never starts: no `scansub` call, no hit), the diagonal `s - q`, no flag, and
// `last_hit` set after every extension, also one below the cutoff. The extension is
// `s_BlastAaExtendOneHit` with the search's scores (`extension::extend_hit_one_hit`, on
// the two-hit extension's walks), on the context's letters and the chunk's, as the
// two-hit extensions read them. The
// unit path and `scan_chunk_bucketed` are two-hit only: a batch of window 0 is a
// `reference_loop` one for `x_pair_par::decide`, and `scan_chunk` takes this function
// for every chunk.
/// Scans `chunk` of the frame with NCBI's one-hit word finder (`-window_size 0`) and returns the
/// chunk's HSP list. `AROUND`: the chunk's scan reads letters before it (`x_soft_scan_letters`),
/// so its words can lie before it.
fn scan_chunk_one_hit<S: TbxScore, const AROUND: bool>(
    frame: &FrameCx<'_>,
    chunk: SubjectChunk,
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    scan_chunk_size: Option<usize>,
    sc: S,
) -> TblastxChunkScanResult {
    // The subject's values through the frame (`FrameCx::cx`: the batch's, or their copy with the
    // subject's own `cutoff_scores`, `WordFinderCx`).
    let cx = frame.cx;
    let &WordFinderCx {
        wordsize,
        diag_mask,
        contexts_ref,
        cutoff_scores,
        x_dropoff_per_context,
        x_scan_bf,
        offset_array_size,
        x_read_failure,
        ..
    } = cx;
    let &FrameCx {
        lookup_ref,
        s_idx,
        s_len,
        s_f_idx,
        s_frame,
        subject_full,
        x_soft_calls,
        ..
    } = frame;
    // The chunk's letters and its empty lists; a chunk shorter than a word ends here
    // (`begin_chunk_scan`, as in `scan_chunk_plain`: aa_ungapped.c:762-766, `scan_range[2] <
    // scan_range[1]`, no scan; blast_extend.c:167-173). The one-hit extension takes no sentinel
    // walks (`tbx_sent` false).
    let ChunkStart {
        mut init_hsps,
        stats,
        subject,
        ..
    } = match begin_chunk_scan(cx, frame, chunk, diag_array, diag_offset, false) {
        Ok(start) => start,
        Err(short) => return short,
    };
    // The scan ranges of `scan_chunk_plain` (masksubj.inl:43-58; `-lcase_masking`:
    // `x_soft_chunk`, and a chunk without ranges is not searched).
    let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
    let x_soft = x_soft_calls.map(|calls| x_soft_chunk(calls, &chunk));
    if let Some(None) = x_soft {
        return TblastxChunkScanResult {
            chunk,
            hits: Vec::new(),
            stats,
        };
    }
    let x_soft = x_soft.flatten();
    // The scan's letters from `x_lead` letters before the chunk, and the ranges on them
    // (`x_soft_scan_letters`, `x_scan_ranges`), as in `scan_chunk_plain`.
    let (x_scan_letters, x_lead) =
        x_soft_scan_letters(x_soft, subject, subject_full, &chunk, wordsize as usize);
    debug_assert_eq!(AROUND, x_lead > 0, "the caller's choice of the loop");
    let x_around = ChunkLetters {
        letters: x_scan_letters,
        lead: x_lead,
        length: chunk.length,
    };
    let scan_interiors = tblastx_scan_interiors(
        chunk.length,
        if x_soft.is_some() {
            None
        } else {
            scan_chunk_size
        },
    );
    let wordsize_usize = wordsize as usize;
    for (interior_start, interior_end) in scan_interiors.into_iter() {
        let seq_ranges = match x_soft {
            Some(ranges) => x_scan_ranges(ranges, x_lead),
            None => clip_tblastx_seq_ranges_for_scan_interior(
                &base_seq_ranges,
                interior_start,
                interior_end,
                wordsize_usize,
                subject.len(),
            ),
        };
        if seq_ranges.is_empty() {
            continue;
        }
        // [C] scan_range[1] = subject->seq_ranges[0].left;
        // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
        let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
        // [C] while (scan_range[1] <= scan_range[2])
        while scan_range[1] <= scan_range[2] {
            let prev_scan_left = scan_range[1];
            // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
            let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                x_scan_bf,
                lookup_ref,
                x_scan_letters,
                &seq_ranges,
                offset_pairs,
                offset_array_size,
                &mut scan_range,
            );
            if hits == 0 && scan_range[1] == prev_scan_left {
                // The guard of `scan_chunk_plain`.
                break;
            }
            // The hits at the chunk's offsets (NCBI's, below 0 before it, as Uint4).
            if AROUND {
                x_pairs_at_chunk(&mut offset_pairs[..hits as usize], x_lead);
            }
            // [C] for (i = 0; i < hits; ++i)
            for pair in &offset_pairs[..hits as usize] {
                let query_offset = pair.q_off;
                let subject_offset = pair.s_off;
                // [C] diag_coord = (subject_offset - query_offset) & diag_mask;
                // (Uint4 offsets: the difference wraps.)
                let diag_coord =
                    (subject_offset.wrapping_sub(query_offset) & (diag_mask as u32)) as usize;
                let diag_entry = &mut diag_array[diag_coord];
                // [C] diff = subject_offset -
                // [C]     (diag_array[diag_coord].last_hit - diag_offset);
                // (Uint4 minus Int4 in unsigned arithmetic, read as the Int4 `diff`.)
                let diff = subject_offset
                    .wrapping_sub(diag_entry.last_hit().wrapping_sub(*diag_offset) as u32)
                    as i32;
                // [C] if (diff >= 0)
                if diff < 0 {
                    continue;
                }
                // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                let ctx_idx = lookup_ref.get_context_idx(query_offset as i32);
                let ctx = &contexts_ref[ctx_idx];
                let q_raw = query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
                // The context's letters (`query->sequence` from the context's offset,
                // blast_query_info.c:311-315), as the two-hit extensions read them.
                let query = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                // What the scan does with the extension's result, written once for the two passes
                // below (the first pass's result stays in its registers).
                let mut record = |hit: crate::algorithm::tblastx::extension::OneHit| {
                    // [C] if (score >= cutoffs->cutoff_score)
                    // [C]     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                    // [C]                      query_offset, subject_offset, hsp_len, score);
                    if hit.score >= cutoff_scores[ctx_idx] {
                        // NCBI reference (598d8ae6):
                        // c++/src/algo/blast/core/aa_ungapped.c:1054,1083
                        // A word whose running sum is reset after its best part (a low
                        // -threshold allows it) gives a length of 0 or less, saved as any
                        // other: its ends are at or before its starts, and below 0 near a
                        // frame's start. NCBI links it with the others and its
                        // re-evaluation deletes it; LOSAT keeps its `Int4` ends up to the
                        // first linking call and drops it there, or without sum statistics
                        // before the re-evaluation (`chaining.rs`, the impl of
                        // `UngappedHit`; `subject_link::link_subject_hits_int4`,
                        // `subject_link::drop_unread_without_linking`).
                        let hsp_q_absolute = ctx.frame_base + hit.q_start;
                        init_hsps.push(InitHSP {
                            q_start_absolute: hsp_q_absolute,
                            q_end_absolute: hsp_q_absolute + hit.length,
                            s_start: hit.s_start,
                            s_end: hit.s_start + hit.length,
                            q_seed_absolute: query_offset as i32,
                            s_seed: subject_offset as i32,
                            score: hit.score,
                            ctx_idx,
                            s_f_idx,
                            q_idx: ctx.q_idx,
                            s_idx: s_idx as u32,
                            q_frame: ctx.frame,
                            s_frame: s_frame.frame,
                            q_orig_len: ctx.orig_len,
                            s_orig_len: s_len,
                        });
                    }
                    // [C] diag_array[diag_coord].last_hit =
                    // [C]     s_last_off - (wordsize - 1) + diag_offset;
                    diag_entry.set_last_hit(hit.s_last_off - (wordsize - 1) + *diag_offset);
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:785-790 and
                // :1020-1083
                // A word before the chunk (a range of masks stored out of order) is extended on
                // the frame's letters around the chunk, as NCBI's extension reads them
                // (`extend_hit_one_hit_outside`); its HSP can start before the chunk.
                if AROUND && (subject_offset as i32) < 0 {
                    record(
                        crate::algorithm::tblastx::extension::extend_hit_one_hit_outside(
                            sc,
                            query,
                            frame.cx.query_batch.at(ctx_idx, frame.s_idx),
                            &x_around,
                            q_raw,
                            subject_offset as i32,
                            x_dropoff_per_context[ctx_idx],
                        ),
                    );
                    continue;
                }
                // [C] score = s_BlastAaExtendOneHit(matrix, subject, query, ...);
                // On the context alone; where a walk could go on into the query batch past the
                // context's ends, as NCBI's do, again on the batch (`one_hit_in_batch`;
                // extension/query_batch.rs; aa_ungapped.c:841, 897).
                let again = match crate::algorithm::tblastx::extension::extend_hit_one_hit_alone(
                    sc,
                    query,
                    subject,
                    q_raw,
                    subject_offset as usize,
                    x_dropoff_per_context[ctx_idx],
                ) {
                    Some(OneHitAlone::Hit(hit)) => {
                        record(hit);
                        continue;
                    }
                    Some(OneHitAlone::Again(again)) => one_hit_in_batch(
                        sc,
                        frame,
                        ctx,
                        ctx_idx,
                        query_offset,
                        subject,
                        subject_offset,
                        x_dropoff_per_context,
                        again,
                    ),
                    None => None,
                };
                let Some(hit) = again else {
                    // A word past the letters, which the scan does not give.
                    x_read_failure.fail(
                        s_idx,
                        anyhow::anyhow!(
                            "TBLASTX -window_size 0: the word at query offset {} and subject offset {} is past the letters (internal error)",
                            query_offset,
                            subject_offset
                        ),
                    );
                    continue;
                };
                record(hit);
            }
        }
    }
    // [C] Blast_ExtendWordExit(ewp, subject->length);
    // aa_ungapped.c:234 `Blast_InitHitListSortByScore`, then the conversion of
    // `scan_chunk_plain` (blast_engine.c:561-584; `finish_chunk_scan`).
    // The list keeps NCBI's `Int4` ends of an HSP of length 0 or less
    // (`get_ungapped_hsp_list_int4`, `adjust_tblastx_chunk_subject_offsets`).
    finish_chunk_scan(
        cx,
        frame,
        chunk,
        init_hsps,
        stats,
        diag_array,
        diag_offset,
        get_ungapped_hsp_list_int4,
    )
}

// EXPERIMENT (SEEDBUCKET): the scan order (`scan_chunk_plain`) or the bucketed order
// (`scan_chunk_bucketed`, only for tables of at least LOSAT_X_SEEDBUCKET_MIN_CELLS
// cells).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// Dispatch point of SEEDBUCKET: `scan_chunk_plain` (the port of the C loop), or
// `scan_chunk_bucketed` for a table of at least the minimum cells. A traced or debugged
// run always uses `scan_chunk_plain`.
/// Scans `chunk` of the frame with the search's word finder and returns the chunk's HSP list: the
/// one-hit one for `-window_size 0`, else the two-hit one in the scan order or SEEDBUCKET's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-236 (`BlastAaWordFinder`,
/// called as `aux_struct->WordFinder`, c++/src/algo/blast/core/blast_engine.c:494-498)
pub(super) fn scan_chunk<S: TbxScore>(
    frame: &FrameCx<'_>,
    chunk: SubjectChunk,
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    scan_chunk_size: Option<usize>,
    sc: S,
) -> TblastxChunkScanResult {
    // The subject's values through the frame (`FrameCx::cx`: the batch's, or their copy with the
    // subject's own `cutoff_scores`, `WordFinderCx`).
    let cx = frame.cx;
    let &WordFinderCx {
        window,
        scan_debug_range,
        ..
    } = cx;
    let &FrameCx { x_soft_calls, .. } = frame;
    // The one-hit word finder (`-window_size 0`) has a loop of its own, compiled out of
    // line (`x_out_of_line`) so that this function keeps the two-hit scan's code.
    if window == 0 {
        return x_out_of_line(|| {
            // A chunk whose ranges start before it (masks stored out of order,
            // `common::translated_chunks::scan_lead`) is the only one whose words can lie before
            // it: only its loop tests for them (`AROUND`). The choice is made in this out-of-line
            // call, so that the two-hit scan below keeps its code (inlined here, it cost
            // `-db_soft_mask` two-hit searches 2.9 % more instructions).
            let around = x_soft_calls
                .and_then(|calls| x_soft_chunk(calls, &chunk))
                .is_some_and(|ranges| {
                    crate::algorithm::common::translated_chunks::scan_lead(ranges) > 0
                });
            if around {
                scan_chunk_one_hit::<S, true>(
                    frame,
                    chunk,
                    offset_pairs,
                    diag_array,
                    diag_offset,
                    scan_chunk_size,
                    sc,
                )
            } else {
                scan_chunk_one_hit::<S, false>(
                    frame,
                    chunk,
                    offset_pairs,
                    diag_array,
                    diag_offset,
                    scan_chunk_size,
                    sc,
                )
            }
        });
    }
    // The debug/trace paths exist only in the reference loop.
    let x_debugging = scan_debug_range.is_some() || trace_hsp_target().is_some();
    if crate::algorithm::tblastx::x_seed_bucket::bucketed(diag_array.len())
        && !x_debugging
        && x_soft_calls.is_none()
    {
        scan_chunk_bucketed(
            frame,
            chunk,
            offset_pairs,
            diag_array,
            diag_offset,
            scan_chunk_size,
            sc,
        )
    } else {
        scan_chunk_plain(
            frame,
            chunk,
            offset_pairs,
            diag_array,
            diag_offset,
            scan_chunk_size,
            sc,
        )
    }
}

// No NCBI counterpart: a call boundary; it computes no value NCBI computes.
/// Calls `f` in a function of its own that is never inlined, so that the function that calls it
/// (`scan_chunk`, for the one-hit scan of `-window_size 0`) does not take `f`'s code.
#[inline(never)]
fn x_out_of_line<R>(f: impl FnOnce() -> R) -> R {
    f()
}
