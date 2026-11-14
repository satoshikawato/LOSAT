//! The word finder of TBLASTX's frame loop, which its two chunk scans share:
//! `scan_chunk_plain` (the port of NCBI's loop) and `scan_chunk_bucketed` (SEEDBUCKET's order of
//! the same hits) start a chunk with `begin_chunk_scan`, scan it with `scan_next_hits`, run the
//! two-hit test and extension of each hit with `two_hit_one_hit`, and end the chunk with
//! `finish_chunk_scan`. They read the batch's values (`WordFinderCx`), the subject frame's
//! (`FrameCx`) and the chunk's (`ChunkCx`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:440-619
//! (`s_BlastAaWordFinder_TwoHit`: the scan :502-507, the hit loop :509-608, `Blast_ExtendWordExit`
//! :614) and :200-236 (`BlastAaWordFinder`, which sorts the init hit list, :234).
//!
//! This module is a child of `run_impl` for the chunk types and helpers it keeps private
//! (`SubjectChunk`, `TblastxChunkScanStats`, `TblastxChunkScanResult`,
//! `advance_tblastx_diag_offset`, `adjust_tblastx_chunk_subject_offsets`). Its functions are the
//! statements the two closures had in common, moved here in their order; each closure keeps what
//! only it does (plain: the `-lcase_masking` ranges, the order check and the scan dump; bucketed:
//! the buckets and their flushes) and the order in which it calls them.

use super::{
    adjust_tblastx_chunk_subject_offsets, advance_tblastx_diag_offset, atomic_max_i32,
    atomic_max_usize, atomic_min_i32, convert_coords, extend_hit_two_hit_with,
    get_ungapped_hsp_list, sort_init_hsps_by_score_ncbi, trace_hsp_target, trace_init_hsp_if_match,
    trace_match_target, x_tbx_sent_subject, AtomicOrdering, AtomicU64, DiagStruct,
    DiagnosticCounters, HashSet, InitHSP, Instant, OffsetPair, QueryContext, QueryFrame,
    SubjectChunk, TblastxChunkScanResult, TblastxChunkScanStats, TbxScore,
};
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;

/// The batch's values the word finder reads, made once per batch by `search_query_batch` under the
/// names the frame loop read them by: the diagonal table's window and the word size, the table's
/// mask and size, the diagnostics (`diag_enabled`, `diagnostics`), LOSAT_TIMING (`timing_enabled`
/// and the clocks of the scan and the extension), LOSAT_DEBUG_EXTENSION
/// (`extension_debug_enabled`), the query contexts with their cutoffs and X-drops, the scan's
/// LOSAT_X_TBXSCANBF (`x_scan_bf`) and offset-array size, and LOSAT_X_TBXUNGAPSENT (`tbx_sent`).
/// `Sync`: every worker's `process_subject` reads it.
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
}

/// The subject's and the subject frame's values the word finder reads, made per frame by the
/// frame loop of `process_subject` under the names it reads them by: the batch's lookup table as
/// the subject's seed stage lends it (`lookup_ref`, borrowed from `x_seed_stage`, so the frame
/// loop's body holds this struct), the subject's index, length and listed frames, whether it
/// collects the HSP-saving statistics (`is_long_sequence`, `collect_hsp_saving_stats`), and the
/// frame (`s_f_idx`, `s_frame`, its letters with and without the sentinels).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812 (the subject frame of
/// the frame loop of `s_BlastSearchEngineCore`).
pub(super) struct FrameCx<'a> {
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
}

/// The chunk's values a hit reads, made once per chunk by each scan closure before its scan: the
/// chunk, its letters, whether its extensions take the sentinel walks (LOSAT_X_TBXUNGAPSENT), the
/// diagonal offset (fixed during the chunk's scan: only `begin_chunk_scan`, for a chunk shorter
/// than a word, and `finish_chunk_scan` advance it), and the scores of the search's matrix.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:473 (`diag_offset =
/// diag->offset;`, read once per call of the word finder).
pub(super) struct ChunkCx<'s, S> {
    pub(super) chunk: SubjectChunk,
    pub(super) subject: &'s [u8],
    pub(super) subject_sent: bool,
    pub(super) diag_offset: i32,
    pub(super) sc: S,
}

/// The start of a chunk's scan, which the scan closures read under the same names: the chunk's
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
/// take the sentinel walks. A chunk shorter than a word has no hit: its diagonal offset is
/// advanced as after a scan, and its empty result is the `Err`, which the closure returns.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:243-256 (the chunk of
/// `s_GetNextSubjectChunk`) and c++/src/algo/blast/core/aa_ungapped.c:614 (`Blast_ExtendWordExit`,
/// c++/src/algo/blast/core/blast_extend.c:162-185).
pub(super) fn begin_chunk_scan<'s>(
    cx: &WordFinderCx<'_>,
    frame: &FrameCx<'s>,
    chunk: SubjectChunk,
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
) -> Result<ChunkStart<'s>, TblastxChunkScanResult> {
    let &WordFinderCx {
        window,
        wordsize,
        tbx_sent,
        ..
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
        && x_tbx_sent_subject(tbx_sent, subject_full, chunk.offset, chunk.length);

    if subject.len() < wordsize as usize {
        // NCBI still advances the diagonal table offset even when no hits can be found.
        // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:499-500 (the clamp of the empty scan
        // range) and :614 (the advance),
        // ncbi-blast/c++/src/algo/blast/core/blast_extend.c:167-173
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
/// where it started: the closure's scan ends there.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-507 (`scansub` in the scan
/// loop of `s_BlastAaWordFinder_TwoHit`).
// `#[inline(always)]` (SW7 DECISIONS 5): out of line, the scan it calls left the closure of the
// bucketed order, which inlines it: +0.36 % on its frame-loop run (callgrind n1, SW7 tx-5).
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
            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607 (the body of the hit
/// loop of `s_BlastAaWordFinder_TwoHit`, :440-619; `s_BlastAaExtendTwoHit`, :1089)
/// ```c
/// diag_coord = (query_offset - subject_offset) & diag_mask;
/// ...
/// if (diag_array[diag_coord].flag) {
/// ...
///     last_hit = diag_array[diag_coord].last_hit - diag_offset;
///     diff = subject_offset - last_hit;
/// ...
///     if (diff >= window) {
/// ...
///     score = s_BlastAaExtendTwoHit(matrix, subject, query,
/// ...
///     if (score >= cutoffs->cutoff_score)
/// ...
///     if (right_extend) {
///         diag_array[diag_coord].flag = 1;
///         diag_array[diag_coord].last_hit =
///             s_last_off - (wordsize - 1) + diag_offset;
/// ```
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
    } = chunk_cx;
    let diag_ptr = diag_array.as_mut_ptr();

    // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
    // NCBI uses Uint4 for offsets; apply unsigned wrapping.
    // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
    //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:516
    let diag_coord = (query_offset.wrapping_sub(subject_offset) & (diag_mask as u32)) as usize;

    // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
    let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

    // [C] if (diag_array[diag_coord].flag)
    // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:519-531
    if diag_entry.flag() != 0 {
        // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
        let subject_plus_offset = subject_offset.wrapping_add(diag_offset as u32);
        if subject_plus_offset < diag_entry.last_hit() as u32 {
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
        // [C] diff = subject_offset - last_hit;
        // NCBI uses Uint4 for subject_offset; compute with unsigned wrap.
        // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
        //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:535-536
        let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

        // [C] if (diff >= window)
        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:538-544
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

        // [C] if (diff < wordsize)
        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:549-551
        if diff < wordsize {
            if cx.diag_enabled {
                cx.diagnostics
                    .base
                    .seeds_second_hit_overlap
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            return None;
        }

        // [C] curr_context = BSearchContextInfo(query_offset, query_info);
        // NCBI passes Uint4 query_offset into BSearchContextInfo (Int4).
        // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
        //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:560
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
        // Reference: blast_query_info.c:311-315, blast_util.c:112-116.
        let query_full = &ctx.aa_seq;
        let query = &query_full[1..query_full.len() - 1];

        // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:566-573
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
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:686-688
        // ```c
        // Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
        // ```
        // SAFETY: `ctx_idx < contexts_ref.len()` (the `ctx` comment above),
        // and `cutoff_scores` is `vec![0; contexts_ref.len()]`.
        let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx) };
        // [C] cutoffs->x_dropoff (per-context x_dropoff)
        // Reference: aa_ungapped.c:579
        // SAFETY: `ctx_idx < contexts_ref.len()` (the `ctx` comment above),
        // and `x_dropoff_per_context` maps `contexts`, one entry each.
        let x_dropoff = unsafe { *x_dropoff_per_context.get_unchecked(ctx_idx) };

        // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
        //                                   last_hit + wordsize, subject_offset, query_offset, ...)
        // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
        let t0 = if timing_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
            if !subject_sent {
                extend_hit_two_hit_with(
                    sc,
                    query,
                    subject,
                    (last_hit + wordsize) as usize,
                    subject_offset as usize,
                    q_raw,
                    x_dropoff,
                    // NCBI aa_ungapped.c:576-582 (call above):
                    // score = s_BlastAaExtendTwoHit(...);
                    // Search-local LOSAT diagnostic flag only.
                    extension_debug_enabled,
                )
            } else {
                // SAFETY: the contexts passed `x_tbx_sent_contexts` (tbx_sent)
                // and `x_tbx_sent_subject` accepted the chunk
                // (`subject_sent`).
                unsafe {
                    sc.sent_two_hit(
                        query_full,
                        subject_full,
                        (last_hit + wordsize) as usize,
                        subject_offset as usize,
                        q_raw,
                        x_dropoff,
                    )
                }
            };
        if let Some(t0) = t0 {
            ungapped_ns.fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
            ungapped_calls.fetch_add(1, AtomicOrdering::Relaxed);
        }

        let hsp_q: i32 = hsp_q_u as i32;
        let hsp_s: i32 = hsp_s_u as i32;
        let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
        let s_last_off: i32 = s_last_off_u as i32;

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
        // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:596-606
        // if (right_extend) {
        //     diag_array[diag_coord].flag = 1;
        //     diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
        // } else {
        //     diag_array[diag_coord].last_hit = subject_offset + diag_offset;
        // }
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

        // [C] if (score >= cutoffs->cutoff_score)
        // NCBI reference: aa_ungapped.c:575-591 (Extension後のcutoffチェック)
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
                // Reference: blast_gapalign.c:4756-4768, blast_aascan.c:110-113.
                let q_aa_start = hsp_q_u;
                let q_aa_end = hsp_qe_u;
                let s_aa_start = hsp_s_u + chunk.offset;
                let s_aa_end = _hsp_se_u + chunk.offset;
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
            let hsp_q_absolute = ctx.frame_base + hsp_q;
            let hsp_qe_absolute = ctx.frame_base + (hsp_q + hsp_len);

            let init = InitHSP {
                q_start_absolute: hsp_q_absolute,
                q_end_absolute: hsp_qe_absolute,
                s_start: hsp_s,
                s_end: hsp_s + hsp_len,
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:589-592
                // ```c
                // BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                //                  query_offset, subject_offset, hsp_len,
                //                  score);
                // ```
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
/// score and converts them to the chunk's HSP list in the frame's subject offsets.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:614 (`Blast_ExtendWordExit`)
/// and :234 (`Blast_InitHitListSortByScore` in `BlastAaWordFinder`, :200-236),
/// c++/src/algo/blast/core/blast_engine.c:561-562 (`BLAST_GetUngappedHSPList`) and :578
/// (`Blast_HSPListAdjustOffsets`).
pub(super) fn finish_chunk_scan(
    cx: &WordFinderCx<'_>,
    frame: &FrameCx<'_>,
    chunk: SubjectChunk,
    mut init_hsps: Vec<InitHSP>,
    stats: TblastxChunkScanStats,
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
) -> TblastxChunkScanResult {
    let &WordFinderCx {
        window,
        contexts_ref,
        ..
    } = cx;
    let &FrameCx {
        s_frames_listed, ..
    } = frame;
    // [C] Blast_ExtendWordExit(ewp, subject->length);
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:614
    // ```c
    // /* increment the offset in the diagonal array */
    // Blast_ExtendWordExit(ewp, subject->length);
    // ```
    advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:234-235
    // ```c
    // Blast_InitHitListSortByScore(init_hitlist);
    // return status;
    // ```
    // BlastAaWordFinder sorts the chunk's init hit list by
    // score_compare_match (blast_extend.c:274-310) before
    // BLAST_GetUngappedHSPList sees it.
    sort_init_hsps_by_score_ncbi(&mut init_hsps);
    let mut hits = if init_hsps.is_empty() {
        Vec::new()
    } else {
        // NCBI: BLAST_GetUngappedHSPList equivalent - per chunk conversion,
        // then Blast_HSPListAdjustOffsets and Blast_HSPListsMerge.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:561-584
        // ```c
        // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
        //         hit_params->options, &hsp_list);
        // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
        // status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
        //      kHspNumMax, &(backup.offset), INT4_MIN, overlap, ...);
        // ```
        get_ungapped_hsp_list(init_hsps, contexts_ref, s_frames_listed)
    };
    adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
    TblastxChunkScanResult { chunk, hits, stats }
}
