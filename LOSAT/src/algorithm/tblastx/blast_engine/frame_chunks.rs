//! The chunks of one TBLASTX subject frame: the chunk driver, which splits the frame into NCBI's
//! subject chunks, scans each with the frame loop's `scan_chunk` and merges its HSPs into the
//! frame's list (`search_frame_chunks`, each chunk's result through `add_chunk_result`), and the
//! append of the frame's list to the subject's (`append_frame_hits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-593 (the chunk loop of
//! `s_BlastSearchEngineOneContext`, :423-602) and :844 (`Blast_HSPListAppend` in the frame loop of
//! `s_BlastSearchEngineCore`, :805-855; c++/src/algo/blast/core/blast_hits.c:2809-2855).
//!
//! This module is a child of `run_impl` for the chunk helpers it keeps private
//! (`SubjectSplitState`, `merge_tblastx_subject_chunk_hits`, `advance_tblastx_diag_offset`). Its
//! functions are statements of the frame loop of `process_subject` moved here in their order. The
//! scan closures stay in the frame loop, which lends `scan_chunk`, the worker's offset pairs and
//! diagonal table, and what the chunks' results are added to.

use super::{
    merge_one_hit_subject_chunk_hits, merge_tblastx_subject_chunk_hits,
    sort_ungapped_hits_by_score_ncbi, tblastx_scan_chunk_size_for_run,
    ungapped_hits_is_sorted_by_score_ncbi, DiagStruct, OffsetPair, OneHitChunkMerge, SubjectChunk,
    SubjectChunkStatus, SubjectSplitState, TblastxChunkScanResult, TbxScore, UngappedHit,
    DBSEQ_CHUNK_OVERLAP,
};
// The parallel chunks (builds with threads only).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    super::advance_tblastx_diag_offset,
    rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator},
};

/// How the chunks of a subject frame are searched, made per frame by `process_subject` under the
/// names it reads them by: LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS (`use_serial_scan_chunks`, with the
/// search's `num_threads` for the scan chunks' size), the parallel chunks (`use_parallel_chunks`),
/// the longest chunk (`max_dbseq_len`), and for the parallel chunks' own tables the window, the
/// diagonal table's size and the number of offset pairs (builds with threads only).
pub(super) struct ChunkDrive {
    pub(super) use_serial_scan_chunks: bool,
    pub(super) use_parallel_chunks: bool,
    pub(super) num_threads: usize,
    pub(super) max_dbseq_len: usize,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) window: i32,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) diag_array_size: i32,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) offset_pairs_len: usize,
}

/// What each chunk's result is added to: the frame's HSP list and the subject's HSP-saving
/// statistics (`process_subject`'s locals of the same names), and for the one-hit word finder
/// (`-window_size 0`, `one_hit`) where its chunk merge reports a refusal.
pub(super) struct ChunkSums<'s> {
    pub(super) frame_ungapped_hits: &'s mut Vec<UngappedHit>,
    pub(super) stats_hsp_saved: &'s mut usize,
    pub(super) stats_hsp_filtered_by_cutoff: &'s mut usize,
    pub(super) stats_score_distribution: &'s mut Vec<i32>,
    pub(super) one_hit: Option<OneHitChunkMerge<'s>>,
}

/// Searches the chunks of one subject frame (`split_state`) with `scan_chunk`, which scans a chunk
/// with the worker's `offset_pairs`, `diag_array` and `diag_offset` (the parallel chunks with
/// tables of their own, after which `diag_offset` walks the chunks in order), and adds each
/// chunk's result to `sums` in chunk order: serially in scan chunks
/// (LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS), on the pool (the parallel chunks), or one after the other.
/// `scan_chunk` is `Sync` for the parallel chunks.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-593 (the chunk loop of
/// `s_BlastSearchEngineOneContext`: `s_GetNextSubjectChunk` :221, the word finder :494-498,
/// `Blast_HSPListsMerge` :581-584).
pub(super) fn search_frame_chunks<S: TbxScore, F>(
    drive: &ChunkDrive,
    split_state: &mut SubjectSplitState,
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    mut sums: ChunkSums<'_>,
    scan_chunk: &F,
    sc: S,
) where
    F: Fn(
            SubjectChunk,
            &mut [OffsetPair],
            &mut [DiagStruct],
            &mut i32,
            Option<usize>,
            S,
        ) -> TblastxChunkScanResult
        + Sync,
{
    let &ChunkDrive {
        use_serial_scan_chunks,
        use_parallel_chunks,
        num_threads,
        max_dbseq_len,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        window,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        diag_array_size,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        offset_pairs_len,
    } = drive;
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
    // while (scan_range[1] <= scan_range[2]) { hits = scansub(...); }
    if use_serial_scan_chunks {
        loop {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
            // ```c
            // while (TRUE) {
            //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
            //                                   dbseq_chunk_overlap);
            //    if (status == SUBJECT_SPLIT_DONE) break;
            //    BlastInitHitListReset(init_hitlist);
            //    aux_struct->WordFinder(..., init_hitlist, ...);
            //    BLAST_GetUngappedHSPList(..., &hsp_list);
            //    Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
            //    status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
            // }
            // ```
            let chunk = match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                SubjectChunkStatus::Done => break,
                SubjectChunkStatus::Ok(chunk) => chunk,
            };

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-123
            // ```c
            // for (s = s_first; s <= s_last; s++) {
            //     ...
            //     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
            // }
            // ```
            //
            // Experimental scan chunks split only the scan walk inside this
            // NCBI search unit. The reducer still extends against `subject`,
            // the full real chunk, and `Blast_ExtendWordExit` runs once, at the end
            // of the chunk's scan.
            let scan_size = tblastx_scan_chunk_size_for_run(chunk.length, num_threads);
            // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
            // while (scan_range[1] <= scan_range[2]) { hits = scansub(...); }
            crate::utils::threading::report_stage(
                "tblastx",
                "serial_scan_chunks",
                chunk.length.div_ceil(scan_size),
                false,
            );
            let result = scan_chunk(
                chunk,
                offset_pairs,
                diag_array,
                diag_offset,
                Some(scan_size),
                sc,
            );
            add_chunk_result(&mut sums, result);
        }
    } else if use_parallel_chunks {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
            // ```c
            // while (TRUE) {
            //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
            //                                   dbseq_chunk_overlap);
            //    if (status == SUBJECT_SPLIT_DONE) break;
            //    ...
            //    status = Blast_HSPListsMerge(...);
            // }
            // ```
            let mut chunks = Vec::new();
            loop {
                match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                    SubjectChunkStatus::Done => break,
                    SubjectChunkStatus::Ok(chunk) => chunks.push(chunk),
                }
            }
            // NCBI reference: c++/src/algo/blast/core/blast_engine.c:452-584
            // status = s_GetNextSubjectChunk(subject, &backup, ...);
            crate::utils::threading::report_stage(
                "tblastx",
                "subject_chunks",
                chunks.len(),
                chunks.len() > 1,
            );
            if chunks.len() <= 1 {
                for chunk in chunks {
                    let result = scan_chunk(chunk, offset_pairs, diag_array, diag_offset, None, sc);
                    add_chunk_result(&mut sums, result);
                }
            } else {
                let mut chunk_results: Vec<TblastxChunkScanResult> =
                    Vec::with_capacity(chunks.len());
                chunks
                    .par_iter()
                    .map_init(
                        || {
                            (
                                vec![OffsetPair::default(); offset_pairs_len],
                                vec![DiagStruct::default(); diag_array_size as usize],
                                sc,
                            )
                        },
                        |state, chunk| {
                            let (offset_pairs, diag_array, sc) = state;
                            for diag in diag_array.iter_mut() {
                                *diag = DiagStruct::default();
                            }
                            let mut chunk_diag_offset = window;
                            scan_chunk(
                                *chunk,
                                offset_pairs,
                                diag_array,
                                &mut chunk_diag_offset,
                                None,
                                *sc,
                            )
                        },
                    )
                    .collect_into_vec(&mut chunk_results);
                chunk_results.sort_by_key(|result| result.chunk.offset);
                for result in chunk_results {
                    // Keep the canonical per-subject diagonal offset moving in the
                    // same chunk order as NCBI before the next subject frame starts.
                    // The experimental workers above use local diagonal tables; the
                    // ordered merge below remains the NCBI `Blast_HSPListsMerge` path.
                    // References: blast_engine.c:561-584, blast_extend.c:167-173.
                    advance_tblastx_diag_offset(
                        diag_offset,
                        diag_array,
                        window,
                        result.chunk.length,
                    );
                    add_chunk_result(&mut sums, result);
                }
            }
        }
        #[cfg(any(
            not(feature = "parallel"),
            all(target_arch = "wasm32", not(feature = "wasm-threads"))
        ))]
        {
            unreachable!("parallel chunk mode is disabled for this target");
        }
    } else {
        loop {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
            // ```c
            // while (TRUE) {
            //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
            //                                   dbseq_chunk_overlap);
            //    if (status == SUBJECT_SPLIT_DONE) break;
            //    BlastInitHitListReset(init_hitlist);
            //    aux_struct->WordFinder(..., init_hitlist, ...);
            //    BLAST_GetUngappedHSPList(..., &hsp_list);
            //    Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
            //    status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
            // }
            // ```
            let chunk = match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                SubjectChunkStatus::Done => break,
                SubjectChunkStatus::Ok(chunk) => chunk,
            };

            // NCBI: BlastInitHitListReset(init_hitlist) - reset per chunk.
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
            let result = scan_chunk(chunk, offset_pairs, diag_array, diag_offset, None, sc);
            add_chunk_result(&mut sums, result);
        }
    }
}

/// Adds the result of one chunk to `sums`: its HSP-saving statistics, and its HSPs merged into
/// the frame's list where they meet the list's HSPs in the chunks' overlap.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:581-584
/// ```c
///         status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
///                      kHspNumMax, &(backup.offset), INT4_MIN,
///                      overlap, score_options->gapped_calculation,
///                      Blast_ProgramIsMapping(program_number));
/// ```
fn add_chunk_result(sums: &mut ChunkSums<'_>, result: TblastxChunkScanResult) {
    *sums.stats_hsp_saved += result.stats.hsp_saved;
    *sums.stats_hsp_filtered_by_cutoff += result.stats.hsp_filtered_by_cutoff;
    sums.stats_score_distribution
        .extend(result.stats.score_distribution);
    if !result.hits.is_empty() {
        // The lists of the one-hit word finder hold NCBI's `Int4` ends
        // (`merge_one_hit_subject_chunk_hits`).
        if let Some(one_hit) = sums.one_hit {
            merge_one_hit_subject_chunk_hits(
                sums.frame_ungapped_hits,
                result.hits,
                result.chunk.offset,
                result.chunk.overlap,
                one_hit,
            );
        } else {
            merge_tblastx_subject_chunk_hits(
                sums.frame_ungapped_hits,
                result.hits,
                result.chunk.offset,
                result.chunk.overlap,
            );
        }
    }
}

/// Appends a frame's HSP list to the subject's (`combined_ungapped_hits`), sorted by score when
/// the append leaves it out of order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844 (`Blast_HSPListAppend`,
/// c++/src/algo/blast/core/blast_hits.c:2809-2855, which combines the lists with
/// `s_BlastHSPListsCombineByScore`, :2749).
pub(super) fn append_frame_hits(
    combined_ungapped_hits: &mut Vec<UngappedHit>,
    frame_ungapped_hits: Vec<UngappedHit>,
) {
    if !frame_ungapped_hits.is_empty() {
        // NCBI: Blast_HSPListAppend merges per-frame subject HSP lists
        // into the combined translated-subject list, then sorts the
        // combined list by score through s_BlastHSPListsCombineByScore.
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:804-845
        // ```c
        // for (context=first_context; context<=last_context; context++) {
        //     subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
        //     ...
        //     Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2758-2766
        // ```c
        // for (index=combined_hsp_list->hspcnt, index1=0;
        //      index1<hsp_list->hspcnt; index1++) {
        //    combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
        // }
        // combined_hsp_list->hspcnt = new_hspcnt;
        // Blast_HSPListSortByScore(combined_hsp_list);
        // ```
        combined_ungapped_hits.extend(frame_ungapped_hits);
        if !ungapped_hits_is_sorted_by_score_ncbi(combined_ungapped_hits) {
            sort_ungapped_hits_by_score_ncbi(combined_ungapped_hits);
        }
    }
}
