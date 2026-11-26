//! The frames of one TBLASTX subject and their chunks: the frame loop (`search_subject_frames`),
//! the chunk driver, which splits a frame into NCBI's subject chunks, scans each with the word
//! finder (`word_finder::scan_chunk`) and merges its HSPs into the frame's list
//! (`search_frame_chunks`, each chunk's result through `add_chunk_result`), and the append of the
//! frame's list to the subject's (`append_frame_hits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:804-855 (the frame loop of
//! `s_BlastSearchEngineCore`), :478-593 (the chunk loop of `s_BlastSearchEngineOneContext`,
//! :423-602) and :844 (`Blast_HSPListAppend`; c++/src/algo/blast/core/blast_hits.c:2809-2855).
//!
//! This module is a child of `run_impl` for the chunk helpers it keeps private
//! (`SubjectSplitState`, `merge_tblastx_subject_chunk_hits`, `advance_tblastx_diag_offset`). Its
//! functions are the frame loop of `process_subject` moved here in its order; `process_subject`
//! lends the worker's offset pairs, diagonal table and diagonal offset, and what the frames' lists
//! and statistics are added to.

use super::word_finder::{self, WordFinderCx};
use super::{
    merge_tblastx_subject_chunk_hits, sort_ungapped_hits_by_score_ncbi,
    tblastx_max_dbseq_len_for_run, tblastx_scan_chunk_size_for_run,
    ungapped_hits_is_sorted_by_score_ncbi, x_db_mask_calls, x_soft_subject_calls, DiagStruct,
    OffsetPair, OneHitChunkMerge, QueryFrame, SubjectChunk, SubjectChunkStatus, SubjectSplitState,
    Subjects, TblastxChunkScanResult, TbxScore, UngappedHit, DBSEQ_CHUNK_OVERLAP,
    TWO_HIT_CHUNK_MERGE_IS_DEFINED,
};
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
// The parallel chunks (builds with threads only).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    super::advance_tblastx_diag_offset,
    rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator},
};

/// The subject's values the frame loop reads, made per subject by `process_subject` under the names
/// it reads them by: the search's subjects and whether they may have lower-case letters
/// (`-lcase_masking`), how the chunks are searched (LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS, the parallel
/// chunks, the search's threads; the size of the parallel chunks' offset arrays in builds with
/// threads), the batch's lookup table as the subject's seed stage lends it (`lookup_ref`, borrowed
/// from `x_seed_stage`: the struct is a temporary of the call, so the borrow ends before the
/// stage), the subject's index and length, its preliminary and listed frames, and whether it
/// collects the HSP-saving statistics.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:702-855 (what
/// `s_BlastSearchEngineCore` reads of the subject in its frame loop, :804-855).
pub(super) struct SubjectScan<'a> {
    pub(super) subjects: Subjects<'a>,
    pub(super) x_soft_subjects: bool,
    /// A masked `-db` search (`-db_soft_mask`, `-db_hard_mask`; series NP, M4b PR 2): the
    /// subject's ranges in bases (`subject_prep::SubjectFrames::x_db_mask`).
    pub(super) x_db_mask: Option<&'a [(i32, i32)]>,
    /// Whether this `-db_soft_mask` subject has a chunk whose scan reads letters before it
    /// (`scan_lead` above zero, series NP, M4b-T3): all six frames run their chunks in order on
    /// the worker's single diagonal table initialized to NCBI's one-thread start.
    pub(super) x_lead_subject: bool,
    /// Whether this `-db_hard_mask` subject has a start from the database's overlapping-mask walk
    /// (`XDbStream::hard_start`, series NP, M4b-3c part 2): all six frames run their chunks in
    /// order on the worker's single diagonal table, and two-hit chunks shorter than the word go
    /// through `common::hard_chunk_scan::OneWord`.
    pub(super) x_hard_start: bool,
    pub(super) use_serial_scan_chunks: bool,
    pub(super) use_parallel_chunks: bool,
    pub(super) num_threads: usize,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) offset_pairs_len: usize,
    pub(super) lookup_ref: &'a BlastAaLookupTable,
    pub(super) s_idx: usize,
    pub(super) s_len: usize,
    pub(super) s_frames_preliminary: &'a [QueryFrame],
    pub(super) s_frames_listed: &'a [QueryFrame],
    pub(super) is_long_sequence: bool,
    pub(super) collect_hsp_saving_stats: bool,
}

/// What the frame loop adds to: the subject's HSP-saving statistics and its HSP list
/// (`process_subject`'s locals of the same names, which it reads after the loop).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844 (`Blast_HSPListAppend` of
/// each frame's list to the subject's, `hsp_list_out`).
pub(super) struct SubjectSums<'s> {
    pub(super) stats_hsp_saved: &'s mut usize,
    pub(super) stats_hsp_filtered_by_cutoff: &'s mut usize,
    pub(super) stats_score_distribution: &'s mut Vec<i32>,
    pub(super) combined_ungapped_hits: &'s mut Vec<UngappedHit>,
}

/// The frame loop of a subject: for each preliminary frame, its `-lcase_masking` ranges and its
/// `word_finder::FrameCx`, the search of its chunks with the word finder (`search_frame_chunks`
/// with `word_finder::scan_chunk`), and the append of the frame's list to the subject's. The
/// worker's `offset_pairs`, `diag_array` and `diag_offset` carry over from frame to frame, as
/// NCBI's diagonal table does.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:804-855 (the frame loop of
/// `s_BlastSearchEngineCore`)
pub(super) fn search_subject_frames<S: TbxScore>(
    cx: &WordFinderCx<'_>,
    subject: &SubjectScan<'_>,
    offset_pairs: &mut [OffsetPair],
    diag_array: &mut [DiagStruct],
    diag_offset: &mut i32,
    sums: SubjectSums<'_>,
    sc: S,
) {
    let &WordFinderCx {
        window,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        diag_array_size,
        x_read_failure,
        subject_besthit,
        ..
    } = cx;
    let &SubjectScan {
        subjects,
        x_soft_subjects,
        x_db_mask,
        x_lead_subject,
        x_hard_start,
        use_serial_scan_chunks,
        use_parallel_chunks,
        num_threads,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        offset_pairs_len,
        lookup_ref,
        s_idx,
        s_len,
        s_frames_preliminary,
        s_frames_listed,
        is_long_sequence,
        collect_hsp_saving_stats,
    } = subject;
    let SubjectSums {
        stats_hsp_saved,
        stats_hsp_filtered_by_cutoff,
        stats_score_distribution,
        combined_ungapped_hits,
    } = sums;
    for (s_f_idx, s_frame) in s_frames_preliminary.iter().enumerate() {
        // NCBI: Diagonal state is NOT reset between subject frames.
        // NCBI shares ewp (diag_table) across all 6 subject frame iterations.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-855
        // The per-subject reset (`process_subject`, before the unit path) matches NCBI's
        // Blast_ExtendWordNew calloc.
        // Blast_ExtendWordExit at the end of each subject chunk increments diag_offset.

        let subject_full = &s_frame.aa_seq;
        // NCBI subject->sequence points past the leading NULLB sentinel, so offsets
        // are 0-based from the first residue; see blast_engine.c:811-812 and
        // blast_util.c:112-116.
        let subject_all = &subject_full[1..subject_full.len() - 1];
        let s_aa_len = s_frame.aa_len;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:804-855
        let mut frame_ungapped_hits: Vec<UngappedHit> = Vec::new();
        let mut split_state = SubjectSplitState::new(s_aa_len);
        let max_dbseq_len = tblastx_max_dbseq_len_for_run();
        // `-lcase_masking`: the scan ranges of this frame's chunks, `None` for a subject
        // without lower-case letters. A masked `-db` search: with `-db_soft_mask` the scan ranges
        // of the frame's chunks, with `-db_hard_mask` its chunks, each scanned whole
        // (`SubjectSplitState::listed`).
        let x_soft_calls = match x_db_mask {
            Some(ranges) => {
                let hard = matches!(subjects, Subjects::Stream(stream)
                    if stream.mask().is_some_and(|mask| mask.hard()));
                let calls = x_db_mask_calls(ranges, hard, s_frame, s_frames_preliminary);
                if hard {
                    split_state = SubjectSplitState::listed(s_aa_len, &calls);
                    None
                } else {
                    Some(calls)
                }
            }
            None if x_soft_subjects => {
                x_soft_subject_calls(subjects, s_idx, s_frame, s_frames_preliminary)
            }
            None => None,
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-498 and
        // c++/src/algo/blast/core/aa_ungapped.c:516-547
        // A chunk whose scan reads letters before it (a range of masks stored out of order,
        // `common::translated_chunks::scan_lead`) gives hits at offsets below 0, whose tests read
        // the cells earlier chunks left within reach. For a lead subject (`x_lead_subject`), a
        // later lead frame whose earlier chunks were skipped as `SUBJECT_SPLIT_NO_RANGE` can also
        // reach cells left by an earlier non-lead frame, so the whole lead subject runs its chunks
        // in order on the worker's single diagonal table (no parallel chunks with private tables).
        let x_chunks_in_order = x_lead_subject
            || x_hard_start
            || x_soft_calls.as_ref().is_some_and(|calls| {
                calls.iter().any(|(_, ranges)| {
                    crate::algorithm::common::translated_chunks::scan_lead(ranges) > 0
                })
            });

        // What the word finder reads of this subject and frame (`word_finder::FrameCx`).
        let frame_cx = word_finder::FrameCx {
            cx,
            lookup_ref,
            s_idx,
            s_len,
            s_frames_listed,
            is_long_sequence,
            collect_hsp_saving_stats,
            s_f_idx,
            s_frame,
            subject_full,
            subject_all,
            x_soft_calls: x_soft_calls.as_deref(),
            x_hard_start,
        };

        // The word finder of a chunk (`word_finder::scan_chunk`), which the chunk driver calls.
        let scan_chunk = |chunk: SubjectChunk,
                          offset_pairs: &mut [OffsetPair],
                          diag_array: &mut [DiagStruct],
                          diag_offset: &mut i32,
                          scan_chunk_size: Option<usize>,
                          sc: S|
         -> TblastxChunkScanResult {
            word_finder::scan_chunk(
                &frame_cx,
                chunk,
                offset_pairs,
                diag_array,
                diag_offset,
                scan_chunk_size,
                sc,
            )
        };

        // The chunks of this frame, each scanned by `scan_chunk` and merged into the frame's
        // list (`search_frame_chunks`).
        search_frame_chunks(
            &ChunkDrive {
                use_serial_scan_chunks,
                use_parallel_chunks: use_parallel_chunks && !x_chunks_in_order,
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
            },
            &mut split_state,
            offset_pairs,
            diag_array,
            diag_offset,
            ChunkSums {
                frame_ungapped_hits: &mut frame_ungapped_hits,
                stats_hsp_saved,
                stats_hsp_filtered_by_cutoff,
                stats_score_distribution,
                one_hit: (window <= 0).then_some(OneHitChunkMerge {
                    failure: x_read_failure,
                    s_idx,
                    window,
                }),
                hard_merge: match subjects {
                    Subjects::Stream(stream) if x_hard_start => Some(HardChunkMerge {
                        failure: x_read_failure,
                        s_idx,
                        oid: stream.seqdb().oid(s_idx as u32),
                    }),
                    _ => None,
                },
                subject_besthit,
            },
            &scan_chunk,
            sc,
        );

        // The frame's list joins the subject's (`append_frame_hits`).
        append_frame_hits(combined_ungapped_hits, frame_ungapped_hits);
    } // End of subject frame loop
}

/// How the chunks of a subject frame are searched, made per frame by the frame loop under the
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

/// Refusal target for a hard-masked subject whose chunk ends the last hard range (NCBI reads the
/// merge overlap past its list) and would join HSPs across the split (`add_chunk_result`).
#[derive(Clone, Copy)]
pub(super) struct HardChunkMerge<'s> {
    pub(super) failure: &'s crate::algorithm::tblastx::subjects::XReadFailure,
    pub(super) s_idx: usize,
    pub(super) oid: crate::seqdb::Oid,
}

/// What each chunk's result is added to: the frame's HSP list and the subject's HSP-saving
/// statistics (`SubjectSums`), and for the one-hit word finder (`-window_size 0` or less,
/// `one_hit`) where a chunk merge NCBI leaves undefined is refused.
pub(super) struct ChunkSums<'s> {
    pub(super) frame_ungapped_hits: &'s mut Vec<UngappedHit>,
    pub(super) stats_hsp_saved: &'s mut usize,
    pub(super) stats_hsp_filtered_by_cutoff: &'s mut usize,
    pub(super) stats_score_distribution: &'s mut Vec<i32>,
    pub(super) one_hit: Option<OneHitChunkMerge<'s>>,
    /// Present for a `-db_hard_mask` subject on the overlapping-mask walk (`x_hard_start`):
    /// refuses joins when NCBI reads the merge overlap past its list of hard ranges.
    pub(super) hard_merge: Option<HardChunkMerge<'s>>,
    /// `-subject_besthit`: the frame's list is filtered after each merge.
    pub(super) subject_besthit: bool,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
    // Repeatedly invoke `scansub` while the active scan range remains non-empty.
    if use_serial_scan_chunks {
        loop {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:452-584
            let chunk = match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                SubjectChunkStatus::Done => break,
                SubjectChunkStatus::Ok(chunk) => chunk,
            };

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:83-123
            //
            // Experimental scan chunks split only the scan walk inside this
            // NCBI search unit. The reducer still extends against `subject`,
            // the full real chunk, and `Blast_ExtendWordExit` runs once, at the end
            // of the chunk's scan.
            let scan_size = tblastx_scan_chunk_size_for_run(chunk.length, num_threads);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
            // Repeatedly invoke `scansub` while the active scan range remains non-empty.
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:452-584
            let mut chunks = Vec::new();
            loop {
                match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                    SubjectChunkStatus::Done => break,
                    SubjectChunkStatus::Ok(chunk) => chunks.push(chunk),
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:452-584
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
                    // same chunk order as NCBI before the next subject frame starts
                    // (skipping `SUBJECT_SPLIT_NO_RANGE` chunks that had no exit).
                    // The experimental workers above use local diagonal tables; the
                    // ordered merge below remains the NCBI `Blast_HSPListsMerge` path.
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-483,561-584, c++/src/algo/blast/core/blast_extend.c:167-173.
                    if result.exited {
                        advance_tblastx_diag_offset(
                            diag_offset,
                            diag_array,
                            window,
                            result.chunk.ncbi_length,
                        );
                    }
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:452-584
            let chunk = match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                SubjectChunkStatus::Done => break,
                SubjectChunkStatus::Ok(chunk) => chunk,
            };

            // Reset the initial hit list for each subject chunk (NCBI `BlastInitHitListReset`).
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
            let result = scan_chunk(chunk, offset_pairs, diag_array, diag_offset, None, sc);
            add_chunk_result(&mut sums, result);
        }
    }
}

/// Adds the result of one chunk to `sums`: its HSP-saving statistics, and its HSPs merged into
/// the frame's list where they meet the list's HSPs in the chunks' overlap.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:581-584
fn add_chunk_result(sums: &mut ChunkSums<'_>, result: TblastxChunkScanResult) {
    *sums.stats_hsp_saved += result.stats.hsp_saved;
    *sums.stats_hsp_filtered_by_cutoff += result.stats.hsp_filtered_by_cutoff;
    sums.stats_score_distribution
        .extend(result.stats.score_distribution);
    if !result.hits.is_empty() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:576-584
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2873-2911
        // On the last hard range NCBI reads the merge overlap past its list of ranges; refuse if
        // the list and incoming HSPs would join across that undefined overlap.
        if let Some(merge) = sums.hard_merge {
            if result.chunk.overlap_past_list
                && !sums.frame_ungapped_hits.is_empty()
                && crate::algorithm::common::hard_chunk_scan::merge_joins(
                    &crate::core::subject_split::HardChunk {
                        offset: result.chunk.offset,
                        length: result.chunk.ncbi_length,
                        residual: 0,
                        full_range: false,
                        range: 0,
                        overlap: result.chunk.overlap,
                        overlap_past_list: true,
                    },
                    sums.frame_ungapped_hits.iter().map(|h| h.s_aa_end as i32),
                    result.hits.iter().map(|h| h.s_aa_start as i32),
                )
            {
                merge.failure.fail(
                    merge.s_idx,
                    crate::algorithm::common::hard_chunk_scan::short_chunk_merge_refusal(
                        merge.oid, "TBLASTX",
                    ),
                );
                return;
            }
        }
        if let Err(undefined) = merge_tblastx_subject_chunk_hits(
            sums.frame_ungapped_hits,
            result.hits,
            result.chunk.offset,
            result.chunk.overlap,
        ) {
            // Only an HSP of the one-hit word finder of length 0 or less makes a merge undefined.
            sums.one_hit
                .expect(TWO_HIT_CHUNK_MERGE_IS_DEFINED)
                .refuse(undefined);
        }
        if sums.subject_besthit {
            super::frame_subject_besthit(sums.frame_ungapped_hits);
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:804-845
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766
        combined_ungapped_hits.extend(frame_ungapped_hits);
        if !ungapped_hits_is_sorted_by_score_ncbi(combined_ungapped_hits) {
            sort_ungapped_hits_by_score_ncbi(combined_ungapped_hits);
        }
    }
}
