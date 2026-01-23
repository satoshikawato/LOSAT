//! BLASTN search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c
//!            ncbi-blast/c++/src/algo/blast/core/blast_engine.c
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process:
//! - Query and subject sequence reading
//! - Lookup table construction
//! - Seed finding with two-hit extension
//! - HSP filtering and output generation
//!
//! # Module Structure
//!
//! - `run` - Main run() function for BLASTN search

mod run;

// Re-export the main run function
pub use run::run;

use rustc_hash::FxHashMap;
use crate::common::Hit;
use crate::stats::KarlinParams;
use super::alignment::build_blastna_matrix;
use super::filtering::{purge_hsps_with_common_endpoints_ex, reevaluate_hsp_with_ambiguities_gapped, subject_best_hit};
use super::hsp::{BlastnHsp, score_compare_hsps as score_compare_blastn_hsps};
use super::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};

// Re-export for backward compatibility
pub use crate::algorithm::common::evalue::calculate_evalue_database_search as calculate_evalue;

/// Filter HSPs to remove redundant overlapping hits.
/// This function applies overlap filtering to match NCBI BLAST's behavior of outputting individual HSPs.
/// Chaining/clustering functionality has been removed to match NCBI BLAST blastn behavior.
///
/// NCBI reference: blast_traceback.c:637-669 (2-phase endpoint purge with trimming)
/// Phase 1: purge=FALSE - trim overlapping HSPs using edit scripts
/// Phase 2: Re-evaluate trimmed HSPs with Blast_HSPReevaluateWithAmbiguitiesGapped
/// Phase 3: purge=TRUE - delete remaining duplicates
/// Phase 4: score re-sort and interval tree containment purge
pub fn filter_hsps(
    hits: Vec<Hit>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:651-659
    // ```c
    // query = query_blk->sequence +
    //     query_info->contexts[hsp->context].query_offset;
    // query_length = query_info->contexts[hsp->context].query_length;
    // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query,
    //              query_length, subject, subject_length, hit_params,
    //              score_params, sbp);
    // ```
    encoded_queries_blastna: &[Vec<u8>],
    encoded_subjects_blastna: &[Vec<u8>],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    _db_len_total: usize,
    _db_num_seqs: usize,
    _params: &KarlinParams,
    _use_dp: bool,
    verbose: bool,
    query_lengths: &[usize], // Query index -> length mapping
    cutoff_scores: &[i32],   // Query index -> cutoff score (NCBI: hit_params->cutoff_score_min)
) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    let total_start = std::time::Instant::now();

    // Always output individual HSPs (NCBI BLAST behavior)
    // NCBI BLAST's blastn does not use chaining/clustering - all HSPs are saved individually
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
    // } BlastHSPList;
    // ```
    let mut result_hits: Vec<BlastnHsp> = hits.into_iter().map(BlastnHsp::from_hit).collect();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1052-1127 (BlastScoreBlkNuclMatrixCreate)
    let score_matrix = build_blastna_matrix(reward, penalty);

    // NCBI blastn does NOT apply culling by default (culling_limit = 0)
    // Culling is only applied when -culling_limit is explicitly set
    // Reference: blast_options.c:1496 - culling_limit defaults to 0
    // For NCBI parity, we skip the domination filtering entirely
    let filter_start = std::time::Instant::now();
    let result_hits_count = result_hits.len();

    // Sort by raw_score (descending) to match NCBI BLAST's sorting order
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1354
    // ```c
    // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
    //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
    //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
    //     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
    //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    // }
    // ```
    result_hits.sort_by(score_compare_blastn_hsps);

    // NOTE: NCBI does NOT have explicit dedup for identical coords+score.
    // NCBI only uses Blast_HSPListPurgeHSPsWithCommonEndpoints which handles
    // HSPs with common start OR end points (not both). Any exact duplicates
    // are handled by the endpoint purge below.
    // Reference: blast_hits.c:2455-2535

    // NCBI BLAST 2-phase endpoint purge with traceback trimming
    // NCBI reference: blast_traceback.c:637-669
    //
    // Int4 extra_start = Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    // for (index=extra_start; index < hsp_list->hspcnt; index++) {
    //     hsp = hsp_array[index];
    //     if (!hsp) continue;
    //     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, ...);
    //     if (delete_hsp) hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    // }
    // Blast_HSPListPurgeNullHSPs(hsp_list);
    // if(program_number == eBlastTypeBlastn) {
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // }

    let before_purge = result_hits.len();

    // Phase 1: purge=FALSE - trim overlapping HSPs using edit scripts
    // NCBI reference: blast_traceback.c:638
    let (mut result_hits, extra_start) = purge_hsps_with_common_endpoints_ex(result_hits, false);
    if verbose {
        eprintln!("[INFO]   Endpoint purge phase 1 (trim): {} -> {} hits, extra_start={}",
                  before_purge, result_hits.len(), extra_start);
    }

    // Phase 2: Re-evaluate trimmed HSPs
    // NCBI reference: blast_traceback.c:647-665
    // for (index=extra_start; index < hsp_list->hspcnt; index++) {
    //     ...
    //     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query, query_length,
    //                      subject, subject_length, hit_params, score_params, sbp);
    //     ...
    // }
    let mut deleted_count = 0;
    for i in extra_start..result_hits.len() {
        let hit = &mut result_hits[i];

        // Get cutoff score for this query
        // NCBI reference: blast_traceback.c:654 - uses hit_params->cutoff_score_min
        // The cutoff_score is pre-computed in run() using compute_blastn_cutoff_score()
        let q_idx = hit.q_idx as usize;
        let s_idx = hit.s_idx as usize;
        let cutoff_score = cutoff_scores.get(q_idx).copied().unwrap_or(22);

        // Get sequences for re-evaluation
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:651-659
        // ```c
        // query = query_blk->sequence +
        //     query_info->contexts[hsp->context].query_offset;
        // query_length = query_info->contexts[hsp->context].query_length;
        // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query,
        //              query_length, subject, subject_length, hit_params,
        //              score_params, sbp);
        // ```
        let context_idx = q_idx * 2 + if hit.query_frame < 0 { 1 } else { 0 };
        let q_blastna = match encoded_queries_blastna.get(context_idx) {
            Some(seq) => seq.as_slice(),
            None => continue,
        };
        let s_blastna = match encoded_subjects_blastna.get(s_idx) {
            Some(seq) => seq.as_slice(),
            None => continue,
        };
        let delete = reevaluate_hsp_with_ambiguities_gapped(
            hit,
            q_blastna,
            s_blastna,
            reward,
            penalty,
            gap_open,
            gap_extend,
            cutoff_score,
            &score_matrix,
        );
        if delete {
            hit.raw_score = i32::MIN; // Mark for deletion
            deleted_count += 1;
        }
    }

    // Remove marked HSPs (equivalent to Blast_HSPListPurgeNullHSPs)
    // NCBI reference: blast_traceback.c:666
    result_hits.retain(|h| h.raw_score != i32::MIN);
    if verbose && deleted_count > 0 {
        eprintln!("[INFO]   Re-evaluation deleted {} trimmed HSPs", deleted_count);
    }

    // Phase 3: purge=TRUE for BLASTN
    // NCBI reference: blast_traceback.c:668
    // if(program_number == eBlastTypeBlastn) {
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // }
    let before_phase3 = result_hits.len();
    let (mut result_hits, _) = purge_hsps_with_common_endpoints_ex(result_hits, true);
    if verbose && before_phase3 != result_hits.len() {
        eprintln!("[INFO]   Endpoint purge phase 3 (delete): {} -> {} hits", before_phase3, result_hits.len());
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:671-688
    // ```c
    // Blast_HSPListSortByScore(hsp_list);
    // Blast_IntervalTreeReset(tree);
    // for (index = 0; index < hsp_list->hspcnt; index++) {
    //     if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
    //                                      hit_options->min_diag_separation)) {
    //         hsp_array[index] = Blast_HSPFree(hsp);
    //     } else {
    //         BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
    //     }
    // }
    // ```
    result_hits.sort_by(score_compare_blastn_hsps);
    let mut final_hits: Vec<BlastnHsp> = Vec::with_capacity(result_hits.len());
    let min_diag_separation: i32 = 0;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1482-1495
    // ```c
    // if (min_diag_separation)
    //     options->min_diag_separation = min_diag_separation;
    // ```
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
    // } BlastHSPList;
    // ```
    let mut trees: FxHashMap<(u32, u32), BlastIntervalTree> = FxHashMap::default();
    for hit in result_hits {
        let q_idx = hit.q_idx as usize;
        let s_idx = hit.s_idx as usize;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:651-655
        // ```c
        // query = query_blk->sequence +
        //     query_info->contexts[hsp->context].query_offset;
        // query_length = query_info->contexts[hsp->context].query_length;
        // ```
        let context_idx = q_idx * 2 + if hit.query_frame < 0 { 1 } else { 0 };
        let q_seq = match encoded_queries_blastna.get(context_idx) {
            Some(seq) => seq,
            None => {
                final_hits.push(hit);
                continue;
            }
        };
        let s_seq = match encoded_subjects_blastna.get(s_idx) {
            Some(seq) => seq,
            None => {
                final_hits.push(hit);
                continue;
            }
        };
        let query_len = query_lengths
            .get(q_idx)
            .copied()
            .unwrap_or_else(|| hit.query_length.max(q_seq.len()));
        let subject_len = s_seq.len();
        if query_len == 0 || subject_len == 0 {
            final_hits.push(hit);
            continue;
        }
        let key = (hit.q_idx, hit.s_idx);
        let tree = trees.entry(key).or_insert_with(|| {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:537-550
            // ```c
            // query_start = s_GetQueryStrandOffset(query_info, hsp->context);
            // region_start = query_start + hsp->query.offset;
            // region_end = query_start + hsp->query.end;
            // ```
            BlastIntervalTree::new(
                0,
                (2 * query_len + 1) as i32,
                0,
                (subject_len + 1) as i32,
            )
        });
        // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
        // ```c
        // query_info->contexts[0].query_offset = 0;
        // query_info->contexts[1].query_offset = kStrandLength + 1;
        // ```
        let query_context_offset = if hit.query_frame < 0 {
            (query_len + 1) as i32
        } else {
            0
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
        // ```c
        // if (hsp->query.frame != hsp->subject.frame) {
        //    *q_end = query_length - hsp->query.offset;
        //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
        // }
        // ```
        let (q_offset_0, q_end_0) = if hit.query_frame < 0 {
            (
                query_len.saturating_sub(hit.q_end),
                query_len
                    .saturating_sub(hit.q_start)
                    .saturating_add(1),
            )
        } else {
            (hit.q_start.saturating_sub(1), hit.q_end)
        };
        let tree_hsp = TreeHsp {
            query_offset: q_offset_0 as i32,
            query_end: q_end_0 as i32,
            subject_offset: hit.s_start.min(hit.s_end).saturating_sub(1) as i32,
            subject_end: hit.s_start.max(hit.s_end) as i32,
            score: hit.raw_score,
            query_frame: hit.query_frame,
            query_length: query_len as i32,
            query_context_offset,
            subject_frame_sign: 1,
        };
        if tree.contains_hsp(&tree_hsp, query_context_offset, min_diag_separation) {
            continue;
        }
        tree.add_hsp(tree_hsp, query_context_offset, IndexMethod::QueryAndSubject);
        final_hits.push(hit);
    }

    // NCBI BLAST: Subject Best Hit filtering (OPTIONAL - not applied by default)
    // Reference: blast_hits.c:2537-2606 Blast_HSPListSubjectBestHit
    // Called from blast_engine.c:586-591 ONLY when subject_besthit_opts is set
    // (i.e., when user specifies -subject_besthit on command line)
    //
    // NOTE: Subject best hit filtering is disabled by default to match NCBI behavior.
    // The subject_best_hit() function is available but not called here.
    // To enable, add a command-line option and call subject_best_hit() when enabled.
    let _ = subject_best_hit; // Silence unused import warning

    let filter_time = filter_start.elapsed();

    if verbose {
        eprintln!("[INFO] filter_hsps summary:");
        eprintln!(
            "[INFO]   Total time: {:.3}s",
            total_start.elapsed().as_secs_f64()
        );
        eprintln!(
            "[INFO]   Filtering: {:.3}s ({} -> {} hits)",
            filter_time.as_secs_f64(),
            result_hits_count,
            final_hits.len()
        );
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
    let mut output_hits: Vec<Hit> = Vec::with_capacity(final_hits.len());
    for hit in final_hits {
        output_hits.push(Hit {
            identity: hit.identity,
            length: hit.length,
            mismatch: hit.mismatch,
            gapopen: hit.gapopen,
            q_start: hit.q_start,
            q_end: hit.q_end,
            s_start: hit.s_start,
            s_end: hit.s_end,
            e_value: hit.e_value,
            bit_score: hit.bit_score,
            query_frame: hit.query_frame,
            query_length: hit.query_length,
            q_idx: hit.q_idx,
            s_idx: hit.s_idx,
            raw_score: hit.raw_score,
            gap_info: hit.gap_info,
        });
    }
    output_hits
}
